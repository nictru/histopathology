#!/usr/bin/env python3

import base64
import json
import os
from html import escape
from importlib.metadata import version
from io import BytesIO
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
os.environ["NUMBA_NUM_THREADS"] = "${task.cpus}"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
from lazyslide import open_wsi, pl, pp, tl
from scipy.sparse import csr_matrix
from sklearn.preprocessing import normalize
from spatialdata import SpatialData, read_zarr
from spatialdata._io import write_shapes
from threadpoolctl import threadpool_limits
from zarr.storage import FsspecStore


def smooth_features(wsi, features):
    tiles = wsi.shapes["tiles"]
    spec = wsi.tile_spec("tiles")
    pp.tile_graph(wsi, n_neighs=min(8, len(tiles) - 1))
    graph = wsi.tables["tiles_graph"].obsp["spatial_connectivities"].tocoo()
    coords = tiles.bounds[["minx", "miny"]].to_numpy()
    dx, dy = np.abs(coords[graph.row] - coords[graph.col]).T
    tissue = tiles["tissue_id"].to_numpy()
    # Limit smoothing to adjacent grid positions within each tissue piece.
    local = (
        (dx <= spec.base_stride_width + 1)
        & (dy <= spec.base_stride_height + 1)
        & (tissue[graph.row] == tissue[graph.col])
    )
    adjacency = csr_matrix(
        (graph.data[local], (graph.row[local], graph.col[local])),
        shape=graph.shape,
        dtype=np.result_type(features.dtype, np.float32),
    )
    adjacency = adjacency.maximum(adjacency.T)
    adjacency.setdiag(1.0)
    # Keep the self and neighbor averages sparse for large tile grids.
    return normalize(adjacency, norm="l1") @ features


def assign_domains(wsi, feature_key, resolution):
    table = wsi.tables[feature_key]
    tiles = wsi.shapes["tiles"]
    if len(tiles) < 3 or table.n_vars < 2:
        tiles["domain"] = pd.Categorical(["0"] * len(tiles))
        return (
            "Too few tiles or feature dimensions for clustering. Assigned one domain."
        )

    smoothed = smooth_features(wsi, table.X)
    if np.allclose(smoothed, smoothed[0]):
        tiles["domain"] = pd.Categorical(["0"] * len(tiles))
        return "Smoothed tile embeddings are constant. Assigned one domain."

    table.layers["spatial_features"] = smoothed
    tl.spatial_domain(
        wsi, feature_key, layer_key="spatial_features", resolution=resolution
    )
    labels = wsi.shapes["tiles"]["domain"]
    wsi.shapes["tiles"]["domain"] = pd.Categorical(labels.astype(str))
    return ""


def overlay(wsi):
    figure, ax = plt.subplots(figsize=(8, 6), layout="constrained")
    pl.tiles(
        wsi,
        color="domain",
        ax=ax,
        alpha=0.4,
        in_bounds=False,
        show_contours=False,
        mark_origin=False,
        scalebar=False,
        linewidth=0,
        target_dpi=135,
    )
    ax.set_title("")
    ax.set_axis_off()
    output = BytesIO()
    figure.savefig(
        output, format="jpeg", dpi=135, pil_kwargs={"quality": 88, "optimize": True}
    )
    plt.close(figure)
    image = base64.b64encode(output.getvalue()).decode("ascii")
    return (
        '<div class="mqc-custom-content-image">'
        f'<img src="data:image/jpeg;base64,{image}" alt="Tissue domains overlaid on the original slide" /></div>'
    )


sample = "${meta.id}"
tile_model = "${tile_model}"
prefix = "${prefix}"
resolution = float("${domain_resolution}")
store = read_zarr("embeddings.zarr")
tile_info = store.attrs["tile_embeddings"]
tile_count = tile_info["tile_count"]
if tile_info["model"] != tile_model:
    raise ValueError(
        f"Expected {tile_model} tile embeddings, found {tile_info['model']}"
    )

shapes = {}
attrs = {}
domain_count = 0
image = ""
note = "No eligible tiles."
if tile_count:
    feature_key = f"{tile_model}_tiles"
    table = store.tables[feature_key]
    tiles = store.shapes["tiles"].set_index("tile_id", drop=False)
    ids = pd.Index(table.obs["tile_id"])
    if (
        len(tiles) != tile_count
        or len(ids) != tile_count
        or not tiles.index.is_unique
        or not ids.is_unique
        or not ids.isin(tiles.index).all()
    ):
        raise ValueError(
            f"Tile IDs do not identify {tile_count} unique tile embeddings"
        )

    # LazySlide joins the geometries and feature rows by position.
    wsi = open_wsi(read_zarr("source.zarr"), image_key="wsi")
    wsi.shapes["tiles"] = tiles.loc[ids].reset_index(drop=True)
    wsi.tables[feature_key] = table
    wsi.attrs["tile_spec"] = store.attrs["tile_spec"]
    with threadpool_limits(limits=int("${task.cpus}")):
        note = assign_domains(wsi, feature_key, resolution)
        image = overlay(wsi)
    shapes["tiles"] = wsi.shapes["tiles"]
    domain_count = int(shapes["tiles"]["domain"].nunique())
    attrs["tile_spec"] = store.attrs["tile_spec"]

attrs["tissue_domains"] = {
    "sample_id": sample,
    "tile_model": tile_model,
    "resolution": resolution,
    "tile_count": tile_count,
    "domain_count": domain_count,
    "note": note,
}
output = f"{prefix}.{tile_model}.tissue_domains"
spatial_data = SpatialData(shapes=shapes, attrs=attrs)
output_path = Path(f"{output}.zarr")
zarr_store = FsspecStore.from_url(str(output_path.resolve()), storage_options={"auto_mkdir": True})
root = zarr.create_group(store=zarr_store, zarr_format=3)
spatial_data.write_attrs(zarr_group=root)
for shape_name, output_shapes in spatial_data.shapes.items():
    shape_path = f"shapes/{shape_name}"
    root.require_group(shape_path)
    # SpatialData requires a local group handle for the GeoParquet sidecar.
    write_shapes(output_shapes, zarr.open_group(output_path / shape_path, mode="r+"))
zarr.consolidate_metadata(zarr_store)
domain_label = "domain" if domain_count == 1 else "domains"
tile_label = "tile" if tile_count == 1 else "tiles"
summary = f"{domain_count:,} {domain_label} across {tile_count:,} {tile_label}. {note}".strip()
payload = {
    "id": output,
    "parent_id": "histopathology_tissue_domains",
    "parent_name": "Tissue domains",
    "parent_description": (
        "Clusters of spatially smoothed tile embeddings overlaid on the original slides. "
        "Domain labels are assigned independently for each slide and encoder."
    ),
    "section_name": f"{sample} · {tile_model}",
    "plot_type": "image" if image else "html",
    "data": f"{image}<p>{escape(summary)}</p>",
}
Path(f"{output}_mqc.json").write_text(json.dumps(payload), encoding="utf-8")
packages = ("igraph", "lazyslide", "matplotlib", "scanpy", "spatialdata", "wsidata")
versions = {package: version(package) for package in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
