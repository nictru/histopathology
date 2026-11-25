#!/usr/bin/env python3

import base64
import json
import os
from html import escape
from importlib.metadata import version
from io import BytesIO
from pathlib import Path

os.environ["NUMBA_NUM_THREADS"] = "${task.cpus}"
os.environ["NUMBA_CACHE_DIR"] = str(Path(".numba-cache").resolve())
os.environ["MPLBACKEND"] = "Agg"

import dask
import matplotlib
import numpy as np
import pandas as pd
from anndata import AnnData
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from PIL import Image
from spatialdata import SpatialData, read_zarr, transform
from spatialdata.models import TableModel
from spatialdata.transformations import Identity, get_transformation
from threadpoolctl import threadpool_limits
from xarray import DataTree

CPUS = int("${task.cpus}")
CELL_TYPE = "cell_type"
SECTION = "histopathology_cell_niches"


def classify_neighborhoods(adata, radius, resolution):
    from squidpy import gr

    labels = pd.Series(pd.NA, index=adata.obs_names, dtype="string")
    degrees = np.zeros(adata.n_obs, dtype=int)
    classes = adata.obs[CELL_TYPE].cat.categories
    profiles = np.zeros((adata.n_obs, len(classes)), dtype=np.float32)
    if adata.n_obs < 2:
        return (
            labels,
            degrees,
            profiles,
            "Too few classified instances for a neighborhood.",
        )

    gr.spatial_neighbors_radius(adata, radius=radius, n_jobs=CPUS)
    graph = adata.obsp["spatial_connectivities"]
    degrees = graph.getnnz(axis=1)
    connected = degrees > 0
    n_connected = int(connected.sum())
    if n_connected == 0:
        return (
            labels,
            degrees,
            profiles,
            "No classified neighbors within the selected radius.",
        )

    if len(classes) == 1:
        profiles[connected, 0] = 1
    elif n_connected < 3:
        members = pd.get_dummies(adata.obs[CELL_TYPE], dtype=np.float32).to_numpy()
        profiles[connected] = (graph @ members)[connected] / degrees[connected, None]
        return (
            labels,
            degrees,
            profiles,
            "Too few instances with neighbors for clustering.",
        )
    else:
        gr.calculate_niche_neighborhood(
            adata,
            groups=CELL_TYPE,
            resolutions=resolution,
            scale=False,
            n_neighbors=min(15, n_connected - 1),
            cluster_mask=pd.Series(connected, index=adata.obs_names),
            rng=0,
        )
        profiles = adata.obsm.pop("niche_embedding").astype(np.float32, copy=False)
        labels = adata.obs[f"nhood_niche_res={resolution}"].astype("string")
        labels = labels.mask(labels == "not_a_niche")

    connected_profiles = profiles[connected]
    if np.array_equal(connected_profiles.min(axis=0), connected_profiles.max(axis=0)):
        labels.loc[connected] = "0"
        return (
            labels,
            degrees,
            profiles,
            "Connected instances share one neighborhood composition.",
        )
    return labels, degrees, profiles, ""


def image_tag(jpeg, alt):
    data = base64.b64encode(jpeg.getvalue()).decode("ascii")
    return f'<img src="data:image/jpeg;base64,{data}" alt="{escape(alt, quote=True)}" style="max-width:100%;height:auto">'


def report_images(image, adata, groups, colors, crop_px):
    levels = (
        [level.image for level in image.children.values()]
        if isinstance(image, DataTree)
        else [image]
    )
    base = levels[0]
    height, width = base.sizes["y"], base.sizes["x"]
    preview = next(
        (level for level in levels if max(level.sizes["y"], level.sizes["x"]) <= 960),
        levels[-1],
    )
    step = max(1, int(np.ceil(max(preview.sizes["y"], preview.sizes["x"]) / 960)))
    pixels = preview.isel(y=slice(None, None, step), x=slice(None, None, step))
    figure = Figure(figsize=(7, 6), layout="constrained")
    ax = figure.subplots()
    ax.imshow(pixels.transpose("y", "x", "c").values, extent=(0, width, height, 0))
    xy = adata.obsm["spatial"]
    for name, positions in groups.items():
        ax.scatter(*xy[positions].T, c=[colors[name]], s=4, linewidths=0, alpha=0.65)
    ax.legend(
        handles=[Patch(facecolor=colors[name], label=name) for name in groups],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        ncols=min(5, len(groups)),
        frameon=False,
    )
    ax.set(xlim=(0, width), ylim=(height, 0))
    ax.set_axis_off()
    output = BytesIO()
    figure.savefig(
        output, format="jpeg", dpi=120, pil_kwargs={"quality": 85, "optimize": True}
    )
    overlay = image_tag(output, "Spatial cell niches overlaid on the original slide")

    crops = []
    crop_w, crop_h = min(crop_px, width), min(crop_px, height)
    for name, positions in groups.items():
        profiles = adata.X[positions]
        closest = np.argmin(np.sum((profiles - profiles.mean(axis=0)) ** 2, axis=1))
        x, y = xy[positions[closest]]
        x = int(np.clip(x - crop_w // 2, 0, width - crop_w))
        y = int(np.clip(y - crop_h // 2, 0, height - crop_h))
        crop = base.isel(y=slice(y, y + crop_h), x=slice(x, x + crop_w))
        thumbnail = Image.fromarray(crop.transpose("y", "x", "c").values).convert("RGB")
        thumbnail.thumbnail((220, 220), Image.Resampling.LANCZOS)
        output = BytesIO()
        thumbnail.save(output, format="JPEG", quality=85, optimize=True)
        crop_html = image_tag(output, f"Representative crop for {name}")
        color = matplotlib.colors.to_hex(colors[name])
        crops.append(
            f'<figure>{crop_html}<figcaption style="border-bottom:4px solid {color}">'
            f"{name} · x={x:,}, y={y:,}</figcaption></figure>"
        )
    return overlay + '<div class="d-flex flex-wrap gap-3">' + "".join(crops) + "</div>"


def write_report(prefix, title, note, image, adata, radius, mpp):
    common = {
        "parent_id": SECTION,
        "parent_name": "Spatial cell niches",
        "parent_description": (
            "Squidpy clusters instances by the cell-type composition around each instance. "
            "Niche labels are fitted independently for each slide and cell-typing method."
        ),
    }
    groups = {}
    if adata is not None:
        groups = {
            name: np.flatnonzero(adata.obs["niche"].eq(name).to_numpy())
            for name in adata.obs["niche"].cat.categories
        }
    overview = {
        **common,
        "id": f"{prefix}.cell_niches.overview",
        "section_name": f"{title} · Map and crops",
        "plot_type": "html",
        "data": f"<p>{escape(note)}</p>" if note else "",
    }
    if groups:
        assigned = adata.obs["niche"].notna().sum()
        overview["data"] += (
            f"<p>{assigned:,} of {adata.n_obs:,} classified instances assigned to a niche.</p>"
        )
        palette = (
            matplotlib.colormaps["tab20"].colors[::2]
            + matplotlib.colormaps["tab20"].colors[1::2]
        )
        if len(groups) > len(palette):
            palette = matplotlib.colormaps["turbo"](np.linspace(0, 1, len(groups)))
        colors = dict(zip(groups, palette))
        crop_px = int(np.clip(4 * radius / mpp, 64, 1024))
        overview["data"] += report_images(image, adata, groups, colors, crop_px)
        heatmap_id = f"{prefix}.cell_niches.composition"
        composition = {
            **common,
            "id": heatmap_id,
            "section_name": f"{title} · Composition",
            "description": f"Mean local cell-type proportions within {radius:g} µm. Parentheses show the number of instances assigned to each niche.",
            "plot_type": "heatmap",
            "xcats": adata.var_names.tolist(),
            "ycats": [
                f"{name} ({len(positions):,})" for name, positions in groups.items()
            ],
            "pconfig": {
                "id": f"{heatmap_id}_heatmap",
                "title": "Spatial cell niche composition",
                "xlab": "Neighbor cell type",
                "ylab": "Niche",
                "zlab": "Mean neighbor proportion (%)",
                "min": 0,
                "max": 100,
                "xcats_samples": False,
                "ycats_samples": False,
                "cluster_rows": False,
                "cluster_cols": False,
            },
            "data": [
                np.round(adata.X[positions].mean(axis=0, dtype=float) * 100, 2).tolist()
                for positions in groups.values()
            ],
        }
        Path(f"{heatmap_id}_mqc.json").write_text(
            json.dumps(composition, allow_nan=False), encoding="utf-8"
        )
    Path(f"{prefix}.cell_niches.overview_mqc.json").write_text(
        json.dumps(overview), encoding="utf-8"
    )


prefix = "${prefix}"
sample = "${meta.id}"
radius = float("${radius}")
resolution = float("${leiden_resolution}")
store = read_zarr("${instances}")
shape_key = "cells" if "cells" in store.shapes else "nuclei"
instances = store.shapes.get(shape_key)
segmenter = store.attrs["segmentation"]["model"]
typing = store.attrs.get("cell_typing", {})
classifier = typing.get("model", "segmentation")
organ = typing.get("organ")
column = "cytoformer_cell_type" if classifier == "cytoformer" else "class"
title = f"{sample} · {segmenter}"
if classifier == "cytoformer":
    title += f" + CytoFormer ({organ})"

slide = read_zarr("slide.zarr")
mpp = float(slide.attrs["ingestion"]["mpp"])
if not np.isfinite(mpp) or mpp <= 0:
    raise ValueError("Spatial cell niches require a finite positive slide pixel size.")

table = None
note = "No classified instances."
if instances is not None and column in instances:
    labels = instances[column].astype("string").str.strip()
    keep = labels.notna() & ~labels.str.casefold().isin(
        ["", "background", "unclassified"]
    )
    if keep.any():
        classified = instances.loc[keep]
        if not isinstance(
            get_transformation(classified, to_coordinate_system="global"), Identity
        ):
            classified = transform(classified, to_coordinate_system="global")
        centers = classified.geometry.centroid
        coordinates = np.column_stack((centers.x, centers.y))
        if not np.isfinite(coordinates).all():
            raise ValueError(
                "Classified instances must have finite global centroid coordinates."
            )
        obs = pd.DataFrame(
            {
                "instance_id": instances.index[keep].to_numpy(),
                "region": pd.Categorical([shape_key] * len(coordinates)),
                CELL_TYPE: pd.Categorical(labels[keep]),
            },
            index=instances.index[keep].astype(str),
        )
        adata = AnnData(obs=obs, obsm={"spatial": coordinates * mpp})
        with threadpool_limits(limits=CPUS):
            labels, degrees, profiles, note = classify_neighborhoods(
                adata, radius, resolution
            )
        niche_names = {
            value: f"Niche {int(value) + 1}"
            for value in sorted(labels.dropna().unique(), key=int)
        }
        obs["niche"] = pd.Categorical(
            labels.map(niche_names), categories=list(niche_names.values())
        )
        obs["classified_neighbors"] = degrees
        table = AnnData(
            profiles,
            obs=obs,
            var=pd.DataFrame(index=obs[CELL_TYPE].cat.categories),
            obsm={"spatial": coordinates},
        )
        TableModel.parse(
            table, region=shape_key, region_key="region", instance_key="instance_id"
        )

metadata = {
    "sample_id": sample,
    "segmentation_model": segmenter,
    "classifier": classifier,
    "organ": organ,
    "radius_um": radius,
    "resolution": resolution,
    "classified_instances": table.n_obs if table is not None else 0,
    "assigned_instances": int(table.obs["niche"].notna().sum())
    if table is not None
    else 0,
    "niche_count": len(table.obs["niche"].cat.categories) if table is not None else 0,
    "note": note,
}
SpatialData(
    shapes={shape_key: instances} if instances is not None else {},
    tables={"cell_niches": table} if table is not None else {},
    attrs={"cell_niches": metadata},
).write(f"{prefix}.cell_niches.zarr")
with dask.config.set(scheduler="threads", num_workers=CPUS):
    write_report(prefix, title, note, slide.images["wsi"], table, radius, mpp)
versions = {
    name: version(name)
    for name in ("squidpy", "igraph", "spatialdata", "anndata", "matplotlib")
}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
