#!/usr/bin/env python3

import base64
import json
import os
from html import escape
from importlib.metadata import version
from io import BytesIO
from math import ceil, isfinite
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"

import dask
import numpy as np
import shapely
from lazyslide import open_wsi, tl
from matplotlib.figure import Figure
from spatialdata import SpatialData, read_zarr, to_polygons, transform
from spatialdata.models import ShapesModel
from xarray import DataTree


def make_regions(slide, tiles, mask):
    wsi = open_wsi(slide, image_key="wsi")
    wsi.shapes["tiles"] = transform(tiles, to_coordinate_system="global")
    tl.tile_shaper(wsi, groupby="domain", key_added="tissue_regions")

    label = mask.labels["wsi_tissue"]
    if isinstance(label, DataTree):
        label = label[f"scale{len(label.children) - 1}"].image
    contours = to_polygons(label)
    # Use the same pixel-edge convention and tissue contours as the tiling step.
    contours.geometry = contours.geometry.translate(xoff=0.5, yoff=0.5)
    contours = transform(contours, to_coordinate_system="global")
    height, width = wsi.properties.shape
    foreground = contours.geometry.union_all().intersection(
        shapely.box(0, 0, width, height)
    )

    regions = wsi.shapes["tissue_regions"].clip(foreground).explode(ignore_index=True)
    regions = regions.loc[regions.geom_type.eq("Polygon") & regions.area.gt(0)].copy()
    bounds = regions.bounds
    order = (
        bounds.assign(domain=regions["domain"].astype(str))
        .sort_values(["domain", "miny", "minx", "maxy", "maxx"], kind="stable")
        .index
    )
    regions = regions.loc[order].reset_index(drop=True)
    regions["region_id"] = [f"region_{i}" for i in range(1, len(regions) + 1)]
    return regions


def area_without_overlap(regions):
    # Overlap is ambiguous for assignment, so exclude it from each density denominator.
    geometry = regions.geometry.to_numpy()
    shared_area = np.zeros(len(regions))
    index = regions.sindex
    for i, polygon in enumerate(geometry):
        neighbors = index.query(polygon, predicate="intersects")
        others = neighbors[neighbors != i]
        if len(others):
            shared = shapely.intersection(polygon, geometry[others])
            shared_area[i] = shapely.area(shapely.union_all(shared))
    area = shapely.area(geometry)
    return np.maximum(area - shared_area, 0), shared_area


def overlay(slide, regions):
    image = slide.images["wsi"]
    levels = (
        [level.image for level in image.children.values()]
        if isinstance(image, DataTree)
        else [image]
    )
    base = levels[0]
    width, height = base.sizes["x"], base.sizes["y"]
    preview = next(
        (level for level in levels if max(level.sizes["x"], level.sizes["y"]) <= 1024),
        levels[-1],
    )
    step = ceil(max(preview.sizes["x"], preview.sizes["y"]) / 1024)
    pixels = preview.isel(y=slice(None, None, step), x=slice(None, None, step))
    figure = Figure(figsize=(7, 6), layout="constrained")
    ax = figure.subplots()
    ax.imshow(pixels.transpose("y", "x", "c").values, extent=(0, width, height, 0))
    regions.plot(
        ax=ax, column="domain", categorical=True, cmap="tab20", alpha=0.45, legend=True
    )
    ax.set(xlim=(0, width), ylim=(height, 0))
    ax.set_axis_off()
    stream = BytesIO()
    figure.savefig(
        stream, format="jpeg", dpi=120, pil_kwargs={"quality": 85, "optimize": True}
    )
    image = base64.b64encode(stream.getvalue()).decode("ascii")
    return (
        '<div class="mqc-custom-content-image">'
        f'<img src="data:image/jpeg;base64,{image}" alt="Tissue regions on the original slide" /></div>'
    )


tile_model = "${tile_model}"
sample = "${meta.id}"
domain_store = read_zarr("domains.zarr")
domains = domain_store.attrs["tissue_domains"]
if domains["tile_model"] != tile_model or domains["sample_id"] != sample:
    raise ValueError(
        "The tissue domain map does not match the selected slide and tile encoder."
    )

slide = read_zarr("slide.zarr")
mpp = float(slide.attrs["ingestion"]["mpp"])
if not isfinite(mpp) or mpp <= 0:
    raise ValueError("Tissue region measurements require a positive slide pixel size.")

mask = read_zarr("tissue.zarr")
mask_info = mask.attrs["tissue_detection"]
tiles = domain_store.shapes.get("tiles")
shapes = {}
html = ""
summary = "No tissue regions available."
with dask.config.set(scheduler="threads", num_workers=int("${task.cpus}")):
    if tiles is not None and len(tiles) and mask_info["has_tissue"]:
        regions = make_regions(slide, tiles, mask)
        if len(regions):
            area, shared_area = area_without_overlap(regions)
            regions["tissue_area_mm2"] = area * (mpp / 1000) ** 2
            regions["overlap_area_mm2"] = shared_area * (mpp / 1000) ** 2
            shapes["tissue_regions"] = ShapesModel.parse(regions)
            html = overlay(slide, regions)
            area_mm2 = regions["tissue_area_mm2"].sum()
            summary = f"{len(regions):,} regions from {regions['domain'].nunique():,} domains. Tissue area available for assignment is {area_mm2:.3f} mm²."

output = "${prefix}.${tile_model}.tissue_regions"
metadata = {
    "sample_id": sample,
    "tile_model": tile_model,
    "domain_resolution": domains["resolution"],
    "tissue_method": mask_info["method"],
    "mpp": mpp,
    "region_count": len(shapes.get("tissue_regions", ())),
}
SpatialData(shapes=shapes, attrs={"tissue_regions": metadata}).write(f"{output}.zarr")
payload = {
    "id": output,
    "parent_id": "histopathology_tissue_regions",
    "parent_name": "Tissue regions",
    "parent_description": "Connected regions of similar tile morphology within the primary tissue mask.",
    "section_name": f"{sample} · {tile_model}",
    "plot_type": "image" if html else "html",
    "data": f"{html}<p>{escape(summary)}</p>",
}
Path(f"{output}_mqc.json").write_text(json.dumps(payload), encoding="utf-8")
packages = ("geopandas", "lazyslide", "matplotlib", "shapely", "spatialdata", "wsidata")
versions = {package: version(package) for package in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
