#!/usr/bin/env python3

import json
from importlib.metadata import version
from pathlib import Path

import dask
from lazyslide import open_wsi, pp
from spatialdata import SpatialData, read_zarr, to_polygons, transform
from wsidata.io import add_tissues
from xarray import DataTree

slide = open_wsi(read_zarr("source.zarr"), image_key="wsi")
mask = read_zarr("tissue.zarr")
if mask.attrs["tissue_detection"]["has_tissue"]:
    label = mask.labels["wsi_tissue"]
    if isinstance(label, DataTree):
        label = label[f"scale{len(label.children) - 1}"].image
    with dask.config.set(scheduler="threads", num_workers=int("${task.cpus}")):
        contours = to_polygons(label)
    # Shift from the vectorizer's pixel centers to the tiler's pixel edges.
    contours.geometry = contours.geometry.translate(xoff=0.5, yoff=0.5)
    contours = transform(contours, to_coordinate_system="global")
    if len(contours):
        add_tissues(slide, "tissues", contours.geometry.explode(ignore_index=True))
        pp.tile_tissues(
            slide,
            int("${pixels}"),
            mpp=float("${resolution}"),
            overlap=float("${tile_overlap}"),
            background_fraction=float("${background}"),
        )

shapes = {"tiles": slide.shapes["tiles"]} if "tiles" in slide.shapes else {}
attrs = {"tile_spec": slide.attrs.get("tile_spec", {})}
SpatialData(shapes=shapes, attrs=attrs).write("${prefix}.tiles.zarr")

versions = {name: version(name) for name in ("lazyslide", "spatialdata", "wsidata")}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
