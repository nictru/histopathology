#!/usr/bin/env python3

import json
from importlib.metadata import version
from math import ceil
from pathlib import Path

import dask
import dask.array as da
import numpy as np
import zarr
from dask_image.ndinterp import affine_transform
from lazyslide import open_wsi, pp
from rasterio.features import rasterize
from spatialdata import SpatialData, read_zarr
from spatialdata._io import write_labels
from spatialdata.models import Image2DModel, Labels2DModel
from spatialdata.transformations import get_transformation
from xarray import DataTree

IMAGE = "wsi"
TISSUE = "wsi_tissue"
PREVIEW_EDGE = 2048


def detect_preview(levels):
    preview = next(
        (level for level in levels if max(level.sizes["x"], level.sizes["y"]) <= PREVIEW_EDGE),
        levels[-1],
    )
    factors = {axis: ceil(preview.sizes[axis] / PREVIEW_EDGE) for axis in ("y", "x")}
    if max(factors.values()) > 1:
        preview = preview.coarsen(factors, boundary="pad").mean().round().astype(np.uint8)
    pixels = preview.data.compute()
    shape = pixels.shape[1:]
    if np.ptp(pixels, axis=(1, 2)).any():
        preview_wsi = open_wsi(SpatialData(images={IMAGE: Image2DModel.parse(pixels)}), image_key=IMAGE)
        pp.find_tissues(
            preview_wsi, level=0, method="otsu", to_hsv=True, threshold=None,
            filter_artifacts=False, in_bounds=False,
        )
        if "tissues" in preview_wsi.shapes:
            return rasterize(preview_wsi.shapes["tissues"].geometry, out_shape=shape, dtype="uint8")
    return np.zeros(shape, dtype=np.uint8)


def labels_for_level(mask, level):
    shape = (level.sizes["y"], level.sizes["x"])
    chunks = level.data.chunksize[-2:]
    if not mask.any():
        aligned = da.zeros(shape, chunks=chunks, dtype=np.uint8)
    elif mask.shape == shape:
        aligned = da.from_array(mask, chunks=chunks)
    else:
        ratio = np.asarray(mask.shape) / shape
        aligned = affine_transform(
            da.from_array(mask), np.diag(ratio), offset=(ratio - 1) / 2,
            output_shape=shape, output_chunks=chunks, order=0, mode="nearest",
        )
    return Labels2DModel.parse(
        aligned, dims=("y", "x"), transformations=get_transformation(level, get_all=True),
    )


image = read_zarr("source.zarr").images[IMAGE]
levels = [level.image for level in image.children.values()] if isinstance(image, DataTree) else [image]
with dask.config.set(scheduler="threads", num_workers=int("${task.cpus}")):
    mask = detect_preview(levels)
    labels = labels_for_level(mask, levels[0])
    if len(levels) > 1:
        # Squidpy samples at the coarsest level and applies the mask at full resolution.
        labels = DataTree.from_dict({
            "scale0": labels.to_dataset(name="image"),
            "scale1": labels_for_level(mask, levels[-1]).to_dataset(name="image"),
        })
    attrs = {"tissue_detection": {"method": "otsu", "channel": "saturation", "has_tissue": bool(mask.any())}}
    output = "${prefix}.tissue.zarr"
    SpatialData(attrs=attrs).write(output, consolidate_metadata=False)
    root = zarr.open_group(output, mode="r+")
    write_labels(labels, root, TISSUE, storage_options={"chunks": (256, 256), "shards": (2048, 2048)})
    zarr.consolidate_metadata(root.store)

versions = {name: version(name) for name in ("lazyslide", "spatialdata", "wsidata")}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
