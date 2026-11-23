#!/usr/bin/env python3

import json
from dataclasses import asdict
from importlib.metadata import version
from math import ceil
from pathlib import Path

import dask
import dask.array as da
import numpy as np
from dask_image.ndinterp import affine_transform
from spatialdata import SpatialData, read_zarr
from spatialdata.models import Labels2DModel
from spatialdata.transformations import get_transformation
from squidpy.experimental import im
from xarray import DataTree

IMAGE = "wsi"
TISSUE = "wsi_tissue"
PREVIEW_EDGE = 2048


def tissue_preview(image):
    levels = [level.image for level in image.children.values()] if isinstance(image, DataTree) else [image]
    levels.sort(key=lambda level: level.sizes["x"] * level.sizes["y"], reverse=True)
    preview = next(
        (level for level in levels if max(level.sizes["x"], level.sizes["y"]) <= PREVIEW_EDGE),
        levels[-1],
    )
    factors = {axis: ceil(preview.sizes[axis] / PREVIEW_EDGE) for axis in ("y", "x")}
    if max(factors.values()) > 1:
        preview = preview.coarsen(factors, boundary="trim").mean(keep_attrs=True).round().astype(np.uint8)
    preview = preview.compute()
    slide = SpatialData(images={IMAGE: preview.chunk()})
    if not np.ptp(preview.data, axis=(1, 2)).any():
        return slide, np.zeros(preview.shape[1:], dtype=bool)
    return slide, im.detect_tissue(slide, image_key=IMAGE, inplace=False) > 0


def tissue_labels(image, mask):
    return Labels2DModel.parse(
        mask.astype(np.uint8),
        dims=("y", "x"),
        transformations=get_transformation(image, get_all=True),
    )


def align_mask(mask, base):
    shape = (base.sizes["y"], base.sizes["x"])
    chunks = base.data.chunksize[-2:]
    if mask.shape == shape:
        return da.from_array(mask, chunks=chunks)
    ratio = np.asarray(mask.shape) / shape
    return affine_transform(
        da.from_array(mask), np.diag(ratio), offset=(ratio - 1) / 2,
        output_shape=shape, output_chunks=chunks, order=0, mode="nearest",
    )


slide = read_zarr("source.zarr")
image = slide.images[IMAGE]
method = "${method}"
prefix = "${prefix}"
fit_reference = json.loads("${fit_reference}")

with dask.config.set(scheduler="threads", num_workers=int("${task.cpus}")):
    preview, mask = tissue_preview(image)
    if fit_reference:
        if not mask.any():
            raise ValueError("The reference slide '${sample_id}' has no detected tissue.")
        preview.labels[TISSUE] = tissue_labels(preview.images[IMAGE], mask)
        fit = im.fit_stain_reference(preview, image_key=IMAGE, method=method)
        record = {"reference_sample": "${sample_id}", "fit": asdict(fit)}
        Path(f"{prefix}.stain_reference.json").write_text(
            json.dumps(record, default=np.ndarray.tolist, allow_nan=False) + chr(10)
        )
    else:
        record = json.loads(Path("reference.json").read_text())
        fit = im.StainFit(**record["fit"])
        if fit.method != method:
            raise ValueError(f"The reference was fitted for {fit.method}, expected {method}.")
        if mask.any():
            if isinstance(image, DataTree) and len(image.children) > 1:
                # Squidpy materializes the coarse mask for statistics and streams the full-resolution mask.
                levels = (image["scale0"].image, image[f"scale{len(image.children) - 1}"].image)
                slide.labels[TISSUE] = DataTree.from_dict({
                    f"scale{index}": tissue_labels(level, align_mask(mask, level)).to_dataset(name="image")
                    for index, level in enumerate(levels)
                })
            else:
                base = image["scale0"].image if isinstance(image, DataTree) else image
                slide.labels[TISSUE] = tissue_labels(base, align_mask(mask, base))
            fit.transform(slide, image_key=IMAGE, image_key_added="normalized", output_dtype=np.uint8)
            image = slide.images["normalized"]
        else:
            record["skip_reason"] = "no_tissue"
        attrs = {**slide.attrs, "normalization": {"method": method, **record}}
        SpatialData(images={IMAGE: image}, attrs=attrs).write(f"{prefix}.normalized.zarr")

versions = {name: version(name) for name in ("squidpy", "spatialdata")}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
