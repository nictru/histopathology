#!/usr/bin/env python3

import json
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path

import dask
import numpy as np
from spatialdata import SpatialData, read_zarr
from squidpy.experimental import im

IMAGE = "wsi"
TISSUE = "wsi_tissue"

slide = read_zarr("source.zarr")
mask = read_zarr("tissue.zarr")
slide.labels[TISSUE] = mask.labels[TISSUE]
has_tissue = mask.attrs["tissue_detection"]["has_tissue"]
method = "${method}"
prefix = "${prefix}"
fit_reference = json.loads("${fit_reference}")

with dask.config.set(scheduler="threads", num_workers=int("${task.cpus}")):
    if fit_reference:
        if not has_tissue:
            raise ValueError("The reference slide '${sample_id}' has no detected tissue.")
        fit = im.fit_stain_reference(slide, image_key=IMAGE, method=method)
        record = {"reference_sample": "${sample_id}", "fit": asdict(fit)}
        Path(f"{prefix}.stain_reference.json").write_text(
            json.dumps(record, default=np.ndarray.tolist, allow_nan=False) + chr(10)
        )
    else:
        record = json.loads(Path("reference.json").read_text())
        fit = im.StainFit(**record["fit"])
        if fit.method != method:
            raise ValueError(f"The reference was fitted for {fit.method}, expected {method}.")
        image = slide.images[IMAGE]
        if has_tissue:
            fit.transform(slide, image_key=IMAGE, image_key_added="normalized", output_dtype=np.uint8)
            image = slide.images["normalized"]
        else:
            record["skip_reason"] = "no_tissue"
        attrs = {**slide.attrs, "normalization": {"method": method, **record}}
        SpatialData(images={IMAGE: image}, attrs=attrs).write(f"{prefix}.normalized.zarr")

versions = {name: version(name) for name in ("squidpy", "spatialdata")}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
