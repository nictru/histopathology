#!/usr/bin/env python3

import json
from importlib.metadata import version
from pathlib import Path

import torch
from lazyslide import open_wsi, tl
from lazyslide_models import TimmModel
from spatialdata import SpatialData, read_zarr

name = "${model}"
device = "${task.accelerator ? 'cuda' : 'cpu'}"
(checkpoint,) = Path("checkpoint").iterdir()
grid = read_zarr("tiles.zarr")
tile_count = len(grid.shapes.get("tiles", ()))
tables = {}

if tile_count:
    slide = read_zarr("source.zarr")
    slide.shapes.update(grid.shapes)
    slide.attrs["tile_spec"] = grid.attrs["tile_spec"]
    wsi = open_wsi(slide, image_key="wsi")
    torch.set_num_threads(int("${task.cpus}"))
    encoder = name
    if name == "resnet50":
        encoder = TimmModel(
            "resnet50.a1_in1k", pretrained_cfg_overlay={"file": str(checkpoint)}
        )
    key = f"{name}_tiles"
    tl.feature_extraction(
        wsi,
        model=encoder,
        model_path=checkpoint,
        tile_key="tiles",
        key_added=key,
        batch_size=int("${batch}"),
        device=device,
        pbar=False,
    )
    tables[key] = wsi.tables[key]

attrs = {
    **grid.attrs,
    "tile_embeddings": {
        "model": name,
        "checkpoint": checkpoint.name,
        "device": device,
        "tile_count": tile_count,
    },
}
SpatialData(shapes=dict(grid.shapes), tables=tables, attrs=attrs).write(
    "${prefix}.${model}.embeddings.zarr"
)

packages = ("lazyslide", "lazyslide-models", "spatialdata", "timm", "torch", "wsidata")
versions = {package: version(package) for package in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
