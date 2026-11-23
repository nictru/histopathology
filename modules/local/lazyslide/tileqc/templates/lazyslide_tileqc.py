#!/usr/bin/env python3

import json
from importlib.metadata import version
from pathlib import Path

import torch
from lazyslide import open_wsi, tl
from lazyslide_models import MODEL_REGISTRY
from lazyslide_models.tile_prediction.cv_features import (
    Brightness,
    Contrast,
    CVCompose,
    Saturation,
    Sharpness,
)
from spatialdata import read_zarr


def tile_model(method):
    if method == "tileqc":
        return CVCompose(Brightness(), Contrast(), Sharpness(), Saturation())

    # Preserve the registry model's transforms and scoring with the staged export.
    class CheckpointModel(MODEL_REGISTRY[method]):
        def __init__(self):
            self.model = torch.export.load("checkpoint.pt2").module()

    return CheckpointModel()


def score_tiles(tiles, method):
    slide = read_zarr("source.zarr")
    slide.shapes["tiles"] = tiles.shapes["tiles"]
    slide.attrs["tile_spec"] = tiles.attrs["tile_spec"]
    wsi = open_wsi(slide, image_key="wsi")
    torch.set_num_threads(int("${task.cpus}"))
    tl.tile_prediction(wsi, tile_model(method), tile_key="tiles", device="cpu", pbar=False)
    return wsi.shapes["tiles"]


method = "${method}"
prefix = "${prefix}"
tiles = read_zarr("tiles.zarr")
if "tiles" in tiles.shapes and len(tiles.shapes["tiles"]):
    tiles.shapes["tiles"] = score_tiles(tiles, method)
tiles.attrs["tile_qc"] = {"method": method}
tiles.write(f"{prefix}.{method}.zarr")

versions = {name: version(name) for name in ("lazyslide", "lazyslide-models", "spatialdata", "wsidata")}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
