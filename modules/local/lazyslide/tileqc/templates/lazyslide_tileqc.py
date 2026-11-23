#!/usr/bin/env python3

import json
import os
from importlib.metadata import version
from pathlib import Path

from spatialdata import read_zarr


def score_tiles(tiles, method):
    if method != "tileqc":
        os.environ["HF_HUB_CACHE"] = str(Path("hf-cache").absolute())
        os.environ["HF_HUB_OFFLINE"] = "1"

    import torch
    from lazyslide import open_wsi, tl

    model = method
    if method == "tileqc":
        from lazyslide_models.tile_prediction.cv_features import (
            Brightness,
            Contrast,
            CVCompose,
            Saturation,
            Sharpness,
        )

        model = CVCompose(Brightness(), Contrast(), Sharpness(), Saturation())

    slide = read_zarr("source.zarr")
    slide.shapes["tiles"] = tiles.shapes["tiles"]
    slide.attrs["tile_spec"] = tiles.attrs["tile_spec"]
    wsi = open_wsi(slide, image_key="wsi")
    torch.set_num_threads(int("${task.cpus}"))
    tl.tile_prediction(wsi, model, tile_key="tiles", device="cpu", pbar=False)
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
