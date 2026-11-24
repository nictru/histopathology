#!/usr/bin/env python3

import json
import os
from importlib.metadata import version
from pathlib import Path

os.environ.setdefault("HF_HOME", str(Path(".huggingface").resolve()))
os.environ.setdefault("HF_HUB_DISABLE_XET", "1")
os.environ.setdefault("HF_HUB_DISABLE_IMPLICIT_TOKEN", "0")
os.environ.setdefault("TORCH_HOME", str(Path(".torch").resolve()))

import torch
from lazyslide import open_wsi, tl
from lazyslide_models import MODEL_REGISTRY, list_models
from spatialdata import SpatialData, read_zarr

name = "${model}"
device = "${task.accelerator ? 'cuda' : 'cpu'}"
encoder_class = MODEL_REGISTRY.get(name)
if name != "resnet50" and (
    name not in list_models("vision") + list_models("multimodal")
    or not callable(getattr(encoder_class, "encode_image", None))
):
    raise ValueError(f"Unsupported tile embedding model: {name}")
grid = read_zarr("tiles.zarr")
tile_count = len(grid.shapes.get("tiles", ()))
tables = {}

if tile_count:
    slide = read_zarr("source.zarr")
    slide.shapes.update(grid.shapes)
    slide.attrs["tile_spec"] = grid.attrs["tile_spec"]
    wsi = open_wsi(slide, image_key="wsi")
    torch.set_num_threads(int("${task.cpus}"))
    encoder = "resnet50.a1_in1k" if name == "resnet50" else name
    key = f"{name}_tiles"
    # Hugging Face reads HF_TOKEN from the task environment. Passing it to
    # the timm-backed registry models would save it to disk via login().
    tl.feature_extraction(
        wsi,
        model=encoder,
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
        "weight_source": "timm" if name == "resnet50" else "lazyslide_registry",
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
