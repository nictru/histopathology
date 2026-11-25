#!/usr/bin/env python3

import json
import os
from importlib.metadata import version
from pathlib import Path

os.environ.setdefault("HF_HOME", str(Path(".huggingface").resolve()))
os.environ.setdefault("HF_HUB_DISABLE_XET", "1")

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from lazyslide import open_wsi, tl
from lazyslide_models import MODEL_REGISTRY
from spatialdata import read_zarr

sample = "${meta.id}"
tile_model = "${tile_model}"
method = "${method}"
store = read_zarr("embeddings.zarr")
tile_info = store.attrs["tile_embeddings"]
tile_count = tile_info["tile_count"]
if tile_info["model"] != tile_model:
    raise ValueError(f"Expected {tile_model} tile embeddings, found {tile_info['model']}")

aggregated = np.empty((0, 0), dtype=np.float32)
if tile_count:
    key = f"{tile_model}_tiles"
    table = store.tables[key]
    tiles = store.shapes["tiles"].set_index("tile_id", drop=False)
    ids = pd.Index(table.obs["tile_id"])
    if (
        len(tiles) != tile_count
        or len(ids) != tile_count
        or not tiles.index.is_unique
        or not ids.is_unique
    ):
        raise ValueError(f"Tile IDs do not identify {tile_count} unique tile embeddings")

    # The LazySlide aggregation API reads coordinates and feature rows positionally.
    wsi = open_wsi(read_zarr("source.zarr"), image_key="wsi")
    wsi.shapes["tiles"] = tiles.loc[ids].reset_index(drop=True)
    wsi.tables[key] = table
    wsi.attrs["tile_spec"] = store.attrs["tile_spec"]
    registry_key = "chief-slide-encoder" if method == "chief" else method
    if registry_key in MODEL_REGISTRY:
        expected = MODEL_REGISTRY[registry_key].vision_encoder
        if expected != tile_model:
            raise ValueError(f"{method} requires {expected} tile embeddings")
    torch.set_num_threads(int("${task.cpus}"))
    tl.feature_aggregation(
        wsi,
        feature_key=key,
        encoder=method,
        device="cpu",
        amp=False,
        compile=False,
    )
    aggregated = np.asarray(table.uns["agg_ops"]["agg_slide"]["features"], dtype=np.float32)
    if aggregated.ndim != 2 or aggregated.shape[0] != 1 or not np.isfinite(aggregated).all():
        raise ValueError(f"{method} did not return one finite embedding for {sample}")

result = AnnData(aggregated, obs=pd.DataFrame(index=[sample] if tile_count else []))
result.uns["slide_embedding"] = {
    "sample_id": sample,
    "tile_model": tile_model,
    "method": method,
    "tile_count": tile_count,
}
if tile_count:
    result.uns["tile_spec"] = store.attrs["tile_spec"]["tiles"]
result.write_h5ad("${prefix}.${tile_model}.${method}.slide_embedding.h5ad")

packages = ("anndata", "lazyslide", "lazyslide-models", "numpy", "spatialdata", "torch", "wsidata")
versions = {package: version(package) for package in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
