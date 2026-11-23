#!/usr/bin/env python3

import json
from importlib.metadata import version
from pathlib import Path

import torch
from lazyslide import open_wsi, seg
from lazyslide_models.segmentation import GrandQCArtifact
from spatialdata import SpatialData, read_zarr


class CheckpointGrandQC(GrandQCArtifact):
    def __init__(self):
        self.model = torch.export.load("checkpoint.pt2").module()


def detect_artifacts(tiles, variant):
    slide = read_zarr("source.zarr")
    slide.shapes["tiles"] = tiles.shapes["tiles"]
    slide.attrs["tile_spec"] = tiles.attrs["tile_spec"]
    wsi = open_wsi(slide, image_key="wsi")
    torch.set_num_threads(int("${task.cpus}"))
    seg.artifact(wsi, tile_key="tiles", model=CheckpointGrandQC(), variant=variant, device="cpu", pbar=False)
    return wsi.shapes.get("artifacts")


variant = "${variant}"
tiles = read_zarr("tiles.zarr")
tile_grid = tiles.shapes.get("tiles")
tile_count = len(tile_grid) if tile_grid is not None else 0
shapes = {}
counts = {}
if tile_count:
    artifacts = detect_artifacts(tiles, variant)
    if artifacts is not None:
        shapes["artifacts"] = artifacts
        counts = artifacts["class"].value_counts().to_dict()

attrs = {"grandqc": {"variant": variant, "tile_count": tile_count, "artifact_regions": counts}}
SpatialData(shapes=shapes, attrs=attrs).write("${prefix}.grandqc.zarr")

versions = {name: version(name) for name in ("lazyslide", "lazyslide-models", "spatialdata", "wsidata")}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
