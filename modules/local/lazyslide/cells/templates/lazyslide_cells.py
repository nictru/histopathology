#!/usr/bin/env python3

import json
from importlib.metadata import version
from pathlib import Path

import torch
from lazyslide import open_wsi, seg
from lazyslide_models import MODEL_REGISTRY
from spatialdata import SpatialData, read_zarr


def load_model(name, checkpoint, device):
    # Load staged weights and select the device before initializing the backend.
    class StagedModel(MODEL_REGISTRY[name]):
        def __init__(self):
            if name == "cellpose":
                from cellpose.models import CellposeModel

                self.model = CellposeModel(
                    pretrained_model=str(checkpoint),
                    device=device,
                    use_bfloat16=device.type == "cpu"
                    or torch.cuda.is_bf16_supported(including_emulation=False),
                )
                self.eval_kwargs = {}
            elif name == "instanseg":
                self.model = torch.jit.load(checkpoint, map_location=device).eval()
            elif name == "nulite":
                from torch.export.passes import move_to_device_pass

                self.model = move_to_device_pass(
                    torch.export.load(checkpoint), device
                ).module()
                self.magnification = "20x"
            else:
                from lazyslide_models.segmentation.cellvit_family.histoplus import (
                    HistoPLUSModel,
                    remap_state_dict,
                )

                self.variant = "20x"
                self.model = HistoPLUSModel(backbone_tile_size=224)
                state = torch.load(checkpoint, map_location="cpu", weights_only=True)
                self.model.load_state_dict(remap_state_dict(state))
                self.model.eval()

    return StagedModel()


name = "${model}"
device = torch.device("${task.accelerator ? 'cuda' : 'cpu'}")
checkpoint = next(Path("checkpoint").iterdir())
shape_key = "cells" if name == "cellpose" else "nuclei"
grid = read_zarr("tiles.zarr")
tile_spec = grid.attrs["tile_spec"].get("tiles")
shapes = {}

if tile_spec:
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError(
            "Segmentation requested CUDA, but no CUDA device is visible to the task."
        )
    slide = read_zarr("source.zarr")
    slide.shapes["tiles"] = grid.shapes["tiles"]
    tissue_key = tile_spec["tissue_name"]
    slide.shapes[tissue_key] = grid.shapes[tissue_key]
    slide.attrs["tile_spec"] = grid.attrs["tile_spec"]
    wsi = open_wsi(slide, image_key="wsi")
    torch.set_num_threads(int("${task.cpus}"))
    seg.cells(
        wsi,
        model=load_model(name, checkpoint, device),
        batch_size=1,
        device=device,
        amp=False,
        pbar=False,
        overlap_ownership=True,
        key_added=shape_key,
    )
    if shape_key in wsi.shapes:
        shapes[shape_key] = wsi.shapes[shape_key]

segmentation = {
    "model": name,
    "device": device.type,
    "checkpoint": checkpoint.name,
    "instance_type": "model_defined" if name == "cellpose" else "nucleus",
    "instance_count": len(shapes[shape_key]) if shapes else 0,
    "tile_spec": tile_spec,
}
SpatialData(shapes=shapes, attrs={"segmentation": segmentation}).write(
    "${prefix}.${model}.cells.zarr"
)

packages = ["lazyslide", "lazyslide-models", "spatialdata", "torch", "wsidata"]
if name == "cellpose":
    packages.append("cellpose")
versions = {package: version(package) for package in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
