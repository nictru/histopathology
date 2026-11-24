#!/usr/bin/env python3

import base64
import json
from html import escape
from importlib.metadata import version
from io import BytesIO
from pathlib import Path

import matplotlib
import numpy as np
import torch
from lazyslide import open_wsi, seg
from lazyslide_models import MODEL_REGISTRY
from matplotlib.figure import Figure
from shapely.geometry import box
from spatialdata import SpatialData, read_zarr

matplotlib.use("Agg")


def write_multiqc(
    filename, *, kind, section, description, sample, figure=None, note=""
):
    html = ""
    if figure is not None:
        output = BytesIO()
        figure.savefig(
            output, format="jpeg", dpi=110, pil_kwargs={"quality": 85, "optimize": True}
        )
        image = base64.b64encode(output.getvalue()).decode("ascii")
        alt = escape(f"{section} for {sample}", quote=True)
        html = (
            '<div class="mqc-custom-content-image">'
            f'<img src="data:image/jpeg;base64,{image}" alt="{alt}" /></div>'
        )
    if note:
        html += f"<p>{escape(note)}</p>"
    payload = {
        "id": Path(filename).stem.removesuffix("_mqc"),
        "parent_id": f"histopathology_{kind}",
        "parent_name": section,
        "parent_description": description,
        "section_name": sample,
        "plot_type": "image" if figure is not None else "html",
        "data": html,
    }
    Path(filename).write_text(json.dumps(payload), encoding="utf-8")


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


def plot_instances(wsi, shapes):
    # Center the crop on a detected instance so sparse slides still show a mask.
    center = shapes.geometry.iloc[len(shapes) // 2].representative_point()
    slide_height, slide_width = wsi.properties.shape
    width, height = min(512, slide_width), min(512, slide_height)
    x = max(0, min(int(center.x) - width // 2, slide_width - width))
    y = max(0, min(int(center.y) - height // 2, slide_height - height))
    indices = shapes.sindex.query(
        box(x, y, x + width, y + height), predicate="intersects"
    )
    instances = shapes.iloc[indices].assign(preview_color=np.arange(len(indices)) % 20)

    figure = Figure(figsize=(9, 4.5), layout="constrained")
    original, mask = figure.subplots(1, 2)
    original.imshow(wsi.reader.get_region(x, y, width, height), interpolation="nearest")
    instances.plot(
        ax=mask, column="preview_color", cmap="tab20", vmin=0, vmax=19, linewidth=0
    )
    mask.set_xlim(x, x + width)
    mask.set_ylim(y + height, y)
    for ax, title in ((original, "Slide detail"), (mask, "Instance mask")):
        ax.set_title(title)
        ax.set_axis_off()
    return (
        figure,
        f"{len(shapes):,} instances slide-wide. Crop origin: x={x:,}, y={y:,} (level 0 pixels).",
    )


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

instances = shapes.get(shape_key)
figure, note = (None, "No instances detected.")
if instances is not None and not instances.empty:
    figure, note = plot_instances(wsi, instances)
write_multiqc(
    "${prefix}.${model}.cells_mqc.json",
    kind="segmentation",
    section="Segmentation masks",
    description="Detail crops from the slide used for inference, with filled instance masks.",
    sample="${meta.id}" + f" · {name}",
    figure=figure,
    note=note,
)

packages = [
    "geopandas",
    "lazyslide",
    "lazyslide-models",
    "matplotlib",
    "spatialdata",
    "torch",
    "wsidata",
]
if name == "cellpose":
    packages.append("cellpose")
versions = {package: version(package) for package in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
