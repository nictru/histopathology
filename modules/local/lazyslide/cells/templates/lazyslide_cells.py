#!/usr/bin/env python3

import base64
import json
import os
from html import escape
from importlib.metadata import version
from io import BytesIO
from pathlib import Path

os.environ.setdefault("HF_HOME", str(Path(".huggingface").resolve()))
os.environ.setdefault("HF_HUB_DISABLE_XET", "1")
os.environ.setdefault("HF_HUB_DISABLE_IMPLICIT_TOKEN", "0")
os.environ.setdefault("CELLPOSE_LOCAL_MODELS_PATH", str(Path(".cellpose").resolve()))

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


def load_model(name, device):
    if name == "cellpose":
        from cellpose.models import CellposeModel

        class DeviceCellpose(MODEL_REGISTRY["cellpose"]):
            def __init__(self):
                self.model = CellposeModel(
                    pretrained_model="cpsam_v2",
                    device=device,
                    use_bfloat16=device.type == "cpu"
                    or torch.cuda.is_bf16_supported(including_emulation=False),
                )
                self.eval_kwargs = {}

        return DeviceCellpose()
    if name == "nulite" and device.type == "cuda":
        from huggingface_hub import hf_hub_download
        from torch.export.passes import move_to_device_pass

        class CudaNuLite(MODEL_REGISTRY["nulite"]):
            def __init__(self):
                weights = hf_hub_download(
                    "RendeiroLab/LazySlide-models", "NuLite/NuLite_H_exported.pt2"
                )
                # Move device constants in the exported graph as well as its parameters.
                self.model = move_to_device_pass(
                    torch.export.load(weights), device
                ).module()
                self.magnification = "20x"

        return CudaNuLite()
    return MODEL_REGISTRY[name]()


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

    figure = Figure(figsize=(5, 5), layout="constrained")
    ax = figure.subplots()
    ax.imshow(
        wsi.reader.get_region(x, y, width, height),
        extent=(x, x + width, y + height, y),
        origin="upper",
        interpolation="nearest",
    )
    instances.plot(
        ax=ax,
        column="preview_color",
        cmap="tab20",
        vmin=0,
        vmax=19,
        linewidth=0,
        alpha=0.45,
    )
    ax.set(xlim=(x, x + width), ylim=(y + height, y))
    ax.set_axis_off()
    return (
        figure,
        f"{len(shapes):,} instances slide-wide. Crop origin: x={x:,}, y={y:,} (level 0 pixels).",
    )


name = "${model}"
device = torch.device("${task.accelerator ? 'cuda' : 'cpu'}")
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
        model=load_model(name, device),
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
    "weight_source": "cellpose" if name == "cellpose" else "lazyslide_registry",
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
    original = open_wsi(read_zarr("original.zarr"), image_key="wsi")
    figure, note = plot_instances(original, instances)
write_multiqc(
    "${prefix}.${model}.cells_mqc.json",
    kind="segmentation",
    section="Segmentation masks",
    description="Cell or nucleus masks overlaid on detail crops from the original slides.",
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
