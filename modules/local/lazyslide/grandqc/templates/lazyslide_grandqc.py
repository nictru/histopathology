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

import matplotlib
import matplotlib.pyplot as plt
import torch
from lazyslide import open_wsi, pl, seg
from lazyslide_models.segmentation import GrandQCArtifact
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


def detect_artifacts(wsi, tiles, variant):
    wsi.shapes["tiles"] = tiles.shapes["tiles"]
    wsi.attrs["tile_spec"] = tiles.attrs["tile_spec"]
    torch.set_num_threads(int("${task.cpus}"))
    seg.artifact(
        wsi,
        tile_key="tiles",
        model=GrandQCArtifact(variant=variant),
        variant=variant,
        device="cpu",
        pbar=False,
    )
    return wsi.shapes.get("artifacts")


def plot_artifacts(wsi, artifact_regions, tile_count, variant):
    figure, ax = plt.subplots(figsize=(6, 5), layout="constrained")
    viewer = pl.WSIViewer(
        wsi, in_bounds=False, img_bytes_limit=640 * 640 * 3, target_dpi=110
    )
    viewer.add_image()
    viewer.title = f"GrandQC {variant}"
    if artifact_regions is not None and not artifact_regions.empty:
        viewer.add_polygons("artifacts", color_by="class", alpha=0.55, linewidth=0)
        note = f"Tiles evaluated: {tile_count:,}."
    elif tile_count:
        note = f"No artifact regions detected across {tile_count:,} evaluated tiles."
    else:
        note = "No tiles available for GrandQC."
    viewer.show(ax=ax, axis="off")
    return figure, note


variant = "${variant}"
wsi = open_wsi(read_zarr("source.zarr"), image_key="wsi")
tiles = read_zarr("tiles.zarr")
tile_grid = tiles.shapes.get("tiles")
tile_count = len(tile_grid) if tile_grid is not None else 0
shapes = {}
counts = {}
if tile_count:
    artifacts = detect_artifacts(wsi, tiles, variant)
    if artifacts is not None:
        shapes["artifacts"] = artifacts
        counts = artifacts["class"].value_counts().to_dict()

attrs = {
    "grandqc": {
        "variant": variant,
        "tile_count": tile_count,
        "artifact_regions": counts,
    }
}
SpatialData(shapes=shapes, attrs=attrs).write("${prefix}.grandqc.zarr")

figure, note = plot_artifacts(wsi, shapes.get("artifacts"), tile_count, variant)
write_multiqc(
    "${prefix}.grandqc_mqc.json",
    kind="slide_qc",
    section="Slide QC",
    description="GrandQC artifact classes overlaid on the slide.",
    sample="${meta.id}" + f" · GrandQC {variant}",
    figure=figure,
    note=note,
)

versions = {
    name: version(name)
    for name in (
        "lazyslide",
        "lazyslide-models",
        "matplotlib",
        "spatialdata",
        "wsidata",
    )
}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
