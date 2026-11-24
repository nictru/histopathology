#!/usr/bin/env python3

import base64
import json
from html import escape
from importlib.metadata import version
from io import BytesIO
from math import ceil
from pathlib import Path

import datashader as ds
import matplotlib
import numpy as np
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
from matplotlib.figure import Figure
from spatialdata import read_zarr
from xarray import DataTree

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


def tile_model(method):
    if method == "tileqc":
        return CVCompose(Brightness(), Contrast(), Sharpness(), Saturation())

    # Preserve the registry model's transforms and scoring with the staged export.
    class CheckpointModel(MODEL_REGISTRY[method]):
        def __init__(self):
            self.model = torch.export.load("checkpoint.pt2").module()

    return CheckpointModel()


def score_tiles(wsi, tiles, method):
    wsi.shapes["tiles"] = tiles.shapes["tiles"]
    wsi.attrs["tile_spec"] = tiles.attrs["tile_spec"]
    torch.set_num_threads(int("${task.cpus}"))
    tl.tile_prediction(
        wsi, tile_model(method), tile_key="tiles", device="cpu", pbar=False
    )
    return wsi.shapes["tiles"]


def plot_scores(wsi, tiles):
    if tiles is None or tiles.empty:
        return None, "No tiles available for scoring."

    metrics = (
        tiles.select_dtypes(include="number")
        .drop(columns=["tile_id", "tissue_id"], errors="ignore")
        .replace([np.inf, -np.inf], np.nan)
        .dropna(axis="columns", how="all")
    )
    if metrics.empty:
        return None, "No finite tile scores."

    height, width = wsi.properties.shape
    image = wsi.images["wsi"]
    levels = (
        [node.image for node in image.children.values()]
        if isinstance(image, DataTree)
        else [image]
    )
    level = next(
        (level for level in levels if max(level.sizes["y"], level.sizes["x"]) <= 640),
        levels[-1],
    )
    step = ceil(max(level.sizes["y"], level.sizes["x"]) / 640)
    preview = (
        level.isel(y=slice(None, None, step), x=slice(None, None, step))
        .transpose("y", "x", "c")
        .data.compute()
    )
    canvas = ds.Canvas(
        plot_width=preview.shape[1],
        plot_height=preview.shape[0],
        x_range=(0, width),
        y_range=(0, height),
    )
    heatmaps = canvas.polygons(
        tiles[["geometry"]].join(metrics),
        geometry="geometry",
        agg=ds.summary(**{metric: ds.mean(metric) for metric in metrics}),
    )

    cols = min(len(metrics.columns), 3)
    rows = ceil(len(metrics.columns) / cols)
    figure = Figure(figsize=(4 * cols, 4 * rows), layout="constrained")
    axes = figure.subplots(rows, cols, squeeze=False).ravel()
    for ax in axes:
        ax.set_axis_off()
    for ax, (metric, values) in zip(axes, metrics.items()):
        ax.imshow(preview, interpolation="nearest")
        title = metric.replace("_", " ").capitalize()
        low, high = values.min(), values.max()
        display = ax.imshow(
            heatmaps[metric],
            cmap="viridis",
            vmin=low,
            vmax=high if high > low else low + 1,
            alpha=0.8,
            interpolation="nearest",
        )
        ax.set_title(title if high > low else f"{title} ({low:.3g})")
        if high > low:
            figure.colorbar(display, ax=ax, shrink=0.7)
    return figure, f"Tiles scored: {len(tiles):,}."


method = "${method}"
wsi = open_wsi(read_zarr("source.zarr"), image_key="wsi")
prefix = "${prefix}"
tiles = read_zarr("tiles.zarr")
if "tiles" in tiles.shapes and len(tiles.shapes["tiles"]):
    tiles.shapes["tiles"] = score_tiles(wsi, tiles, method)
tiles.attrs["tile_qc"] = {"method": method}
tiles.write(f"{prefix}.{method}.zarr")

method_names = {
    "tileqc": "Image metrics",
    "focuslitenn": "FocusLiteNN",
    "pathprofilerqc": "PathProfilerQC",
}
figure, note = plot_scores(wsi, tiles.shapes.get("tiles"))
write_multiqc(
    f"{prefix}.{method}_mqc.json",
    kind="tile_qc",
    section="Tile QC",
    description="Spatial maps of tile scores. Overlapping scores are averaged and unscored areas show the slide.",
    sample="${meta.id}" + f" · {method_names[method]}",
    figure=figure,
    note=note,
)

versions = {
    name: version(name)
    for name in (
        "datashader",
        "lazyslide",
        "lazyslide-models",
        "matplotlib",
        "spatialdata",
        "wsidata",
    )
}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
