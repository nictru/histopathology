#!/usr/bin/env python3

import base64
import json
import os
from html import escape
from importlib.metadata import version
from io import BytesIO
from math import ceil
from pathlib import Path

os.environ.setdefault("HF_HOME", str(Path(".huggingface").resolve()))
os.environ.setdefault("HF_HUB_DISABLE_XET", "1")

import dask
import dask.array as da
import matplotlib
import numpy as np
import zarr
from dask_image.ndinterp import affine_transform
from lazyslide import open_wsi, pp, seg
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from rasterio.features import rasterize
from rasterio.transform import Affine
from spatialdata import SpatialData, read_zarr
from spatialdata._io import write_labels
from spatialdata.models import Image2DModel, Labels2DModel
from spatialdata.transformations import get_transformation
from xarray import DataTree
from zarr.storage import FsspecStore

matplotlib.use("Agg")

REPORT_EDGE = 640


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


IMAGE = "wsi"
TISSUE = "wsi_tissue"
PREVIEW_EDGE = 2048


def preview_pixels(levels):
    preview = next(
        (
            level
            for level in levels
            if max(level.sizes["x"], level.sizes["y"]) <= PREVIEW_EDGE
        ),
        levels[-1],
    )
    factors = {axis: ceil(preview.sizes[axis] / PREVIEW_EDGE) for axis in ("y", "x")}
    if max(factors.values()) > 1:
        preview = (
            preview.coarsen(factors, boundary="pad").mean().round().astype(np.uint8)
        )
    return preview.data.compute()


def padded_grandqc():
    from lazyslide_models.segmentation import GrandQCTissue
    from torch.nn.functional import pad

    class PaddedGrandQCTissue(GrandQCTissue):
        def segment(self, image):
            height, width = image.shape[-2:]
            if min(height, width) >= 64:
                return super().segment(image)

            # The exported graph requires at least two 32-pixel blocks per axis.
            extra_h, extra_w = max(64 - height, 0), max(64 - width, 0)
            top, left = extra_h // 2, extra_w // 2
            padded = pad(image, (left, extra_w - left, top, extra_h - top))
            result = super().segment(padded)
            return result._replace(
                probability_map=result.probability_map[
                    ..., top : top + height, left : left + width
                ]
            )

    return PaddedGrandQCTissue()


def detect_mask(slide, full_level, pixels, method):
    shape = pixels.shape[1:]
    if not np.ptp(pixels, axis=(1, 2)).any():
        return np.zeros(shape, dtype=np.uint8)

    transform = Affine.identity()
    if method in ("otsu", "entropy"):
        wsi = open_wsi(
            SpatialData(images={IMAGE: Image2DModel.parse(pixels)}), image_key=IMAGE
        )
        pp.find_tissues(
            wsi,
            level=0,
            method=method,
            to_hsv=True,
            threshold=None,
            filter_artifacts=False,
            in_bounds=False,
        )
    else:
        import torch

        wsi = open_wsi(slide, image_key=IMAGE)
        if not wsi.properties.mpp:
            raise ValueError(
                f"{method} requires a slide pixel size. Supply mpp in the samplesheet."
            )
        torch.set_num_threads(int("${task.cpus}"))
        model = padded_grandqc() if method == "grandqc-tissue" else method
        seg.tissue(wsi, model=model, device="cpu", amp=False, compile=False)
        transform = Affine.scale(
            full_level.sizes["x"] / shape[1], full_level.sizes["y"] / shape[0]
        )

    tissues = wsi.shapes.get("tissues")
    if tissues is None or tissues.empty:
        return np.zeros(shape, dtype=np.uint8)
    return rasterize(
        tissues.geometry, out_shape=shape, transform=transform, dtype="uint8"
    )


def labels_for_level(mask, level):
    shape = (level.sizes["y"], level.sizes["x"])
    chunks = level.data.chunksize[-2:]
    if not mask.any():
        aligned = da.zeros(shape, chunks=chunks, dtype=np.uint8)
    elif mask.shape == shape:
        aligned = da.from_array(mask, chunks=chunks)
    else:
        ratio = np.asarray(mask.shape) / shape
        aligned = affine_transform(
            da.from_array(mask),
            np.diag(ratio),
            offset=(ratio - 1) / 2,
            output_shape=shape,
            output_chunks=chunks,
            order=0,
            mode="nearest",
        )
    return Labels2DModel.parse(
        aligned,
        dims=("y", "x"),
        transformations=get_transformation(level, get_all=True),
    )


method = "${method}"
slide = read_zarr("source.zarr")
image = slide.images[IMAGE]
levels = (
    [level.image for level in image.children.values()]
    if isinstance(image, DataTree)
    else [image]
)
with dask.config.set(scheduler="threads", num_workers=int("${task.cpus}")):
    preview = preview_pixels(levels)
    mask = detect_mask(slide, levels[0], preview, method)
    labels = labels_for_level(mask, levels[0])
    if len(levels) > 1:
        # Squidpy samples at the coarsest level and applies the mask at full resolution.
        labels = DataTree.from_dict(
            {
                "scale0": labels.to_dataset(name="image"),
                "scale1": labels_for_level(mask, levels[-1]).to_dataset(name="image"),
            }
        )
    detection = {"method": method, "has_tissue": bool(mask.any())}
    if method == "otsu":
        detection["channel"] = "saturation"
    attrs = {"tissue_detection": detection}
    output = "${mask_prefix}.tissue.zarr"
    output_path = Path(output)
    zarr_store = FsspecStore.from_url(str(output_path.resolve()), storage_options={"auto_mkdir": True})
    root = zarr.create_group(store=zarr_store, zarr_format=3)
    SpatialData(attrs=attrs).write_attrs(zarr_group=root)
    write_labels(
        labels,
        root,
        TISSUE,
        storage_options={"chunks": (256, 256), "shards": (2048, 2048)},
    )
    zarr.consolidate_metadata(zarr_store)

step = ceil(max(mask.shape) / REPORT_EDGE)
figure = Figure(figsize=(6, 5), layout="constrained")
ax = figure.subplots()
ax.imshow(preview[:, ::step, ::step].transpose(1, 2, 0), interpolation="nearest")
ax.imshow(
    mask[::step, ::step],
    cmap=ListedColormap(["none", "#22a699"]),
    vmin=0,
    vmax=1,
    alpha=0.4,
    interpolation="nearest",
)
ax.set_axis_off()
write_multiqc(
    "${mask_prefix}.tissue_mqc.json",
    kind="tissue",
    section="Tissue masks",
    description="Tissue segmentation masks overlaid on the original slides.",
    sample="${meta.id}" + f" · {method}",
    figure=figure,
    note="" if mask.any() else "No tissue detected.",
)

packages = ["lazyslide", "matplotlib", "spatialdata", "wsidata"]
if method not in ("otsu", "entropy"):
    packages.extend(["lazyslide-models", "torch"])
versions = {name: version(name) for name in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
