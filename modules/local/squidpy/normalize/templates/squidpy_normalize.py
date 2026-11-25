#!/usr/bin/env python3

import base64
import json
from dataclasses import asdict
from html import escape
from importlib.metadata import version
from io import BytesIO
from math import ceil
from pathlib import Path

import dask
import matplotlib
import numpy as np
import zarr
from matplotlib.figure import Figure
from spatialdata import SpatialData, read_zarr
from spatialdata._io import write_image
from squidpy.experimental import im
from xarray import DataTree
from zarr.storage import FsspecStore

matplotlib.use("Agg")

PREVIEW_EDGE = 640


def thumbnail(image):
    levels = (
        [node.image for node in image.children.values()]
        if isinstance(image, DataTree)
        else [image]
    )
    edge = next(
        (
            level
            for level in levels
            if max(level.sizes["y"], level.sizes["x"]) <= PREVIEW_EDGE
        ),
        levels[-1],
    )
    step = ceil(max(edge.sizes["y"], edge.sizes["x"]) / PREVIEW_EDGE)
    return (
        edge.isel(y=slice(None, None, step), x=slice(None, None, step))
        .transpose("y", "x", "c")
        .data.compute()
    )


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

slide = read_zarr("source.zarr")
mask = read_zarr("tissue.zarr")
slide.labels[TISSUE] = mask.labels[TISSUE]
has_tissue = mask.attrs["tissue_detection"]["has_tissue"]
method = "${method}"
prefix = "${prefix}"
fit_reference = json.loads("${fit_reference}")

with dask.config.set(scheduler="threads", num_workers=int("${task.cpus}")):
    if fit_reference:
        if not has_tissue:
            raise ValueError(
                "The reference slide '${sample_id}' has no detected tissue."
            )
        fit = im.fit_stain_reference(slide, image_key=IMAGE, method=method)
        record = {"reference_sample": "${sample_id}", "fit": asdict(fit)}
        Path(f"{prefix}.stain_reference.json").write_text(
            json.dumps(record, default=np.ndarray.tolist, allow_nan=False) + chr(10)
        )
    else:
        record = json.loads(Path("reference.json").read_text())
        fit = im.StainFit(**record["fit"])
        if fit.method != method:
            raise ValueError(
                f"The reference was fitted for {fit.method}, expected {method}."
            )
        original = slide.images[IMAGE]
        image = original
        if has_tissue:
            fit.transform(
                slide,
                image_key=IMAGE,
                image_key_added="normalized",
                output_dtype=np.uint8,
            )
            image = slide.images["normalized"]
        else:
            record["skip_reason"] = "no_tissue"
        attrs = {**slide.attrs, "normalization": {"method": method, **record}}
        output = f"{prefix}.normalized.zarr"
        output_path = Path(output)
        zarr_store = FsspecStore.from_url(str(output_path.resolve()), storage_options={"auto_mkdir": True})
        root = zarr.create_group(store=zarr_store, zarr_format=3)
        SpatialData(attrs=attrs).write_attrs(zarr_group=root)
        # SpatialData.write() does not expose raster storage options.
        write_image(
            image,
            root.require_group(f"images/{IMAGE}"),
            IMAGE,
            storage_options={"chunks": (256, 256), "shards": (2048, 2048)},
        )
        zarr.consolidate_metadata(zarr_store)
        original_preview = thumbnail(original)
        normalized_preview = thumbnail(image) if has_tissue else original_preview
        figure = Figure(figsize=(9, 4.5), layout="constrained")
        axes = figure.subplots(1, 2)
        labels = (
            "Original",
            f"Normalized ({method})" if has_tissue else "Original retained",
        )
        for ax, preview, title in zip(
            axes, (original_preview, normalized_preview), labels
        ):
            ax.imshow(preview)
            ax.set_title(title)
            ax.set_axis_off()
        write_multiqc(
            f"{prefix}.normalized_mqc.json",
            kind="normalization",
            section="Stain normalization",
            description="Matched thumbnails before and after stain normalization.",
            sample="${sample_id}",
            figure=figure,
            note=""
            if has_tissue
            else "Normalization skipped because no tissue was detected.",
        )

versions = {name: version(name) for name in ("matplotlib", "squidpy", "spatialdata")}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
