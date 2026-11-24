#!/usr/bin/env python3

import argparse
import base64
import json
import os
import shlex
from html import escape
from importlib.metadata import version
from io import BytesIO
from pathlib import Path

os.environ.setdefault("HF_HOME", str(Path(".huggingface").resolve()))
os.environ.setdefault("HF_HUB_DISABLE_XET", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import cv2
import numpy as np
import timm
import torch
from huggingface_hub import hf_hub_download
from huggingface_hub.errors import GatedRepoError
from lazyslide import open_wsi
from matplotlib import colormaps
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from PIL import Image
from shapely.geometry import box
from spatialdata import read_zarr
from timm.data import IMAGENET_DEFAULT_MEAN, IMAGENET_DEFAULT_STD
from torch import nn

MODEL_REPO = "zhihuanglab/CytoFormer"
MODEL_REVISION = "6a3572026e491774e7530beb3485e31459fdbee5"
IMAGE_SIZE = 224
FIELD_OF_VIEW_UM = 56
TYPE_COLUMN = "cytoformer_cell_type"
PROBABILITY_COLUMN = "cytoformer_probability"


def load_classifier(organ, device):
    token = os.environ.get("HF_TOKEN")
    try:
        mapping = hf_hub_download(
            MODEL_REPO, "organ_celltype_map.json", revision=MODEL_REVISION, token=token
        )
        taxonomy = json.loads(Path(mapping).read_text(encoding="utf-8"))
        if organ not in taxonomy["organs"]:
            supported_organs = ", ".join(taxonomy["organs"])
            raise ValueError(
                f"Unsupported CytoFormer organ {organ!r}. Choose from {supported_organs}."
            )
        weights = hf_hub_download(
            MODEL_REPO, "checkpoint.pth", revision=MODEL_REVISION, token=token
        )
    except GatedRepoError as exc:
        raise RuntimeError(
            f"Request access to {MODEL_REPO} on Hugging Face and configure the "
            "HF_TOKEN Nextflow secret with a read token for that repository."
        ) from exc

    classes = np.asarray(taxonomy["organ_to_celltypes"][organ])
    head_index = taxonomy["organs"].index(organ)
    checkpoint = torch.load(weights, map_location="cpu", weights_only=True)
    # One task uses one organ head. Evaluation of the upstream dropout is an identity.
    head_key = f"head.heads.{head_index}."
    state = {
        name: value
        for name, value in checkpoint.items()
        if name.startswith(("uni_encoder.", "feat_norm."))
    }
    state.update(
        {
            f"head.{parameter}": checkpoint[head_key + parameter]
            for parameter in ("weight", "bias")
        }
    )
    del checkpoint
    # Match the UNI2-h encoder in the CytoFormer checkpoint without loading foundation weights.
    with torch.device("meta"):
        encoder = timm.create_model(
            "vit_giant_patch14_224",
            pretrained=False,
            img_size=IMAGE_SIZE,
            patch_size=14,
            depth=24,
            num_heads=24,
            init_values=1e-5,
            embed_dim=1536,
            mlp_ratio=2.66667 * 2,
            num_classes=0,
            no_embed_class=True,
            mlp_layer=timm.layers.SwiGLUPacked,
            act_layer=nn.SiLU,
            reg_tokens=8,
            dynamic_img_size=True,
        )
        model = nn.Sequential()
        model.add_module("uni_encoder", encoder)
        model.add_module("feat_norm", nn.LayerNorm(encoder.embed_dim))
        model.add_module("head", nn.Linear(encoder.embed_dim, len(classes)))
    model.load_state_dict(state, strict=True, assign=True)
    return model.to(device).eval(), classes


def cell_patch(wsi, center_x, center_y, size):
    height, width = wsi.properties.shape
    x, y = round(center_x) - size // 2, round(center_y) - size // 2
    left, top = max(x, 0), max(y, 0)
    right, bottom = min(x + size, width), min(y + size, height)
    patch = Image.new("RGB", (size, size), "white")
    if right > left and bottom > top:
        region = wsi.reader.get_region(left, top, right - left, bottom - top)
        patch.paste(Image.fromarray(region).convert("RGB"), (left - x, top - y))
    return cv2.resize(
        np.asarray(patch), (IMAGE_SIZE, IMAGE_SIZE), interpolation=cv2.INTER_AREA
    )


@torch.inference_mode()
def classify(instances, wsi, organ, device, batch_size):
    model, classes = load_classifier(organ, device)
    mean = torch.tensor(IMAGENET_DEFAULT_MEAN, device=device).view(1, 3, 1, 1)
    std = torch.tensor(IMAGENET_DEFAULT_STD, device=device).view(1, 3, 1, 1)
    size = 2 * max(8, round(FIELD_OF_VIEW_UM / (2 * wsi.properties.mpp)))
    centers = instances.geometry.centroid
    predicted_types, probabilities = [], []
    use_bfloat16 = device.type == "cuda" and torch.cuda.is_bf16_supported()
    for start in range(0, len(centers), batch_size):
        batch = centers.iloc[start : start + batch_size]
        images = np.stack([cell_patch(wsi, point.x, point.y, size) for point in batch])
        inputs = torch.from_numpy(images).to(device).permute(0, 3, 1, 2).float() / 255
        with torch.autocast(device.type, dtype=torch.bfloat16, enabled=use_bfloat16):
            logits = model((inputs - mean) / std)
        confidence, labels = logits.float().softmax(dim=1).max(dim=1)
        predicted_types.extend(classes[labels.cpu().numpy()])
        probabilities.extend(confidence.cpu().numpy())
    instances[TYPE_COLUMN] = predicted_types
    instances[PROBABILITY_COLUMN] = probabilities
    return classes


def plot_classes(wsi, instances, classes):
    center = instances.geometry.iloc[len(instances) // 2].representative_point()
    slide_height, slide_width = wsi.properties.shape
    width, height = min(512, slide_width), min(512, slide_height)
    x = max(0, min(int(center.x) - width // 2, slide_width - width))
    y = max(0, min(int(center.y) - height // 2, slide_height - height))
    positions = instances.sindex.query(
        box(x, y, x + width, y + height), predicate="intersects"
    )
    visible = instances.iloc[positions]
    palette = colormaps["tab10"].colors
    colors = dict(zip(classes, palette))

    figure = Figure(figsize=(5, 5), layout="constrained")
    ax = figure.subplots()
    ax.imshow(
        wsi.reader.get_region(x, y, width, height),
        extent=(x, x + width, y + height, y),
        interpolation="nearest",
    )
    visible.plot(
        ax=ax, color=visible[TYPE_COLUMN].map(colors).tolist(), linewidth=0, alpha=0.45
    )
    ax.legend(
        handles=[
            Patch(facecolor=colors[label], label=label)
            for label in sorted(visible[TYPE_COLUMN].unique())
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        ncols=2,
        frameon=False,
    )
    ax.set(xlim=(x, x + width), ylim=(y + height, y))
    ax.set_axis_off()

    output = BytesIO()
    figure.savefig(output, format="jpeg", dpi=110, pil_kwargs={"quality": 85})
    return base64.b64encode(output.getvalue()).decode("ascii")


def write_multiqc(prefix, sample, source_model, organ, count, image=None):
    name = f"{sample} · {source_model} · CytoFormer ({organ})"
    if image is None:
        html = "<p>No instances detected.</p>"
    else:
        alt = escape(f"CytoFormer cell types for {name}", quote=True)
        html = (
            '<div class="mqc-custom-content-image">'
            f'<img src="data:image/jpeg;base64,{image}" alt="{alt}" /></div>'
            f"<p>{count:,} classified instances slide-wide.</p>"
        )
    payload = {
        "id": f"{prefix}.cytoformer",
        "parent_id": "histopathology_cell_types",
        "parent_name": "Cell type segmentation",
        "parent_description": "Instance masks colored by predicted cell type on crops of the original slides.",
        "section_name": name,
        "plot_type": "image" if image is not None else "html",
        "data": html,
    }
    Path(f"{prefix}.cytoformer_mqc.json").write_text(json.dumps(payload), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=16)
    options = parser.parse_args(shlex.split(${groovy.json.JsonOutput.toJson(task.ext.args ?: '')}))
    if options.batch_size < 1:
        parser.error("--batch-size must be positive")

    prefix = "${prefix}"
    organ = "${organ}"
    segmentation = read_zarr("${instances}")
    source_model = segmentation.attrs["segmentation"]["model"]
    shape_key = "nuclei" if "nuclei" in segmentation.shapes else "cells"
    instances = segmentation.shapes.get(shape_key)
    count = len(instances) if instances is not None else 0
    device = torch.device("${task.accelerator ? 'cuda' : 'cpu'}")
    image = None

    if count:
        if device.type == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CytoFormer requested CUDA, but no CUDA device is visible.")
        torch.set_num_threads(int("${task.cpus}"))
        wsi = open_wsi(read_zarr("original.zarr"), image_key="wsi")
        classes = classify(instances, wsi, organ, device, options.batch_size)
        image = plot_classes(wsi, instances, classes)

    segmentation.attrs["cell_typing"] = {
        "model": "cytoformer",
        "organ": organ,
        "device": device.type,
        "checkpoint_revision": MODEL_REVISION,
        "instance_count": count,
    }
    segmentation.write(f"{prefix}.cytoformer.zarr")
    write_multiqc(prefix, "${meta.id}", source_model, organ, count, image)
    versions = {
        "CytoFormer": MODEL_REVISION,
        **{
            name: version(name)
            for name in ("lazyslide", "spatialdata", "timm", "torch", "wsidata")
        },
    }
    Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))


if __name__ == "__main__":
    main()
