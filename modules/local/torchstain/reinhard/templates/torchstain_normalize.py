#!/usr/bin/env python3

import json
import os
import platform
from importlib.metadata import version
from pathlib import Path

os.environ["VIPS_CONCURRENCY"] = "1"

import dask
import numpy as np
import pyvips
import torchstain
import xarray as xr
from spatialdata import SpatialData, read_zarr
from wsidata import open_wsi

BRIGHTNESS_PERCENTILE = 95
LIGHT_INTENSITY = 255
ALPHA = 1
BETA = 0.15
CALIBRATION_GRID = 8
CALIBRATION_PATCH_SIZE = 64


def sample_patches(image):
    height, width = image.shape[1:]
    size = min(CALIBRATION_PATCH_SIZE, width, height)
    xs = np.unique(np.linspace(0, width - size, CALIBRATION_GRID, dtype=int))
    ys = np.unique(np.linspace(0, height - size, CALIBRATION_GRID, dtype=int))
    patches = [
        image[:, y : y + size, x : x + size].transpose(1, 2, 0)
        for y in ys
        for x in xs
    ]
    return np.concatenate(dask.compute(*patches)), len(patches)


def fit_macenko(patches):
    brightness = float(np.percentile(patches, BRIGHTNESS_PERCENTILE))
    scale = LIGHT_INTENSITY / max(brightness, 1.0)
    standardized = np.clip(patches.astype(np.float32) * scale, 0, 255).astype(np.uint8)
    tissue_max = LIGHT_INTENSITY * np.exp(-BETA) - 1
    if np.count_nonzero(standardized.max(axis=-1) <= tissue_max) < 3:
        return scale, None, "insufficient_tissue"

    source = torchstain.normalizers.MacenkoNormalizer(backend="numpy")
    target = torchstain.normalizers.MacenkoNormalizer(backend="numpy")
    try:
        with np.errstate(divide="raise", invalid="raise"):
            source.fit(standardized, Io=LIGHT_INTENSITY, alpha=ALPHA, beta=BETA)
            if np.any(source.maxCRef <= 0) or np.linalg.matrix_rank(source.HERef) < 2:
                return scale, None, "stain_estimation_failed"
            concentrations = np.diag(target.maxCRef / source.maxCRef)
            transform = target.HERef @ concentrations @ np.linalg.pinv(source.HERef)
    except (FloatingPointError, np.linalg.LinAlgError):
        return scale, None, "stain_estimation_failed"
    return scale, transform, None



def reinhard_parameters(query, reference, modified):
    # Use the same query statistics for every tile in the libvips pipeline.
    query_means = np.asarray(query["lab_means"])
    query_stds = np.asarray(query["lab_stds"])
    reference_means = np.asarray(reference["lab_means"])
    reference_stds = np.asarray(reference["lab_stds"])
    if modified:
        ref_lightness_std = reference_stds[0]
        query_lightness_std = query_stds[0]
        contrast = (
            1 - query_lightness_std / ref_lightness_std
            if query_lightness_std < ref_lightness_std
            else 0.05
        )
        gains = np.array([1 + contrast, 1, 1])
        centers = np.array([query_means[0], *reference_means[1:]])
    else:
        # A constant query channel maps to the reference mean.
        gains = np.divide(
            reference_stds,
            query_stds,
            out=np.zeros_like(query_stds),
            where=query_stds > 1e-6,
        )
        centers = reference_means
    return gains.tolist(), (centers - query_means * gains).tolist()



def normalize_block(block, method, parameters):
    rgb = np.ascontiguousarray(block.transpose(1, 2, 0))
    height, width, _ = rgb.shape
    image = pyvips.Image.new_from_memory(rgb.data, width, height, 3, "uchar")
    image = image.copy(interpretation="srgb")
    if method == "macenko":
        scale, transform = parameters
        normalized = image.linear(scale, 0).cast("uchar")
        if transform is not None:
            matrix = pyvips.Image.new_from_array(transform.tolist())
            density = -((normalized.cast("float") + 1) / LIGHT_INTENSITY).log()
            normalized = ((-density.recomb(matrix)).exp() * LIGHT_INTENSITY).cast("uchar")
    else:
        gains, offsets = parameters
        normalized = image.colourspace("lab").linear(gains, offsets).colourspace("srgb")
    return normalized.numpy().transpose(2, 0, 1)


def map_image(image, method, parameters):
    data = image.data.rechunk({0: 3}).map_blocks(
        normalize_block,
        method,
        parameters,
        dtype=np.uint8,
        meta=np.empty((0, 0, 0), dtype=np.uint8),
    )
    return image.copy(data=data)


def normalize_pyramid(image, method, parameters):
    if isinstance(image, xr.DataTree):
        def map_level(dataset):
            if "image" in dataset:
                return dataset.assign(image=map_image(dataset.image, method, parameters))
            return dataset

        return xr.map_over_datasets(map_level, image)
    return map_image(image, method, parameters)


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + chr(10))


prefix = "${prefix}"
sample_id = "${sample_id}"
method = "${normalization_method}"
fit_reference = json.loads("${fit_reference}")
(source,) = Path("source").iterdir()
wsi = open_wsi(read_zarr(source), image_key="wsi")
image = wsi.images["wsi"]
base = image["scale0"].image if isinstance(image, xr.DataTree) else image
if base.dtype != np.uint8 or base.shape[0] != 3:
    raise ValueError("Expected an 8-bit RGB WSIData image.")

with dask.config.set(scheduler="threads", num_workers=int("${task.cpus}")):
    patches, patch_count = sample_patches(base.data)
    report = {
        "source": source.name,
        "mpp": wsi.properties.mpp,
        "level_shapes_yx": wsi.properties.level_shape,
        "method": method,
        "backend": "numpy/libvips",
        "calibration_patch_count": patch_count,
    }
    if method == "macenko":
        scale, transform, fallback = fit_macenko(patches)
        parameters = (scale, transform)
        report.update(
            reference="torchstain_default",
            brightness_percentile=BRIGHTNESS_PERCENTILE,
            brightness_scale=scale,
            light_intensity=LIGHT_INTENSITY,
            alpha=ALPHA,
            beta=BETA,
            fallback=fallback,
        )
    else:
        normalizer = torchstain.normalizers.ReinhardNormalizer(backend="numpy")
        normalizer.fit(patches)
        statistics = {
            "lab_means": normalizer.target_means.tolist(),
            "lab_stds": normalizer.target_stds.tolist(),
            "calibration_patch_count": patch_count,
        }
        if fit_reference:
            write_json(f"{prefix}.reinhard_reference.json", {"sample": sample_id, "method": method, **statistics})
        else:
            reference = json.loads(Path("reference.json").read_text())
            parameters = reinhard_parameters(statistics, reference, method == "reinhardmodified")
            report.update(sample=sample_id, reference_sample=reference["sample"], **statistics)

    if not fit_reference:
        normalized = normalize_pyramid(image, method, parameters)
        output = open_wsi(
            SpatialData(images={"wsi": normalized}, attrs=wsi.attrs),
            image_key="wsi",
        )
        output.attrs["normalization"] = report
        output.write(f"{prefix}.normalized.zarr", overwrite=False)
        write_json(f"{prefix}.normalization.json", report)

versions = {
    "python": platform.python_version(),
    "libvips": ".".join(str(pyvips.version(i)) for i in range(3)),
    **{name: version(name) for name in ("torchstain", "numpy", "pyvips", "wsidata", "spatialdata")},
}
write_json("versions.yml", {"${task.process}": versions})
