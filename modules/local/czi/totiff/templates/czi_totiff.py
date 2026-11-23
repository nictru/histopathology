#!/usr/bin/env python3

import json
import math
import os
import platform
import tempfile
import xml.etree.ElementTree as ET
from importlib.metadata import version
from pathlib import Path

import numpy as np

os.environ["VIPS_CONCURRENCY"] = "${task.cpus}"

import openslide
import pyvips
import tifffile
from pylibCZIrw import czi

TILE = 512


def select_scene(document, requested_scene):
    scenes = document.scenes_bounding_rectangle_no_pyramid
    if requested_scene is not None:
        scene = requested_scene
        if scene not in scenes:
            raise ValueError(
                f"CZI scene {scene} is unavailable. Available scenes: {sorted(scenes)}."
            )
    elif len(scenes) > 1:
        raise ValueError(
            f"Specify a scene for this CZI. Available scenes: {sorted(scenes)}."
        )
    else:
        scene = next(iter(scenes), None)
    bounds = (
        scenes[scene]
        if scene is not None
        else document.total_bounding_rectangle_no_pyramid
    )
    return scene, bounds


def select_plane(document):
    dimensions = {
        key: limits
        for key, limits in document.total_bounding_box.items()
        if key not in ("X", "Y")
    }
    if any(stop - start != 1 for start, stop in dimensions.values()):
        raise ValueError(f"Expected a single-plane brightfield CZI: {dimensions}")
    plane = {key: limits[0] for key, limits in dimensions.items()}
    if document.get_channel_pixel_type(plane.get("C", 0)).lower() != "bgr24":
        raise ValueError("Expected a CZI with native 8-bit BGR pixels.")
    return plane


def read_mpp(document):
    root = ET.fromstring(document.raw_metadata)
    calibration = {}
    for distance in root.findall(".//{*}Scaling/{*}Items/{*}Distance"):
        axis = distance.get("Id")
        value = distance.findtext("{*}Value")
        if axis in {"X", "Y"} and value:
            calibration[axis] = float(value) * 1e6
    if calibration.keys() != {"X", "Y"}:
        raise ValueError("CZI has no X/Y calibration. Supply MPP explicitly.")
    if not math.isclose(calibration["X"], calibration["Y"], rel_tol=1e-3):
        raise ValueError(f"CZI pixel sizes differ by axis: {calibration}")
    return sum(calibration.values()) / 2


def read_rgb(document, scene, plane, roi):
    raw = document.read(
        roi=roi, plane=plane, scene=scene, background_pixel=(1.0, 1.0, 1.0)
    )
    _, _, width, height = roi
    if raw.dtype != np.uint8 or raw.shape != (height, width, 3):
        raise ValueError(f"Unexpected CZI window at {roi}: {raw.dtype} {raw.shape}")
    return raw[:, :, ::-1]


def tiles(document, scene, plane, bounds):
    x0, y0, width, height = bounds
    for y in range(0, height, TILE):
        for x in range(0, width, TILE):
            cols, rows = min(TILE, width - x), min(TILE, height - y)
            roi = (x0 + x, y0 + y, cols, rows)
            tile = np.full((TILE, TILE, 3), 255, dtype=np.uint8)
            tile[:rows, :cols] = read_rgb(document, scene, plane, roi)
            yield tile


def check_output(output, document, scene, plane, bounds, mpp):
    x0, y0, width, height = bounds
    with openslide.OpenSlide(output) as wsi:
        if wsi.dimensions != (width, height):
            raise ValueError("Converted slide dimensions differ from the source.")
        for axis in (openslide.PROPERTY_NAME_MPP_X, openslide.PROPERTY_NAME_MPP_Y):
            reported = wsi.properties.get(axis)
            if reported is None or not math.isclose(float(reported), mpp, rel_tol=1e-5):
                raise ValueError(f"Expected {mpp} MPP, read {reported} for {axis}.")
        size = min(64, width, height)
        positions = (
            (0, 0),
            ((width - size) // 2, (height - size) // 2),
            (width - size, height - size),
        )
        for x, y in positions:
            expected = read_rgb(document, scene, plane, (x0 + x, y0 + y, size, size))
            observed = wsi.read_region((x, y), 0, (size, size)).convert("RGB")
            if not np.array_equal(observed, expected):
                raise ValueError(f"Level-0 pixels changed at ({x}, {y}).")
        levels = [[rows, cols] for cols, rows in wsi.level_dimensions]
        if max(width, height) > 256 and len(levels) < 2:
            raise ValueError("The reader could not access the generated pyramid.")
        if len(levels) > 1:
            wsi.read_region((0, 0), 1, (1, 1))
        return levels


(source,) = Path("source").iterdir()
prefix = "${prefix}"
mpp_override = "${mpp == null || mpp == '' ? '' : mpp.toString().toBigDecimal()}"
scene_arg = "${scene == null || scene == '' ? '' : scene.toString().toInteger()}"
requested_scene = int(scene_arg) if scene_arg else None
process = "${task.process}"

if source.suffix.lower() != ".czi":
    raise ValueError("Expected a CZI slide.")
output = Path(f"{prefix}.tif")
cache = czi.CacheOptions(max_memory_usage=256 * 1024**2)
with czi.open_czi(str(source), cache_options=cache) as document:
    scene, bounds = select_scene(document, requested_scene)
    x0, y0, width, height = bounds
    if width <= 0 or height <= 0:
        raise ValueError(f"Empty CZI canvas: {bounds}")
    plane = select_plane(document)
    mpp = float(mpp_override) if mpp_override else read_mpp(document)
    if not math.isfinite(mpp) or mpp <= 0:
        raise ValueError("MPP must be a finite positive value.")

    with tempfile.TemporaryDirectory(dir=".") as scratch:
        base = Path(scratch) / "level0.tif"
        tifffile.imwrite(
            base,
            tiles(document, scene, plane, bounds),
            shape=(height, width, 3),
            dtype=np.uint8,
            bigtiff=True,
            tile=(TILE, TILE),
            photometric="rgb",
            compression="deflate",
            maxworkers=int("${task.cpus}"),
            metadata=None,
        )
        image = pyvips.Image.new_from_file(str(base), access="sequential")
        image.tiffsave(
            str(output),
            tile=True,
            tile_width=256,
            tile_height=256,
            pyramid=True,
            bigtiff=True,
            compression="deflate",
            resunit="cm",
            xres=1000 / mpp,
            yres=1000 / mpp,
        )
    levels = check_output(output, document, scene, plane, bounds, mpp)

report = {
    "source": source.name,
    "format": "pyramidal-bigtiff",
    "reader": "openslide",
    "scene": scene,
    "source_origin_xy": [x0, y0],
    "background_rgb": [255, 255, 255],
    "mpp": mpp,
    "mpp_source": "override" if mpp_override else "czi_scaling",
    "level_shapes_yx": levels,
}
Path(f"{prefix}.ingest.json").write_text(json.dumps(report, indent=2) + chr(10))
versions = {
    "python": platform.python_version(),
    "libvips": ".".join(str(pyvips.version(i)) for i in range(3)),
    **{
        name: version(name)
        for name in (
            "pylibCZIrw",
            "pyvips",
            "tifffile",
            "openslide-bin",
            "openslide-python",
        )
    },
}
Path("versions.yml").write_text(json.dumps({process: versions}, indent=2) + chr(10))
