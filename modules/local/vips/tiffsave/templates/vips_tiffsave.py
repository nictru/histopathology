#!/usr/bin/env python3

import json
import math
import os
import platform
from importlib.metadata import version
from pathlib import Path

os.environ["VIPS_CONCURRENCY"] = "${task.cpus}"

import openslide
import pyvips
import tifffile


def read_tiff_mpp(path):
    with tifffile.TiffFile(path) as image:
        tags = image.pages[0].tags
        unit = tags.get("ResolutionUnit")
        microns_per_unit = {2: 25_400, 3: 10_000}.get(unit.value if unit else None)
        if microns_per_unit is None:
            raise ValueError(
                "TIFF has no physical resolution unit. Supply MPP explicitly."
            )
        scales = []
        for axis in ("XResolution", "YResolution"):
            tag = tags.get(axis)
            if tag is None:
                raise ValueError(f"Missing TIFF {axis} tag. Supply MPP explicitly.")
            numerator, denominator = tag.value
            if numerator <= 0 or denominator <= 0:
                raise ValueError(f"Invalid TIFF {axis} tag. Supply MPP explicitly.")
            scales.append(microns_per_unit * denominator / numerator)
    if not math.isclose(*scales, rel_tol=1e-3):
        raise ValueError(f"TIFF pixel sizes differ by axis: {scales}")
    return sum(scales) / 2


def check_output(output, image, mpp):
    with openslide.OpenSlide(output) as wsi:
        if wsi.dimensions != (image.width, image.height):
            raise ValueError("Converted slide dimensions differ from the source.")
        for axis in (openslide.PROPERTY_NAME_MPP_X, openslide.PROPERTY_NAME_MPP_Y):
            reported = wsi.properties.get(axis)
            if reported is None or not math.isclose(float(reported), mpp, rel_tol=1e-5):
                raise ValueError(f"Expected {mpp} MPP, read {reported} for {axis}.")
        size = min(64, image.width, image.height)
        positions = (
            (0, 0),
            ((image.width - size) // 2, (image.height - size) // 2),
            (image.width - size, image.height - size),
        )
        for x, y in positions:
            expected = image.crop(x, y, size, size).write_to_memory()
            observed = wsi.read_region((x, y), 0, (size, size)).convert("RGB")
            if observed.tobytes() != expected:
                raise ValueError(f"Level-0 pixels changed at ({x}, {y}).")
        levels = [[rows, cols] for cols, rows in wsi.level_dimensions]
        if max(image.width, image.height) > 256 and len(levels) < 2:
            raise ValueError("The reader could not access the generated pyramid.")
        if len(levels) > 1:
            wsi.read_region((0, 0), 1, (1, 1))
        return levels


(source,) = Path("source").iterdir()
prefix = "${prefix}"
mpp_override = "${mpp == null || mpp == '' ? '' : mpp.toString().toBigDecimal()}"
process = "${task.process}"

suffix = source.suffix.lower()
if suffix not in {".tif", ".tiff", ".jpg", ".jpeg"}:
    raise ValueError("Expected a TIFF or JPEG slide.")
if not mpp_override and suffix in {".jpg", ".jpeg"}:
    raise ValueError("Supply verified MPP for JPEG slides.")
mpp = float(mpp_override) if mpp_override else read_tiff_mpp(source)
if not math.isfinite(mpp) or mpp <= 0:
    raise ValueError("MPP must be a finite positive value.")

image = pyvips.Image.new_from_file(str(source), access="sequential").copy()
if image.format != "uchar" or image.bands != 3 or image.interpretation != "srgb":
    raise ValueError("Expected a slide with 8-bit RGB pixels.")
if image.get_typeof("orientation") and image.get("orientation") != 1:
    raise ValueError("The slide has a non-identity orientation tag.")
if suffix in {".tif", ".tiff"}:
    with tifffile.TiffFile(source) as tiff:
        if len(tiff.pages) != 1 or tiff.pages[0].subifds:
            raise ValueError("Expected a flat, single-image TIFF.")
if image.get_typeof("image-description"):
    image.remove("image-description")

output = Path(f"{prefix}.tif")
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
# The sequential input has already been consumed by tiffsave.
reference = pyvips.Image.new_from_file(str(source), access="random")
levels = check_output(output, reference, mpp)
report = {
    "source": source.name,
    "format": "pyramidal-bigtiff",
    "reader": "openslide",
    "scene": None,
    "source_origin_xy": [0, 0],
    "mpp": mpp,
    "mpp_source": "override" if mpp_override else "tiff_resolution",
    "level_shapes_yx": levels,
}
Path(f"{prefix}.ingest.json").write_text(json.dumps(report, indent=2) + chr(10))
versions = {
    "python": platform.python_version(),
    "libvips": ".".join(str(pyvips.version(i)) for i in range(3)),
    **{
        name: version(name)
        for name in ("pyvips", "tifffile", "openslide-bin", "openslide-python")
    },
}
Path("versions.yml").write_text(json.dumps({process: versions}, indent=2) + chr(10))
