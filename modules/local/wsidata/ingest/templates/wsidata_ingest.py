#!/usr/bin/env python3

import json
import math
import platform
import shutil
from importlib.metadata import version
from pathlib import Path

import dask
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from wsidata import open_wsi


def add_pyramid(wsi):
    height, width = wsi.properties.shape
    factors = []
    while max(height, width) > 256 and min(height, width) >= 2:
        factors.append(2)
        height //= 2
        width //= 2
    if not factors:
        return wsi

    image = Image2DModel.parse(
        wsi.images["wsi"]["scale0"].image,
        scale_factors=factors,
    )
    return open_wsi(SpatialData(images={"wsi": image}, attrs=wsi.attrs), image_key="wsi")


(source,) = Path("source").iterdir()
prefix = "${prefix}"
mpp_arg = "${mpp in [null, []] ? '' : mpp.toString().toBigDecimal()}"
scene_arg = "${scene in [null, []] ? '' : scene.toString().toInteger()}"

suffix = source.suffix.lower()
if suffix == ".czi":
    reader, reader_package = "pylibczi", "pylibCZIrw"
elif suffix in {".tif", ".tiff"}:
    reader, reader_package = "tiffslide", "tiffslide"
elif suffix in {".jpg", ".jpeg"}:
    if not mpp_arg:
        raise ValueError("Supply verified MPP for JPEG slides.")
    import jpype

    # Use the Bio-Formats JAR installed by bftools for offline jobs.
    jar = Path(shutil.which("bfconvert")).with_name("bioformats_package.jar")
    jpype.startJVM("-Djava.awt.headless=true", classpath=[str(jar)])
    reader, reader_package = "bioformats", "scyjava"
else:
    raise ValueError(f"Unsupported slide format: {suffix}")

wsi = open_wsi(
    source,
    store=None,
    reader=reader,
    scene=int(scene_arg) if reader == "pylibczi" and scene_arg else None,
    attach_images=True,
    save_images=True,
)
mpp = float(mpp_arg) if mpp_arg else wsi.properties.mpp
if mpp is None or not math.isfinite(mpp) or mpp <= 0:
    raise ValueError("Supply a finite positive MPP for slides without calibrated metadata.")
wsi.set_mpp(mpp)

scene = None
origin = [0, 0]
if reader == "pylibczi":
    if not scene_arg and wsi.n_scenes > 1:
        raise ValueError(f"Specify a scene for this CZI. Available scenes: {list(range(wsi.n_scenes))}.")
    dimensions = wsi.reader.reader.total_bounding_box
    if any(stop - start != 1 for axis, (start, stop) in dimensions.items() if axis not in {"X", "Y"}):
        raise ValueError("Expected a single-plane brightfield CZI.")
    scene = wsi.scene
    raw = wsi.reader.raw_properties
    origin = [raw["origin_x"], raw["origin_y"]]

if wsi.properties.n_level == 1:
    wsi = add_pyramid(wsi)

report = {
    "source": source.name,
    "format": "wsidata-zarr",
    "reader": reader,
    "scene": scene,
    "source_origin_xy": origin,
    "mpp": mpp,
    "mpp_source": "override" if mpp_arg else "slide_metadata",
    "level_shapes_yx": wsi.properties.level_shape,
}
wsi.attrs["ingestion"] = report
with dask.config.set(scheduler="synchronous"):
    wsi.write(f"{prefix}.zarr", overwrite=False)

versions = {"python": platform.python_version()}
versions.update({name: version(name) for name in ("wsidata", "spatialdata", reader_package)})
if reader == "bioformats":
    versions["bioformats"] = str(jpype.JClass("loci.formats.FormatTools").VERSION)
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}, indent=2) + chr(10))
