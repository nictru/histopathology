#!/usr/bin/env python3

import json
from collections import defaultdict
from importlib.metadata import version
from math import isfinite
from numbers import Real
from pathlib import Path
from uuid import NAMESPACE_URL, uuid5

import geopandas as gpd
import numpy as np
import pandas as pd
from lazyslide.io import export_annotations
from rasterio.features import shapes as raster_shapes
from rasterio.transform import Affine
from shapely import get_parts, is_empty, is_valid, set_precision, transform, union_all
from shapely.geometry import mapping, shape
from spatialdata import read_zarr
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation
from tifffile import TiffWriter
from wsidata import open_wsi
from xarray import DataTree

TILE_SIZE = 256
OBJECTS_PER_FILE = 25000
QUPATH_ROI_GRID_SIZE = 0.01


def levels(element):
    if isinstance(element, DataTree):
        return [level.image for level in element.children.values()]
    return [element]


def affine_xy(element):
    return get_transformation(element, to_coordinate_system="global").to_affine_matrix(
        input_axes=("x", "y"), output_axes=("x", "y")
    )


def tiff_tiles(image):
    for y in range(0, image.shape[0], TILE_SIZE):
        for x in range(0, image.shape[1], TILE_SIZE * 8):
            strip = image[y : y + TILE_SIZE, x : x + TILE_SIZE * 8].compute()
            for offset in range(0, strip.shape[1], TILE_SIZE):
                yield strip[:, offset : offset + TILE_SIZE]


def write_slide(store, target, image_transform):
    pyramid = levels(store.images["wsi"])
    if not np.allclose(affine_xy(pyramid[0]), image_transform):
        raise ValueError(f"Image coordinate system differs for {target}.")
    mpp = float(store.attrs["ingestion"]["mpp"])
    if not isfinite(mpp) or mpp <= 0:
        raise ValueError(f"Invalid pixel calibration for {target}: {mpp}")
    with TiffWriter(target, bigtiff=True, ome=True) as writer:
        for index, level in enumerate(pyramid):
            image = level.transpose("y", "x", "c").data
            writer.write(
                tiff_tiles(image),
                shape=image.shape,
                dtype=image.dtype,
                photometric="rgb",
                tile=(TILE_SIZE, TILE_SIZE),
                compression="deflate",
                subifds=len(pyramid) - 1 if index == 0 else None,
                subfiletype=1 if index else 0,
                metadata={
                    "axes": "YXS",
                    "PhysicalSizeX": mpp,
                    "PhysicalSizeXUnit": "µm",
                    "PhysicalSizeY": mpp,
                    "PhysicalSizeYUnit": "µm",
                }
                if index == 0
                else None,
            )


def mask_polygons(label):
    # Polygonize pixel edges in bounded blocks, then dissolve block seams.
    pixels = label.transpose("y", "x").data
    chunks = []
    for y in range(0, pixels.shape[0], 2048):
        for x in range(0, pixels.shape[1], 2048):
            foreground = pixels[y : y + 2048, x : x + 2048].compute() != 0
            if foreground.any():
                polygons = [
                    shape(geometry)
                    for geometry, _ in raster_shapes(
                        foreground.astype("uint8"),
                        mask=foreground,
                        transform=Affine.translation(x, y),
                    )
                ]
                chunks.append(union_all(polygons))
    return gpd.GeoDataFrame(geometry=list(get_parts(union_all(chunks))))


def add_table(instances, table, source, *, profiles=False):
    values = table.obs.drop(columns=["region"]).set_index("instance_id")
    if not values.index.is_unique or not values.index.isin(instances.index).all():
        raise ValueError(f"Unknown or duplicate instance IDs in {source}.")
    matrix = table.X
    if matrix is not None:
        if hasattr(matrix, "toarray"):
            matrix = matrix.toarray()
        units = table.var.get("unit", pd.Series("", index=table.var_names))
        columns = [
            f"{name} [{unit}]" if unit and unit != "unitless" else name
            for name, unit in units.items()
        ]
        if profiles:
            columns = [f"Neighborhood fraction: {name}" for name in columns]
        values = values.join(pd.DataFrame(matrix, index=values.index, columns=columns))
    return instances.join(values.add_prefix(f"{source}: "))


def table_stores(stores, kind, model):
    for store in stores[kind]:
        attrs = store.attrs[kind]
        if attrs["segmentation_model"] == model and kind in store.tables:
            yield attrs["classifier"], store.tables[kind]


def valid_as_float32(geometries):
    packed = transform(geometries, lambda xy: xy.astype(np.float32).astype(np.float64))
    return is_valid(packed)


def qupath_geometries(geometries):
    # QuPath stores GeoJSON coordinates as float32 before reducing ROI precision.
    # Snap in double precision where this first conversion would break topology.
    geometries = np.asarray(geometries)
    affected = np.flatnonzero(~valid_as_float32(geometries))
    if not affected.size:
        return geometries
    affected = affected[is_valid(geometries[affected])]
    snapped = set_precision(geometries[affected], QUPATH_ROI_GRID_SIZE)
    usable = ~is_empty(snapped) & valid_as_float32(snapped)
    result = geometries.copy()
    result[affected[usable]] = snapped[usable]
    return result


def feature(uuid_key, object_type, labels, instance_id, geometry, fields):
    properties = {"objectType": object_type}
    classification = fields.pop("classification", None)
    if classification is None:
        properties["classification"] = {"names": labels, "color": None}
    else:
        path_class = json.loads(classification)
        properties["classification"] = {
            "names": [*labels, path_class["name"]],
            "color": path_class["color"],
        }
    measurements = {}
    metadata = {"Source": " / ".join(labels), "Instance ID": str(instance_id)}
    for key, value in fields.items():
        if np.isscalar(value) and not pd.isna(value):
            if (
                isinstance(value, Real)
                and not isinstance(value, (bool, np.bool_))
                and not key.endswith("_id")
            ):
                if isfinite(value):
                    measurements[key] = float(value)
            else:
                metadata[key] = str(value)
    if measurements:
        properties["measurements"] = measurements
    properties["metadata"] = metadata
    return {
        "type": "Feature",
        "id": str(uuid5(NAMESPACE_URL, uuid_key)),
        "geometry": mapping(geometry),
        "properties": properties,
    }


class SlideExport:
    def __init__(self, target, original):
        self.target = target
        self.wsi = open_wsi(original, image_key="wsi")
        self.level0 = levels(original.images["wsi"])[0]
        self.image_transform = affine_xy(self.level0)
        self.layers = []

    def layer(self, label, groups, object_type, frame, transform, class_column=None):
        if frame is None or frame.empty:
            self.layers.append({"layer": label, "objects": 0})
            return
        coordinates = np.linalg.solve(self.image_transform, transform)
        a, b, x, d, e, y = coordinates[:2].ravel()
        frame = frame.set_geometry(frame.geometry.affine_transform([a, b, d, e, x, y]))
        if frame.geometry.isna().any() or frame.geometry.is_empty.any():
            raise ValueError(f"Missing or empty geometry in {label}.")
        if class_column and class_column in frame:
            classes = frame[class_column].astype("string").str.strip()
            frame = frame.assign(
                qupath_class=classes.fillna("").replace("", "Unclassified")
            )
            class_column = "qupath_class"
        else:
            class_column = None
        self.wsi.shapes["qupath_export"] = ShapesModel.parse(frame)
        exported = export_annotations(self.wsi, "qupath_export", classes=class_column)
        del self.wsi.shapes["qupath_export"]
        columns = exported.drop(
            columns=[exported.geometry.name, "qupath_class"], errors="ignore"
        )
        rows = zip(
            exported.index,
            columns.itertuples(index=False, name=None),
            strict=True,
        )
        for start in range(0, len(exported), OBJECTS_PER_FILE):
            geometries = qupath_geometries(
                exported.geometry.array[start : start + OBJECTS_PER_FILE]
            )
            filename = self.target / f"{len(self.layers):03d}-{start:09d}.geojson"
            with filename.open("w", encoding="utf-8") as output:
                output.write('{"type":"FeatureCollection","features":[')
                for offset, geometry in enumerate(geometries):
                    instance_id, values = next(rows)
                    item = feature(
                        f"{self.target.name}/{label}/{start + offset}",
                        object_type,
                        groups,
                        instance_id,
                        geometry,
                        dict(zip(columns.columns, values, strict=True)),
                    )
                    if offset:
                        output.write(",")
                    json.dump(item, output, allow_nan=False, separators=(",", ":"))
                output.write("]}")
        self.layers.append({"layer": label, "objects": len(exported)})

    def shapes(
        self, store, key, label, groups, object_type, class_column=None, frame=None
    ):
        original = store.shapes.get(key)
        if original is None or original.empty:
            self.layer(label, groups, object_type, None, None)
        else:
            self.layer(
                label,
                groups,
                object_type,
                original if frame is None else frame,
                affine_xy(original),
                class_column,
            )


def export_results(export, stores):
    typed = {
        store.attrs["segmentation"]["model"]: store for store in stores["cell_types"]
    }
    for store in stores["cells"]:
        model = store.attrs["segmentation"]["model"]
        key = "cells" if "cells" in store.shapes else "nuclei"
        instances = store.shapes.get(key)
        class_column = "class"
        if instances is not None and not instances.empty:
            if not instances.index.is_unique:
                raise ValueError(f"Duplicate instance IDs in {model}.")
            instances = instances.assign(instance_type=key)
            annotated = typed.get(model)
            if annotated is not None:
                classifications = annotated.shapes[key][
                    ["cytoformer_cell_type", "cytoformer_probability"]
                ]
                if not classifications.index.is_unique or set(
                    classifications.index
                ) != set(instances.index):
                    raise ValueError(f"CytoFormer and {model} instance IDs differ.")
                instances = instances.join(classifications)
                class_column = "cytoformer_cell_type"
            for kind in ("cell_niches", "cell_features"):
                for classifier, table in table_stores(stores, kind, model):
                    instances = add_table(
                        instances,
                        table,
                        f"{classifier} {kind}",
                        profiles=kind == "cell_niches",
                    )
        export.shapes(
            store,
            key,
            model,
            ["Segmentation", model],
            "detection",
            class_column,
            instances,
        )

    for store in stores["tissue"]:
        method = store.attrs["tissue_detection"]["method"]
        label = levels(store.labels["wsi_tissue"])[0]
        export.layer(
            f"Tissue mask {method}",
            ["Tissue mask", method],
            "annotation",
            mask_polygons(label),
            affine_xy(label),
        )
    for store in stores["slide_qc"]:
        export.shapes(store, "artifacts", "GrandQC", ["GrandQC"], "annotation", "class")
    for store in stores["regions"]:
        encoder = store.attrs["tissue_regions"]["tile_model"]
        export.shapes(
            store,
            "tissue_regions",
            f"Tissue regions {encoder}",
            ["Tissue region", encoder],
            "annotation",
            "domain",
        )
    for store in stores["domains"]:
        encoder = store.attrs["tissue_domains"]["tile_model"]
        export.shapes(
            store,
            "tiles",
            f"Tissue domains {encoder}",
            ["Tissue domain", encoder],
            "tile",
            "domain",
        )

    (grid,) = stores["tiles"]
    tiles = grid.shapes.get("tiles")
    if tiles is not None and not tiles.empty:
        tiles = tiles.set_index("tile_id", drop=False)
        for score_store in stores["tile_qc"]:
            scored = score_store.shapes.get("tiles")
            if scored is None or scored.empty:
                continue
            method = score_store.attrs["tile_qc"]["method"]
            scores = scored.set_index("tile_id").drop(
                columns=["geometry", "tissue_id"], errors="ignore"
            )
            if not scores.index.is_unique or not scores.index.isin(tiles.index).all():
                raise ValueError(f"Tile IDs do not match the tiling grid for {method}.")
            tiles = tiles.join(scores.add_prefix(f"{method}: "))
    export.shapes(grid, "tiles", "Tile grid and QC", ["Tile grid"], "tile", frame=tiles)


sources = json.loads(r"""${sources}""")
stores = defaultdict(list)
for source in sources:
    stores[source["kind"]].append(read_zarr(source["path"]))
(original,) = stores["original"]
target = Path("${prefix}")
target.mkdir()
export = SlideExport(target, original)
write_slide(original, target / "original.ome.tif", export.image_transform)
for normalized in stores["normalized"]:
    if levels(normalized.images["wsi"])[0].sizes != export.level0.sizes:
        raise ValueError("Normalized and original slide dimensions differ.")
    write_slide(normalized, target / "normalized.ome.tif", export.image_transform)
export_results(export, stores)
(target / "layers.json").write_text(json.dumps(export.layers, indent=2) + chr(10))
packages = ("lazyslide", "rasterio", "spatialdata", "tifffile", "wsidata")
Path("versions.yml").write_text(
    json.dumps({"${task.process}": {name: version(name) for name in packages}})
    + chr(10)
)
