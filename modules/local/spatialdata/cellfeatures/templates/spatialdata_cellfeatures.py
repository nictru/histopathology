#!/usr/bin/env python3

import csv
import json
from html import escape
from importlib.metadata import version
from math import isfinite
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from anndata import AnnData
from spatialdata import SpatialData, read_zarr, transform
from spatialdata.models import TableModel
from zarr import open_group

MEASUREMENTS = {
    "area": ("Area", "µm²"),
    "perimeter": ("Perimeter", "µm"),
    "circularity": ("Circularity", ""),
    "solidity": ("Solidity", ""),
}


def measure(instances, mpp):
    geometry = transform(instances, to_coordinate_system="global").geometry.to_numpy()
    valid = shapely.is_valid(geometry) & (shapely.area(geometry) > 0)
    polygons = geometry[valid]
    area = shapely.area(polygons)
    perimeter = shapely.length(polygons)
    values = np.full((len(instances), len(MEASUREMENTS)), np.nan)
    values[valid] = np.column_stack(
        (
            area * mpp**2,
            perimeter * mpp,
            4 * np.pi * area / perimeter**2,
            area / shapely.area(shapely.convex_hull(polygons)),
        )
    )
    return (
        values,
        gpd.GeoSeries(shapely.point_on_surface(polygons), index=instances.index[valid]),
        valid,
    )


def region_assignment(obs, points, regions):
    region_ids = pd.Series(pd.NA, index=obs.index, dtype="string")
    status = pd.Series("Unassigned", index=obs.index)
    if regions is not None and len(regions) and len(points):
        matches = gpd.GeoDataFrame(geometry=points).sjoin(
            regions[["region_id", "geometry"]], predicate="covered_by"
        )["region_id"]
        multiple = matches.index.duplicated(keep=False)
        region_ids.loc[matches.index[~multiple]] = matches.loc[~multiple].astype(str)
        status.loc[matches.index[~multiple]] = "Assigned"
        status.loc[matches.index[multiple].unique()] = "Ambiguous"
    obs["tissue_region_id"] = pd.Categorical(region_ids)
    obs["region_assignment"] = pd.Categorical(status)


def add_niches(obs, store, sample, segmenter, classifier, organ):
    info = store.attrs["cell_niches"]
    if (
        info["sample_id"],
        info["segmentation_model"],
        info["classifier"],
        info["organ"],
    ) != (sample, segmenter, classifier, organ):
        raise ValueError(
            "The niche store does not match this slide and cell-typing source."
        )
    table = store.tables.get("cell_niches")
    if table is None:
        obs["niche"] = pd.Categorical(["Unassigned"] * len(obs))
        return
    assignments = table.obs.set_index("instance_id")["niche"]
    if (
        not assignments.index.is_unique
        or not assignments.index.isin(obs["instance_id"]).all()
    ):
        raise ValueError("Niche assignments contain duplicate or unknown instance IDs.")
    obs["niche"] = pd.Categorical(
        obs["instance_id"].map(assignments).astype("string").fillna("Unassigned")
    )


def composition_rows(obs, regions, identity):
    regions_by_id = (
        {}
        if regions is None
        else {
            str(row.region_id): (str(row.domain), float(row.tissue_area_mm2))
            for row in regions.itertuples()
        }
    )
    group_ids = (
        obs["tissue_region_id"]
        .astype("string")
        .fillna(obs["region_assignment"].astype(str))
    )
    totals = obs.groupby(group_ids, observed=True).size()
    types = [column for column in ("cell_type", "niche") if column in obs]
    counts = {
        column: obs.groupby([group_ids, column], observed=True).size()
        for column in types
    }
    for region_id in [
        *regions_by_id,
        *(key for key in ("Ambiguous", "Unassigned") if key in totals),
    ]:
        domain, area_mm2 = regions_by_id.get(region_id, (None, None))
        total = int(totals.get(region_id, 0))
        groups = [("total", "", total)]
        for column, values in counts.items():
            if region_id in values.index.get_level_values(0):
                groups.extend(
                    (column, str(value), int(count))
                    for value, count in values.loc[region_id].items()
                )
        for category, label, count in groups:
            yield {
                **identity,
                "region_id": region_id,
                "domain": domain,
                "tissue_area_mm2": area_mm2,
                "category": category,
                "label": label,
                "count": count,
                "proportion": count / total if total else None,
                "instances_per_mm2": count / area_mm2 if area_mm2 else None,
            }


def write_report(output, kind, name, title, plot_type, data, **configuration):
    section = "histopathology_cell_features"
    payload = {
        "id": f"{output}.{kind}",
        "parent_id": section,
        "parent_name": "Cell and nucleus features",
        "section_name": f"{escape(title)} · {name}",
        "plot_type": plot_type,
        "data": data,
        **configuration,
    }
    Path(f"{output}.{kind}_mqc.json").write_text(
        json.dumps(payload, allow_nan=False), encoding="utf-8"
    )


def morphology_report(output, title, table):
    datasets = []
    labels = []
    if table is not None and np.isfinite(table.X).any():
        cell_types = table.obs.groupby("cell_type", observed=True, sort=True).indices
        for idx, (name, unit) in enumerate(MEASUREMENTS.values()):
            data = {}
            for label, positions in cell_types.items():
                values = table.X[positions, idx]
                values = values[np.isfinite(values)]
                if len(values):
                    data[f"{label} ({len(values):,})"] = np.percentile(
                        values, np.arange(101)
                    ).tolist()
            if data:
                axis = f"{name} ({unit})" if unit else name
                datasets.append(data)
                labels.append({"name": name, "ylab": axis})
    if datasets:
        write_report(
            output,
            "morphology",
            "Morphology",
            title,
            "boxplot",
            datasets,
            description="Distributions use 101 equally spaced percentiles. The number in parentheses is the number of measured instances.",
            pconfig={
                "id": f"{output}.morphology_plot",
                "title": "Per-instance morphology",
                "data_labels": labels,
            },
        )
    else:
        write_report(
            output,
            "morphology",
            "Morphology",
            title,
            "html",
            "<p>No valid instance polygons to measure.</p>",
        )


def write_composition(prefix, rows, title, region_model):
    fields = (
        "sample",
        "segmentation_model",
        "classifier",
        "organ",
        "region_id",
        "domain",
        "tissue_area_mm2",
        "category",
        "label",
        "count",
        "proportion",
        "instances_per_mm2",
    )
    with Path(f"{prefix}.region_composition.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    output = f"{prefix}.cell_features"
    if not rows:
        write_report(
            output,
            "regions",
            "Tissue region composition",
            title,
            "html",
            "<p>No tissue regions or segmented instances available.</p>",
        )
        return
    data = {}
    for row in rows:
        region = row["region_id"]
        label = row["label"] or "Total"
        row_name = f"{region} · {row['category']} · {label}"
        data[row_name] = {
            "domain": row["domain"],
            "count": row["count"],
            "tissue_area_mm2": row["tissue_area_mm2"],
            "instances_per_mm2": row["instances_per_mm2"],
            "proportion": row["proportion"] * 100
            if row["proportion"] is not None
            else None,
        }
        data[row_name] = {
            key: value for key, value in data[row_name].items() if value is not None
        }
    write_report(
        output,
        "regions",
        "Tissue region composition",
        title,
        "table",
        data,
        description=f"Regions use {escape(region_model)}. Densities use tissue area after excluding overlaps.",
        pconfig={"id": f"{output}.regions_table", "title": "Regional composition"},
        headers={
            "domain": {"title": "Domain", "scale": False},
            "count": {"title": "Count", "format": "{:,.0f}"},
            "tissue_area_mm2": {"title": "Tissue area (mm²)", "format": "{:.4f}"},
            "instances_per_mm2": {"title": "Instances/mm²", "format": "{:,.1f}"},
            "proportion": {"title": "Proportion (%)", "format": "{:.1f}"},
        },
    )


prefix = "${prefix}"
sample = "${meta.id}"
source = read_zarr("${instances}")
shape_key = "cells" if "cells" in source.shapes else "nuclei"
instances = source.shapes.get(shape_key)
segmentation = source.attrs["segmentation"]
segmenter = segmentation["model"]
typing = source.attrs.get("cell_typing", {})
classifier = typing.get("model", "segmentation")
organ = typing.get("organ")
label_column = "cytoformer_cell_type" if classifier == "cytoformer" else "class"
instance_type = segmentation.get(
    "instance_type", "nucleus" if shape_key == "nuclei" else "cell"
)
title = f"{sample} · {segmenter}"
if classifier == "cytoformer":
    title += f" + CytoFormer ({organ})"
slide = open_group("slide.zarr", mode="r").attrs["ingestion"]
mpp = float(slide["mpp"])
if not isfinite(mpp) or mpp <= 0:
    raise ValueError("Cell morphology requires a positive slide pixel size.")

region_store = read_zarr("regions.zarr") if Path("regions.zarr").exists() else None
regions = None
region_model = None
if region_store is not None:
    region_info = region_store.attrs["tissue_regions"]
    if region_info["sample_id"] != sample or region_info["mpp"] != mpp:
        raise ValueError("The tissue regions do not match this calibrated slide.")
    region_model = region_info["tile_model"]
    regions = region_store.shapes.get("tissue_regions")
    if regions is not None:
        regions = transform(regions, to_coordinate_system="global")

obs = pd.DataFrame({"cell_type": pd.Categorical([])})
table = None
if instances is not None and len(instances):
    if not instances.index.is_unique or not instances.index.astype(str).is_unique:
        raise ValueError("Cell or nucleus polygons must have unique instance IDs.")
    values, points, valid = measure(instances, mpp)
    labels = (
        instances.get(label_column, pd.Series(pd.NA, index=instances.index))
        .astype("string")
        .str.strip()
    )
    labels = labels.fillna("Unclassified").replace("", "Unclassified")
    obs = pd.DataFrame(
        {
            "instance_id": instances.index.to_numpy(),
            "region": pd.Categorical([shape_key] * len(instances)),
            "instance_type": pd.Categorical([instance_type] * len(instances)),
            "cell_type": pd.Categorical(labels),
            "valid_geometry": valid,
        },
        index=instances.index,
    )
    if region_store is not None:
        region_assignment(obs, points, regions)
    if Path("niches.zarr").exists():
        add_niches(obs, read_zarr("niches.zarr"), sample, segmenter, classifier, organ)
    coordinates = np.full((len(instances), 2), np.nan)
    coordinates[valid] = np.column_stack((points.x, points.y))
    obs.index = obs.index.astype(str)
    table = AnnData(
        values,
        obs=obs,
        var=pd.DataFrame(
            {"unit": [unit or "unitless" for _, unit in MEASUREMENTS.values()]},
            index=list(MEASUREMENTS),
        ),
        obsm={"spatial": coordinates},
    )
    TableModel.parse(
        table, region=shape_key, region_key="region", instance_key="instance_id"
    )
elif region_store is not None:
    region_assignment(obs, gpd.GeoSeries(), regions)

identity = {
    "sample": sample,
    "segmentation_model": segmenter,
    "classifier": classifier,
    "organ": organ,
}
attrs = {
    **identity,
    "mpp": mpp,
    "instance_type": instance_type,
    "tile_model": region_model,
    "instance_count": len(obs),
    "measured_instances": int(obs["valid_geometry"].sum()) if table is not None else 0,
}
SpatialData(
    shapes={shape_key: instances} if instances is not None else {},
    tables={"cell_features": table} if table is not None else {},
    attrs={"cell_features": attrs},
).write(f"{prefix}.cell_features.zarr")
morphology_report(f"{prefix}.cell_features", title, table)
if region_store is not None:
    write_composition(
        prefix, list(composition_rows(obs, regions, identity)), title, region_model
    )

packages = ("anndata", "geopandas", "shapely", "spatialdata")
versions = {name: version(name) for name in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
