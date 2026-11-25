#!/usr/bin/env python3

import csv
import json
from collections import Counter
from importlib.metadata import version
from math import isfinite
from pathlib import Path

from dask.array import count_nonzero, from_zarr
from spatialdata import read_zarr
from zarr import open_group

SECTION = "histopathology_cell_type_composition"
DESCRIPTION = "Instance counts and proportions with densities normalized by the primary tissue area."


def tissue_area_mm2(mask):
    slide = open_group("slide.zarr", mode="r").attrs["ingestion"]
    mpp = float(slide["mpp"])
    if not isfinite(mpp) or mpp <= 0:
        raise ValueError("Cell density requires a positive slide resolution in microns per pixel.")
    label_group = mask["labels/wsi_tissue"]
    metadata = label_group.attrs.get("ome", label_group.attrs)
    level_zero = metadata["multiscales"][0]["datasets"][0]["path"]
    labels = label_group[level_zero]
    if list(labels.shape) != slide["level_shapes_yx"][0]:
        raise ValueError("The primary tissue mask and original slide have different level-0 dimensions.")
    pixels = from_zarr(labels, chunks=labels.shards or labels.chunks)
    tissue_pixels = int(count_nonzero(pixels).compute(num_workers=int("${task.cpus}")))
    return tissue_pixels * (mpp / 1000) ** 2


def class_counts(store, column):
    instances = store.shapes.get("cells", store.shapes.get("nuclei"))
    if instances is None or len(instances) == 0:
        return Counter()
    if column not in instances:
        return Counter({"Unclassified": len(instances)})
    labels = instances[column].astype("string").str.strip()
    return Counter(labels.fillna("Unclassified").replace("", "Unclassified"))


def write_multiqc(prefix, section_id, section_name, plot_type, data, **configuration):
    payload = {
        "id": section_id,
        "parent_id": SECTION,
        "parent_name": "Cell-type composition and density",
        "parent_description": DESCRIPTION,
        "section_name": section_name,
        "plot_type": plot_type,
        "data": data,
        **configuration,
    }
    Path(f"{prefix}_mqc.json").write_text(
        json.dumps(payload, allow_nan=False), encoding="utf-8"
    )


prefix = "${prefix}"
sample = "${meta.id}"
instances = read_zarr("${instances}")
segmenter = instances.attrs["segmentation"]["model"]
typing = instances.attrs.get("cell_typing", {})
organ = typing.get("organ")
classifier = typing.get("model", "segmentation")
label_column = "cytoformer_cell_type" if classifier == "cytoformer" else "class"
model_label = f"{segmenter} + CytoFormer ({organ})" if classifier == "cytoformer" else segmenter

mask = open_group("tissue.zarr", mode="r")
area = tissue_area_mm2(mask)
counts = class_counts(instances, label_column)
total = sum(counts.values())
rows = []
for cell_type, count in [(None, total), *sorted(counts.items())]:
    rows.append(
        {
            "sample": sample,
            "segmentation_model": segmenter,
            "classifier": classifier,
            "organ": organ,
            "tissue_method": mask.attrs["tissue_detection"]["method"],
            "row_type": "total" if cell_type is None else "class",
            "cell_type": cell_type,
            "count": count,
            "proportion": count / total if total else None,
            "cells_per_mm2": count / area if area else None,
            "tissue_area_mm2": area,
        }
    )

with Path(f"{prefix}.cell_type_composition.csv").open("w", encoding="utf-8", newline="") as output:
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

label = f"{sample} · {model_label}"
table = {}
for row in rows:
    row_label = "Total" if row["row_type"] == "total" else f"Class · {row['cell_type']}"
    metrics = {
        "count": row["count"],
        "proportion_pct": 100 * row["proportion"] if row["proportion"] is not None else None,
        "cells_per_mm2": row["cells_per_mm2"],
        "tissue_area_mm2": area,
    }
    table[f"{label} · {row_label}"] = {key: value for key, value in metrics.items() if value is not None}
write_multiqc(
    f"{prefix}.cell_type_table",
    f"{SECTION}_table",
    "Counts, proportions, and density",
    "table",
    table,
    headers={
        "count": {"title": "Count", "format": "{:,.0f}"},
        "proportion_pct": {"title": "Proportion (%)", "min": 0, "max": 100, "format": "{:.2f}"},
        "cells_per_mm2": {"title": "Cells / mm²", "format": "{:,.2f}"},
        "tissue_area_mm2": {"title": "Tissue area (mm²)", "format": "{:,.4f}"},
    },
    pconfig={
        "id": f"{SECTION}_metrics",
        "title": "Cell-type counts, proportions, and density",
        "no_violin": True,
    },
)
if total:
    write_multiqc(
        f"{prefix}.cell_type_proportions",
        f"{SECTION}_proportions",
        "Class proportions",
        "bargraph",
        {label: {name: 100 * count / total for name, count in counts.items()}},
        pconfig={
            "id": f"{SECTION}_class_percentages",
            "title": "Class proportions by slide and model",
            "ylab": "Instances (%)",
            "cpswitch": False,
        },
    )

versions = {name: version(name) for name in ("spatialdata", "dask", "zarr")}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
