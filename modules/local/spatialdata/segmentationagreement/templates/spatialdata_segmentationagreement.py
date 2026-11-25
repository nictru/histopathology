#!/usr/bin/env python3

import csv
import json
from dataclasses import dataclass
from hashlib import sha256
from html import escape
from importlib.metadata import version
from itertools import combinations
from pathlib import Path

import numpy as np
import shapely
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.spatial.distance import squareform
from spatialdata import read_zarr, transform
from spatialdata.transformations import Identity, get_transformation

IOU_THRESHOLD = 0.5
QUERY_BATCH_SIZE = 4096
SECTION = "histopathology_segmentation_agreement"
METRICS = {
    "foreground_dice": "Foreground Dice",
    "foreground_iou": "Foreground IoU",
    "instance_f1": "Instance F1",
    "mean_matched_iou": "Mean matched IoU",
    "precision": "Instance precision",
    "recall": "Instance recall",
}


@dataclass
class Segmentation:
    instances: shapely.STRtree
    areas: np.ndarray
    foreground: shapely.STRtree
    foreground_area: float


def read_segmentation(path):
    store = read_zarr(path)
    model = store.attrs["segmentation"]["model"]
    layer = "cells" if "cells" in store.shapes else "nuclei"
    geometry = np.empty(0, dtype=object)
    if layer in store.shapes:
        shapes = store.shapes[layer]
        mapping = get_transformation(shapes, to_coordinate_system="global")
        if not isinstance(mapping, Identity):
            shapes = transform(shapes, to_coordinate_system="global")
        geometry = shapes.geometry.to_numpy(copy=True)
        invalid = ~shapely.is_valid(geometry)
        geometry[invalid] = shapely.make_valid(
            geometry[invalid], method="structure", keep_collapsed=False
        )
    areas = shapely.area(geometry)
    valid = np.isfinite(areas) & (areas > 0)
    foreground = shapely.disjoint_subset_union_all(geometry[valid])
    return model, Segmentation(
        shapely.STRtree(geometry[valid]),
        areas[valid],
        shapely.STRtree(shapely.get_parts(foreground)),
        float(shapely.area(foreground)),
    )


def intersections(reference, prediction):
    geometries = reference.geometries
    for offset in range(0, len(geometries), QUERY_BATCH_SIZE):
        left, right = prediction.query(
            geometries[offset : offset + QUERY_BATCH_SIZE], predicate="intersects"
        )
        left += offset
        areas = shapely.area(
            shapely.intersection(geometries[left], prediction.geometries[right])
        )
        yield left, right, areas


def matched_instances(reference, prediction):
    candidates = []
    for left, right, shared in intersections(reference.instances, prediction.instances):
        iou = shared / (reference.areas[left] + prediction.areas[right] - shared)
        eligible = iou >= IOU_THRESHOLD
        if eligible.any():
            candidates.append((left[eligible], right[eligible], iou[eligible]))
    if not candidates:
        return 0, None

    left, right, iou = (
        np.concatenate(column) for column in zip(*candidates, strict=True)
    )
    # Only vertices with an eligible edge need to appear in the assignment.
    unique_left, rows = np.unique(left, return_inverse=True)
    unique_right, columns = np.unique(right, return_inverse=True)
    n_left, n_right = len(unique_left), len(unique_right)

    # A private dummy column allows each row to remain unmatched. A real edge's
    # bonus makes an extra match worth more than the largest possible IoU gain.
    dummy_rows = np.arange(n_left)
    bonus = min(n_left, n_right) + 2
    rewards = csr_matrix(
        (
            np.concatenate((bonus + iou, np.ones(n_left))),
            (
                np.concatenate((rows, dummy_rows)),
                np.concatenate((columns, n_right + dummy_rows)),
            ),
        ),
        shape=(n_left, n_right + n_left),
    )
    matched_left, matched_right = min_weight_full_bipartite_matching(
        rewards, maximize=True
    )
    real = matched_right < n_right
    matched_ious = np.asarray(rewards[matched_left[real], matched_right[real]]) - bonus
    return int(real.sum()), float(matched_ious.mean())


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def score_pair(reference, prediction):
    n_reference, n_prediction = len(reference.areas), len(prediction.areas)
    if reference is prediction:
        shared_area = reference.foreground_area
        matches = n_reference
        mean_iou = 1.0 if matches else None
    else:
        shared_area = sum(
            float(shared.sum())
            for _, _, shared in intersections(
                reference.foreground, prediction.foreground
            )
        )
        # Floating-point sums may exceed the smaller area by round-off error.
        shared_area = min(
            shared_area, reference.foreground_area, prediction.foreground_area
        )
        matches, mean_iou = matched_instances(reference, prediction)

    total_area = reference.foreground_area + prediction.foreground_area
    symmetric = {
        "matched_instances": matches,
        "foreground_dice": ratio(2 * shared_area, total_area),
        "foreground_iou": ratio(shared_area, total_area - shared_area),
        "instance_f1": ratio(2 * matches, n_reference + n_prediction),
        "mean_matched_iou": mean_iou,
    }
    return symmetric, n_reference, n_prediction


def pairwise_rows(segmentations):
    methods = sorted(segmentations)
    rows = {}
    for left, right in combinations(methods, 2):
        symmetric, n_left, n_right = score_pair(
            segmentations[left], segmentations[right]
        )
        rows[left, right] = ordered_pair(left, right, n_left, n_right, symmetric)
        rows[right, left] = ordered_pair(right, left, n_right, n_left, symmetric)
    for method in methods:
        scores, count, _ = score_pair(segmentations[method], segmentations[method])
        rows[method, method] = ordered_pair(method, method, count, count, scores)
    return rows


def ordered_pair(reference, comparison, n_reference, n_comparison, scores):
    matches = scores["matched_instances"]
    return {
        "reference": reference,
        "comparison": comparison,
        "iou_threshold": IOU_THRESHOLD,
        "reference_instances": n_reference,
        "comparison_instances": n_comparison,
        **scores,
        "precision": ratio(matches, n_comparison),
        "recall": ratio(matches, n_reference),
    }


def cluster_methods(methods, rows):
    detected = [
        method for method in methods if rows[method, method]["reference_instances"]
    ]
    empty = [method for method in methods if method not in detected]
    if len(detected) < 2:
        return detected + empty

    distances = np.zeros((len(detected), len(detected)))
    for i, j in combinations(range(len(detected)), 2):
        distances[i, j] = distances[j, i] = (
            1 - rows[detected[i], detected[j]]["instance_f1"]
        )
    tree = linkage(squareform(distances), method="average", optimal_ordering=True)
    return [detected[index] for index in leaves_list(tree)] + empty


def metric_tabs(group):
    buttons, panels = [], []
    for index, (metric, label) in enumerate(METRICS.items()):
        panel = f"{group}-{metric}-panel"
        active = " active" if index == 0 else ""
        selected = "true" if index == 0 else "false"
        buttons.append(
            f'<button type="button" class="nav-link{active}" id="{panel}-tab" '
            f'data-bs-toggle="tab" data-bs-target="#{panel}" role="tab" '
            f'aria-controls="{panel}" aria-selected="{selected}">{label}</button>'
        )
        panels.append(
            f'<div class="tab-pane{active}" id="{panel}" role="tabpanel" tabindex="0" '
            f'aria-labelledby="{panel}-tab" data-agreement-section="{group}-{metric}"></div>'
        )

    # MultiQC custom content creates one section per heatmap. Move the rendered
    # plots into tabs when the default template exposes its section wrappers.
    return (
        f'<div id="{group}-tabs" hidden>'
        '<div class="nav nav-tabs" role="tablist" aria-label="Segmentation agreement metrics">'
        + "".join(buttons)
        + '</div><div class="tab-content">'
        + "".join(panels)
        + "</div></div>"
        + f'<script>document.addEventListener("DOMContentLoaded", () => {{'
        f'const tabs = document.getElementById("{group}-tabs");'
        """
        const panels = Array.from(tabs.querySelectorAll("[data-agreement-section]"));
        const sections = panels.map(panel =>
            document.getElementById("mqc-section-wrapper-" + panel.dataset.agreementSection)
        );
        const plots = sections.map(section => section?.querySelector(".mqc-section-plot"));
        if (plots.some(plot => !plot)) return;
        panels.forEach((panel, index) => {
            panel.append(plots[index]);
            if (index) {
                sections[index].hidden = true;
                const link = document.querySelector('.mqc-nav a[href="#' +
                    panel.dataset.agreementSection + '"]');
                if (link) link.closest("li").hidden = true;
            }
        });
        tabs.hidden = false;
        tabs.addEventListener("shown.bs.tab", event => {
            const panel = document.getElementById(event.target.getAttribute("aria-controls"));
            const plot = panel.querySelector(".hc-plot.not_rendered");
            if (plot) window.renderPlot?.(plot.id);
        });
        });</script>
        """
    )


def write_heatmaps(prefix, sample, methods, rows):
    group = "agreement-" + sha256(prefix.encode()).hexdigest()[:16]
    ordered = cluster_methods(methods, rows)
    for index, (metric, label) in enumerate(METRICS.items()):
        section_id = f"{group}-{metric}"
        payload = {
            "id": section_id,
            "parent_id": SECTION,
            "parent_name": "Segmentation agreement",
            "parent_description": (
                "Foreground scores compare the union of detected polygons. Instance scores "
                f"use one-to-one matches at IoU ≥ {IOU_THRESHOLD:g}. Both axes use the "
                "same average-linkage order based on instance F1. Methods without "
                "detections appear last. Blank cells represent undefined scores."
            ),
            "section_name": escape(sample)
            if index == 0
            else f"{escape(sample)} · {label}",
            "section_anchor": section_id,
            "plot_type": "heatmap",
            "xcats": ordered,
            "ycats": ordered,
            "pconfig": {
                "id": f"{section_id}-heatmap",
                "title": f"{sample}: {label}",
                "xlab": "Comparison method",
                "ylab": "Reference method",
                "zlab": label,
                "xcats_samples": False,
                "ycats_samples": False,
                "cluster_rows": False,
                "cluster_cols": False,
                "display_values": True,
                "min": 0,
                "max": 1,
                "colstops": [
                    [0, "#440154"],
                    [0.25, "#3b528b"],
                    [0.5, "#21918c"],
                    [0.75, "#5ec962"],
                    [1, "#fde725"],
                ],
            },
            "data": [
                [rows[left, right][metric] for right in ordered] for left in ordered
            ],
        }
        if index == 0:
            payload["description"] = metric_tabs(group)
        Path(f"{prefix}.segmentation_agreement.{metric}_mqc.json").write_text(
            json.dumps(payload, allow_nan=False), encoding="utf-8"
        )


prefix = "${prefix}"
segmentations = {}
for path in sorted(Path("segmentations").glob("*.cells.zarr")):
    model, segmentation = read_segmentation(path)
    if model in segmentations:
        raise ValueError(f"Duplicate segmentation method {model!r}.")
    segmentations[model] = segmentation
if len(segmentations) < 2:
    raise ValueError("Segmentation agreement requires at least two methods.")

methods = sorted(segmentations)
rows = pairwise_rows(segmentations)
with Path(f"{prefix}.segmentation_agreement.csv").open(
    "w", newline="", encoding="utf-8"
) as output:
    writer = csv.DictWriter(output, fieldnames=list(next(iter(rows.values()))))
    writer.writeheader()
    writer.writerows(rows[key] for key in sorted(rows))

write_heatmaps(prefix, "${meta.id}", methods, rows)

packages = ("numpy", "scipy", "shapely", "spatialdata")
versions = {package: version(package) for package in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
