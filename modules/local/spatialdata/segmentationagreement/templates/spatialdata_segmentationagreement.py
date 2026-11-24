#!/usr/bin/env python3

import base64
import csv
import json
from dataclasses import dataclass
from hashlib import sha256
from html import escape
from importlib.metadata import version
from io import BytesIO
from itertools import combinations
from pathlib import Path

import matplotlib
import numpy as np
import shapely
from matplotlib.figure import Figure
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.spatial.distance import squareform
from spatialdata import read_zarr, transform
from spatialdata.transformations import Identity, get_transformation

matplotlib.use("Agg")

IOU_THRESHOLD = 0.5
QUERY_BATCH_SIZE = 4096
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
    layer = "cells" if model == "cellpose" else "nuclei"
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
        return methods, None

    distances = np.zeros((len(detected), len(detected)))
    for i, j in combinations(range(len(detected)), 2):
        distances[i, j] = distances[j, i] = (
            1 - rows[detected[i], detected[j]]["instance_f1"]
        )
    tree = linkage(squareform(distances), method="average", optimal_ordering=True)
    return [detected[index] for index in leaves_list(tree)] + empty, tree


def heatmap_png(methods, rows, metric, label, tree):
    values = np.array(
        [[rows[left, right][metric] for right in methods] for left in methods],
        dtype=float,
    )
    figure = Figure(figsize=(5.8, 4.8), layout="constrained")
    tiers = (1, 4) if tree is not None else (4,)
    layout = figure.add_gridspec(
        len(tiers), 2, height_ratios=tiers, width_ratios=(20, 1)
    )
    ax = figure.add_subplot(layout[-1, 0])
    if tree is not None:
        branch = figure.add_subplot(layout[0, 0], sharex=ax)
        dendrogram(tree, ax=branch, no_labels=True, link_color_func=lambda _: "#52647d")
        branch.set_axis_off()

    edge = 10 * len(methods)
    centers = np.arange(5, edge, 10)
    palette = matplotlib.colormaps["viridis"].with_extremes(bad="#e2e8f0")
    image = ax.imshow(
        values, cmap=palette, vmin=0, vmax=1, extent=(0, edge, edge, 0), aspect="auto"
    )
    ax.set_xticks(centers, methods, rotation=25, ha="right")
    ax.set_yticks(centers, methods)
    ax.set(xlabel="Comparison method", ylabel="Reference method")
    ax.set_xlim(0, edge)
    for (row, column), value in np.ndenumerate(values):
        missing = np.isnan(value)
        text = "—" if missing else f"{value:.2f}"
        color = "#182235" if missing or value >= 0.55 else "white"
        ax.text(
            centers[column], centers[row], text, ha="center", va="center", color=color
        )
    colorbar = figure.add_subplot(layout[-1, 1])
    figure.colorbar(image, cax=colorbar, ticks=(0, 0.25, 0.5, 0.75, 1), label=label)
    output = BytesIO()
    figure.savefig(output, format="png", dpi=110)
    return base64.b64encode(output.getvalue()).decode("ascii")


def report_tabs(prefix, methods, rows):
    ordered, tree = cluster_methods(methods, rows)
    anchor = "agreement-" + sha256(prefix.encode()).hexdigest()[:16]
    buttons, panels = [], []
    for index, (metric, label) in enumerate(METRICS.items()):
        panel = f"{anchor}-{metric}"
        active = " active" if index == 0 else ""
        selected = "true" if index == 0 else "false"
        buttons.append(
            f'<button type="button" class="nav-link{active}" id="{panel}-tab" '
            f'data-bs-toggle="tab" data-bs-target="#{panel}" role="tab" '
            f'aria-controls="{panel}" aria-selected="{selected}">{label}</button>'
        )
        image = heatmap_png(ordered, rows, metric, label, tree)
        panels.append(
            f'<div class="tab-pane{active}" id="{panel}" role="tabpanel" '
            f'aria-labelledby="{panel}-tab" tabindex="0">'
            '<div class="mqc-custom-content-image">'
            f'<img src="data:image/png;base64,{image}" alt="{escape(label)} by segmentation method" '
            'style="max-width:100%;height:auto"></div></div>'
        )
    return (
        '<div class="nav nav-tabs" role="tablist" aria-label="Segmentation agreement metrics">'
        + "".join(buttons)
        + '</div><div class="tab-content">'
        + "".join(panels)
        + "</div>"
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

payload = {
    "id": f"{prefix}.segmentation_agreement",
    "parent_id": "histopathology_segmentation_agreement",
    "parent_name": "Segmentation agreement",
    "parent_description": (
        "Pairwise comparison of segmentation methods on each slide. Foreground scores "
        "compare the union of detected polygons. Instance scores use one-to-one matches "
        f"at IoU ≥ {IOU_THRESHOLD:g}. Both axes share an average-linkage clustering "
        "based on instance F1. Methods without detections follow the clustered methods. "
        "A dash indicates an undefined score."
    ),
    "section_name": "${meta.id}",
    "plot_type": "html",
    "data": report_tabs(prefix, methods, rows),
}
Path(f"{prefix}.segmentation_agreement_mqc.json").write_text(
    json.dumps(payload), encoding="utf-8"
)

packages = ("matplotlib", "numpy", "scipy", "shapely", "spatialdata")
versions = {package: version(package) for package in packages}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
