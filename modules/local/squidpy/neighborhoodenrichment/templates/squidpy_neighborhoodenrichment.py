#!/usr/bin/env python3

import csv
import json
import os
from importlib.metadata import version
from itertools import combinations_with_replacement
from pathlib import Path

os.environ["NUMBA_NUM_THREADS"] = "${task.cpus}"
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(".numba-cache").resolve()))

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialdata import read_zarr
from zarr import open_group

CLUSTER_KEY = "cell_type"
SECTION = "histopathology_cell_neighborhoods"


def classified_instances(store, column):
    instances = store.shapes.get("cells", store.shapes.get("nuclei"))
    if instances is None or column not in instances:
        return pd.Categorical([]), np.empty((0, 2))
    labels = instances[column].astype("string").str.strip()
    keep = labels.notna() & ~labels.str.casefold().isin(["", "background", "unclassified"])
    centers = instances.loc[keep].geometry.centroid
    coordinates = np.column_stack((centers.x, centers.y))
    if not np.isfinite(coordinates).all():
        raise ValueError("Classified instances must have finite centroid coordinates.")
    return pd.Categorical(labels[keep]), coordinates


def enrichment(labels, coordinates, radius_um, permutations):
    from squidpy import gr

    mpp = float(open_group("slide.zarr", mode="r").attrs["ingestion"]["mpp"])
    if not np.isfinite(mpp) or mpp <= 0:
        raise ValueError("Cell neighborhood enrichment requires a positive slide pixel size.")
    adata = AnnData(
        obs=pd.DataFrame({CLUSTER_KEY: labels}, index=[str(i) for i in range(len(labels))])
    )
    adata.obsm["spatial"] = coordinates * mpp
    gr.spatial_neighbors_radius(adata, radius=radius_um, n_jobs=int("${task.cpus}"))
    edge_count = adata.obsp["spatial_connectivities"].nnz // 2
    if not edge_count:
        return edge_count, None

    result = gr.nhood_enrichment(
        adata,
        cluster_key=CLUSTER_KEY,
        n_perms=permutations,
        rng=0,
        n_jobs=int("${task.cpus}"),
        show_progress_bar=False,
        copy=True,
    )
    return edge_count, result


def finite_score(score):
    return float(score) if np.isfinite(score) else None


def write_heatmap(prefix, model_label, axes, scores, radius_um):
    section_id = f"{SECTION}_{prefix}"
    payload = {
        "id": section_id,
        "parent_id": SECTION,
        "parent_name": "Cell-type neighborhood enrichment",
        "parent_description": "Squidpy compares observed neighboring class pairs with shuffled cell-type labels.",
        "section_name": model_label,
        "description": f"Centroids are connected within {radius_um:g} µm. Positive z-scores indicate enrichment and negative scores indicate depletion.",
        "plot_type": "heatmap",
        "xcats": axes,
        "ycats": axes,
        "pconfig": {
            "id": f"{section_id}_zscores",
            "title": "Cell-type neighborhood enrichment",
            "zlab": "Enrichment z-score",
            "min": -5,
            "max": 5,
            "colstops": [[0, "#2166ac"], [0.5, "#f7f7f7"], [1, "#b2182b"]],
        },
        "data": [[finite_score(value) for value in row] for row in scores],
    }
    Path(f"{prefix}.cell_neighborhoods_mqc.json").write_text(
        json.dumps(payload, allow_nan=False), encoding="utf-8"
    )


prefix = "${prefix}"
sample = "${meta.id}"
radius_um = float("${radius}")
permutations = int("${n_perms}")
store = read_zarr("${instances}")
segmenter = store.attrs["segmentation"]["model"]
typing = store.attrs.get("cell_typing", {})
classifier = typing.get("model", "segmentation")
organ = typing.get("organ")
column = "cytoformer_cell_type" if classifier == "cytoformer" else "class"
labels, coordinates = classified_instances(store, column)
classes = labels.categories.tolist()
base = {
    "sample": sample,
    "segmentation_model": segmenter,
    "classifier": classifier,
    "organ": organ,
    "radius_um": radius_um,
    "permutations": permutations,
    "classified_instances": len(labels),
    "graph_edges": 0,
    "status": "insufficient_classes",
}
rows = []
if len(classes) >= 2:
    base["graph_edges"], result = enrichment(labels, coordinates, radius_um, permutations)
    base["status"] = "no_neighbors"
    if result is not None:
        has_scores = np.isfinite(result.zscore).any()
        base["status"] = "ok" if has_scores else "no_permutation_variance"
        for i, j in combinations_with_replacement(range(len(classes)), 2):
            # The undirected graph stores both orientations for within-class edges.
            observed = int(result.counts[i, j]) // (2 if i == j else 1)
            rows.append(
                (classes[i], classes[j], observed, finite_score(result.zscore[i, j]))
            )
        if has_scores:
            model = f"{segmenter} + CytoFormer ({organ})" if classifier == "cytoformer" else segmenter
            write_heatmap(prefix, f"{sample} · {model}", classes, result.zscore, radius_um)

with Path(f"{prefix}.cell_neighborhoods.csv").open("w", encoding="utf-8", newline="") as output:
    writer = csv.writer(output)
    writer.writerow([*base, "class_a", "class_b", "observed_pairs", "enrichment_zscore"])
    for pair in rows or [(None, None, None, None)]:
        writer.writerow([*base.values(), *pair])

versions = {name: version(name) for name in ("squidpy", "spatialdata", "anndata")}
Path("versions.yml").write_text(json.dumps({"${task.process}": versions}) + chr(10))
