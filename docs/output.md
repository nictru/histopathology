# nf-core/histopathology: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:


- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution



### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The report includes these views for each slide.

- **Tissue masks** show a semi-transparent overlay on the original slide for each selected tissue method. Find the masks at `tissue_masks/<sample>.<method>.tissue.zarr`. With the default Otsu-only run, the filename is `<sample>.tissue.zarr`.
- **Slide QC** overlays GrandQC artifact classes on the original slide.
- **Tile QC** displays spatial score maps for the image metrics, FocusLiteNN, and PathProfilerQC. Each map uses its observed score range and averages scores where tiles overlap.
- **Stain normalization** shows matching original and normalized slide thumbnails when normalization runs.
- **Segmentation masks** show semi-transparent cell or nucleus masks on a crop of the original slide for each selected model. Each crop centers on a detected instance.
- **Cell type segmentation** shows the same crop with instances colored by predicted class when the model returns cell-type labels.
- **Cell-type composition and density** shows a sortable metrics table and stacked bars of class proportions for each slide and model.
- **Cell-type neighborhood enrichment** shows one heatmap per slide and model when neighborhood analysis returns defined scores.
- **CytoFormer previews** appear in **Cell type segmentation** when `--cytoformer_organ` is set. Each preview identifies the source segmenter and organ head.
- **Segmentation agreement** appears when two or more segmentation methods run. Switch metric tabs to compare methods on a fixed 0–1 scale. Each heatmap uses the same hierarchical order, based on instance F1. The dendrogram groups methods with detections, and methods without detections appear after them.

The report embeds thumbnails and detail crops. Use the corresponding Zarr outputs for full-resolution images, masks, geometries, and tile scores.

### CytoFormer cell types

Find the classified instances in `cell_types/<sample>.<segmentation_model>.cytoformer.zarr`. Each store preserves the source `nuclei` or `cells` polygons and adds `cytoformer_cell_type` and `cytoformer_probability` to the shapes layer. The `cell_typing` metadata records the organ, checkpoint revision, device, and instance count.

### Cell-type composition and density

Find one CSV per slide and segmentation model under `cell_types/`. CytoFormer adds a separate CSV with `.cytoformer.cell_type_composition.csv` in the filename. Each CSV contains the total instance count and one row per observed class. Instances without a class appear under **Unclassified**. Proportions use all detected instances as the denominator and appear as fractions in CSV and percentages in MultiQC.

Density divides the instance count by the area of the primary tissue mask in mm². The pipeline counts foreground pixels at full resolution and converts their area using the original slide's pixel size. The CSV records the tissue method and area. Proportion is blank when no instances were detected, and density is blank when the tissue area is zero.

### Cell-type neighborhood enrichment

Find per-pair results in `cell_types/<sample>.<model>.cell_neighborhoods.csv`. CytoFormer results include `.cytoformer` after the segmentation model. Each row records the number of undirected edges for one class pair and Squidpy's enrichment z-score. Positive scores indicate more neighboring pairs than expected under shuffled labels, and negative scores indicate fewer. The MultiQC heatmaps use a fixed color range from −5 to +5. The CSV preserves the full scores.

The graph uses classified instance centroids in microns. The CSV records when a result has fewer than two classes, no edges within the radius, or no variation under permutation. Undefined scores are blank.

### Segmentation agreement

Find per-slide scores in `segmentation/<sample>.segmentation_agreement.csv`. The CSV records both directions of each comparison, self-comparisons, instance counts, and the number of matches.

- **Foreground Dice and IoU** compare the areas covered by each method's polygons on the slide.
- **Instance F1, precision, and recall** count one-to-one matches with polygon IoU of at least 0.5. Matching maximizes the number of eligible pairs first and then their summed IoU. Rows are the reference methods and columns are the comparison methods.
- **Mean matched IoU** averages the IoUs of the accepted pairs.

Undefined scores appear as dashes in the heatmaps and empty fields in the CSV.

### Slide embeddings

Find one AnnData file per slide, tile encoder, and aggregation method in `slide_embeddings/<sample>.<tile_encoder>.<method>.slide_embedding.h5ad`. `X` contains one row for the slide. `uns['slide_embedding']` records the sample ID, tile encoder, method, and number of input tiles. `uns['tile_spec']` records the grid when tiles are available. A slide with no eligible tiles has zero rows and a tile count of zero.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
