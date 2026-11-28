# nf-core/histopathology: Output

## Output directories

Paths are relative to `--outdir`. Each directory appears when its stage produces a file. If you reuse an output directory, files from earlier runs remain in their original locations.

| Directory                             | Contents                                                                 |
| ------------------------------------- | ------------------------------------------------------------------------ |
| `01_preprocessing/slides/`            | Ingested whole-slide images                                              |
| `01_preprocessing/tissue_masks/`      | Tissue masks for each selected method                                    |
| `01_preprocessing/normalized_slides/` | Stain-normalized slides                                                  |
| `01_preprocessing/tiles/`             | Tissue tile grids                                                        |
| `02_quality_control/slides/`          | GrandQC artifact maps                                                    |
| `02_quality_control/tiles/`           | Image metrics, FocusLiteNN, and PathProfilerQC scores                    |
| `03_segmentation/`                    | Segmented nuclei or cells and agreement scores                           |
| `04_cell_types/`                      | CytoFormer classifications, composition, density, enrichment, and niches |
| `05_embeddings/tiles/`                | Tile embeddings for each encoder                                         |
| `05_embeddings/slides/`               | Aggregated slide embeddings                                              |
| `06_tissue_mapping/domains/`          | Tissue domain labels for each tile encoder                               |
| `06_tissue_mapping/regions/`          | Connected regions from the selected domain map                           |
| `07_cell_features/`                   | Cell or nucleus morphology and regional composition                      |
| `08_multiqc/`                         | MultiQC report, data, and plots                                          |
| `pipeline_info/`                      | Nextflow execution reports and run metadata                              |

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `08_multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The report includes these views for each slide.

- **Tissue masks** show a semi-transparent overlay on the original slide for each selected tissue method. Find the masks at `01_preprocessing/tissue_masks/<sample>.<method>.tissue.zarr`. With the default Otsu-only run, the filename is `<sample>.tissue.zarr`.
- **Slide QC** overlays GrandQC artifact classes on the original slide.
- **Tile QC** displays spatial score maps for the image metrics, FocusLiteNN, and PathProfilerQC. Each map uses its observed score range and averages scores where tiles overlap.
- **Stain normalization** shows matching original and normalized slide thumbnails when normalization runs.
- **Segmentation masks** show semi-transparent cell or nucleus masks on a crop of the original slide for each selected model. Each crop centers on a detected instance.
- **Cell type segmentation** shows the same crop with instances colored by predicted class when the model returns cell-type labels.
- **Cell-type composition and density** shows a sortable metrics table and stacked bars of class proportions for each slide and model.
- **Cell-type neighborhood enrichment** shows one heatmap per slide and model when neighborhood analysis returns defined scores.
- **Spatial cell niches** shows a niche map over the original slide, representative crops, and an interactive heatmap of the cell-type composition around each niche's members.
- **Tissue domains** shows a semi-transparent overlay of tile clusters on the original slide for each selected tile encoder.
- **Tissue regions** overlays the region polygons from the selected domain map on the original slide.
- **Cell and nucleus features** shows interactive morphology box plots and a sortable composition table when tissue regions are available.
- **CytoFormer previews** appear in **Cell type segmentation** when `--cytoformer_organ` is set. Each preview identifies the source segmenter and organ head.
- **Segmentation agreement** appears when two or more segmentation methods run. Switch metric tabs and hover over the interactive heatmaps to compare scores on a fixed 0–1 scale. Each heatmap uses the same hierarchical order based on instance F1. Methods without detections appear last.

The report embeds thumbnails and detail crops. Use the corresponding Zarr outputs for full-resolution images, masks, geometries, and tile scores.

### CytoFormer cell types

Find the classified instances in `04_cell_types/<sample>.<segmentation_model>.cytoformer.zarr`. Each store preserves the source `nuclei` or `cells` polygons and adds `cytoformer_cell_type` and `cytoformer_probability` to the shapes layer. The `cell_typing` metadata records the organ, checkpoint revision, device, and instance count.

### Cell-type composition and density

Find one CSV per slide and segmentation model under `04_cell_types/`. CytoFormer adds a separate CSV with `.cytoformer.cell_type_composition.csv` in the filename. Each CSV contains the total instance count and one row per observed class. Instances without a class appear under **Unclassified**. Proportions use all detected instances as the denominator and appear as fractions in CSV and percentages in MultiQC.

Density divides the instance count by the area of the primary tissue mask in mm². The pipeline counts foreground pixels at full resolution and converts their area using the original slide's pixel size. The CSV records the tissue method and area. Proportion is blank when no instances were detected, and density is blank when the tissue area is zero.

### Cell-type neighborhood enrichment

Find per-pair results in `04_cell_types/<sample>.<model>.cell_neighborhoods.csv`. CytoFormer results include `.cytoformer` after the segmentation model. Each row records the number of undirected edges for one class pair and Squidpy's enrichment z-score. Positive scores indicate more neighboring pairs than expected under shuffled labels, and negative scores indicate fewer. The MultiQC heatmaps use a fixed color range from −5 to +5. The CSV preserves the full scores.

The graph uses classified instance centroids in microns. The CSV records when a result has fewer than two classes, no edges within the radius, or no variation under permutation. Undefined scores are blank.

### Spatial cell niches

Find one SpatialData store per slide and cell-typing method at `04_cell_types/<sample>.<model>.cell_niches.zarr`. CytoFormer filenames include `.cytoformer` after the segmentation model. The store preserves the original polygons. The linked `cell_niches` table contains one row per classified instance. It records the predicted cell type, niche, and number of classified neighbors. Each column in `X` gives the proportion of neighbors of one cell type. Instance centroids in `obsm["spatial"]` use level-0 slide pixels.

Instances without classified neighbors have an unassigned niche and a zero profile. Slides without classified instances have no niche table. MultiQC reports when no niches can be assigned. Niche labels and colors are fitted independently for each slide and cell-typing method. Each representative crop centers on the instance with the closest profile to its niche's mean profile.

### Segmentation agreement

Find per-slide scores in `03_segmentation/<sample>.segmentation_agreement.csv`. The CSV records both directions of each comparison, self-comparisons, instance counts, and the number of matches.

- **Foreground Dice and IoU** compare the areas covered by each method's polygons on the slide.
- **Instance F1, precision, and recall** count one-to-one matches with polygon IoU of at least 0.5. Matching maximizes the number of eligible pairs first and then their summed IoU. Rows are the reference methods and columns are the comparison methods.
- **Mean matched IoU** averages the IoUs of the accepted pairs.

Undefined scores appear as blank cells in the heatmaps and empty fields in the CSV.

### Slide embeddings

Find one AnnData file per slide, tile encoder, and aggregation method in `05_embeddings/slides/<sample>.<tile_encoder>.<method>.slide_embedding.h5ad`. `X` contains one row for the slide. `uns['slide_embedding']` records the sample ID, tile encoder, method, and number of input tiles. `uns['tile_spec']` records the grid when tiles are available. A slide with no eligible tiles has zero rows and a tile count of zero.

### Tissue domain maps

Find one SpatialData store per slide and tile encoder at `06_tissue_mapping/domains/<sample>.<tile_encoder>.tissue_domains.zarr`. Each tile polygon keeps its original `tile_id` and has a `domain` label. The `tissue_domains` metadata records the clustering resolution and the number of tiles and domains.

The pipeline smooths tile embeddings over neighboring grid positions in the same tissue piece and clusters the smoothed features with LazySlide. Domain labels and colors are assigned independently for each slide and encoder. Slides without eligible tiles produce an empty store and a note in MultiQC.

### Tissue regions and cell features

Find the regions in `06_tissue_mapping/regions/<sample>.<tile_encoder>.tissue_regions.zarr`. Each connected polygon has a slide-local `region_id`, a domain label, a tissue area in mm², and the area it shares with other regions. The tissue area used for density excludes overlaps. Region boundaries follow the primary tissue contours used for tiling.

Find morphology in `07_cell_features/<sample>.<model>.cell_features.zarr`. CytoFormer filenames include `.cytoformer` after the segmentation model. The linked `cell_features` table preserves the original instance IDs and records area in µm², perimeter in µm, and unitless circularity and solidity. Invalid polygons retain their IDs and have missing measurements. MultiQC compares the distributions by predicted class and shows the number of valid measurements for each class.

With regions enabled, the table records the region covering each instance's representative point. Points in multiple regions are marked **Ambiguous**, and points outside the regions or from invalid geometries are **Unassigned**. The corresponding `*.region_composition.csv` reports counts, proportions, and instances per mm² by region and predicted class, with a separate breakdown by niche when available. Proportions use the instance count for the region or assignment status. Densities are blank for ambiguous and unassigned instances and for regions with zero usable tissue area.

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
