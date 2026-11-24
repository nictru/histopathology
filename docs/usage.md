# nf-core/histopathology: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/histopathology/usage](https://nf-co.re/histopathology/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

Pass a CSV file to `--input` with one row per slide. Give each slide a unique sample ID.

| Column | Description |
| ------ | ----------- |
| `sample` | Sample ID used to name the output files. |
| `slide` | Path to a TIFF, Aperio SVS, JPEG, or Zeiss CZI slide. |
| `scene` | Zero-based scene index for a CZI. Leave it blank for TIFF, SVS, JPEG, or a single-scene CZI. |
| `mpp` | Pixel size in microns. Supply it for JPEG. For TIFF, SVS, and CZI, leave it blank to use the source calibration. |

See the [example samplesheet](../assets/samplesheet.csv).

The pipeline uses WSIData to read CZI with pylibCZI, TIFF and SVS with TiffSlide, and JPEG with Bio-Formats. It writes a self-contained `<sample>.zarr` to `<OUTDIR>/slides/`. The store contains the RGB image under the `wsi` key. Its attributes record the pixel calibration and the original CZI scene origin when applicable. The pylibCZI reader fills gaps in a CZI scene with black pixels.

When enabled, normalization reads the ingested Zarr and writes `<sample>.normalized.zarr` to `<OUTDIR>/normalized_slides/`. Reopen either slide store with WSIData. Ingested images, tissue masks, and normalized images use 2048 × 2048 pixel Zarr v3 shards with 256 × 256 pixel chunks. Use a fresh `--outdir` when replacing older, unsharded outputs.

```python
from spatialdata import read_zarr
from wsidata import open_wsi

wsi = open_wsi(read_zarr("sample.normalized.zarr"), image_key="wsi")
```

## Stain normalization

To enable stain normalization, set `--stain_reference` to a sample ID from the samplesheet. Omit it to skip normalization. Reinhard is the default method. Set `--stain_normalization macenko` or `--stain_normalization vahadane` to use another method.

LazySlide detects tissue once per slide using Otsu thresholding on color saturation. The pipeline writes a binary `wsi_tissue` label to `<OUTDIR>/tissue_masks/<sample>.tissue.zarr` for tiling and optional stain normalization.

The pipeline fits the selected reference slide once and uses it to normalize every slide in the samplesheet, including the reference itself. Squidpy uses the precomputed masks to estimate stain statistics from tissue and preserve the background. Each normalized Zarr records the method, reference sample, and fitted statistics in its attributes. Slides with no detected tissue retain their original pixels.

## Tiling and quality control

The pipeline tiles the tissue on each ingested slide and writes `<OUTDIR>/tiles/<sample>.tiles.zarr`. It defaults to 256 × 256 pixels at the slide's native resolution, with no overlap and at most 30% background per tile. Set `--tile_px`, `--tile_mpp`, `--tile_overlap`, or `--tile_background_fraction` to change these values.

The pipeline measures brightness, contrast, sharpness, and saturation on the original slide tiles. It writes `<sample>.tileqc.zarr` to `<OUTDIR>/qc/`. FocusLiteNN, PathProfilerQC, and GrandQC run by default with pinned checkpoints from the LazySlide team on Hugging Face. Set a checkpoint option to use a different `.pt2` file or URL.

| Parameter                     | Output in `<OUTDIR>/qc/`               |
| ----------------------------- | -------------------------------------- |
| `--focuslitenn_checkpoint`    | `<sample>.focuslitenn.zarr`             |
| `--pathprofilerqc_checkpoint`  | `<sample>.pathprofilerqc.zarr`          |
| `--grandqc_checkpoint`        | `<sample>.grandqc.zarr`                |

The pipeline infers the GrandQC variant from the official checkpoint filename. The default checkpoint is for `7x`. To use `5x` or `10x`, supply the corresponding file or URL with `--grandqc_checkpoint`. Set `--grandqc_variant` if your checkpoint has a different filename. The pipeline makes a separate overlapping 512 × 512 grid at 2 MPP for `5x`, 1.5 MPP for `7x`, or 1 MPP for `10x`.

## Nucleus and cell segmentation

Set `--segmentation_model instanseg` to segment nuclei. Segmentation uses the normalized slides when you supply a stain reference and the ingested slides when you omit it. Use a comma-separated list such as `--segmentation_model instanseg,nulite` to run several methods on each slide. The default is `none`.

| Model | Instance polygons | Tile size | Checkpoint parameter |
| --- | --- | --- | --- |
| `instanseg` | Nuclei | 512 × 512 px | `--instanseg_checkpoint` |
| `cellpose` | Objects detected by the Cellpose checkpoint | 256 × 256 px | `--cellpose_checkpoint` |
| `nulite` | Nuclei with cell-type labels | 256 × 256 px | `--nulite_checkpoint` |
| `histoplus` | Nuclei with cell-type labels | 224 × 224 px | `--histoplus_checkpoint` |

Set `--segmentation_mpp` to choose the tile resolution for all selected segmentation models. The default is 0.5 µm/px with 12.5% overlap. Use a resolution compatible with each selected checkpoint and supply calibrated slides scanned at that resolution or finer. Selected models with the same tile size share a grid. The pipeline reuses the regular grid if its size, MPP, and overlap match the model settings. For example, `--tile_px 512 --tile_mpp 0.5 --tile_overlap 0.125` shares the grid between InstanSeg and tile QC at the default segmentation resolution.

InstanSeg, Cellpose, and NuLite use pinned public weights by default. Use the corresponding checkpoint parameter to supply a compatible file or URL.

For HistoPLUS, request access on the [model page](https://huggingface.co/Owkin-Bioptimus/histoplus) and store a Hugging Face read token with access to the repository as the `HF_TOKEN` [Nextflow secret](https://docs.seqera.io/nextflow/secrets). The pipeline downloads the pinned 20× checkpoint once and stages it from the work directory. Use a work directory restricted to authorized users.

For AWS Batch runs launched from Seqera Platform, configure the secret in Platform. For AWS Batch runs launched directly from the CLI, supply an authorized checkpoint file or URI with `--histoplus_checkpoint`.

Find results in `<OUTDIR>/segmentation/<sample>.<model>.cells.zarr`. The `nuclei` or `cells` shapes are in the slide coordinate system. Each store records the model, checkpoint filename, device, tile specification, and instance count. A slide with zero detections has an instance count of zero and no shapes layer.

### GPU acceleration

Add the `gpu` profile to your container profile, for example `-profile docker,gpu`. The segmentation process requests one NVIDIA GPU and uses the CUDA 13.0 environment. Use a compatible NVIDIA driver on the execution host. On AWS Batch, configure a GPU-capable queue for processes with the `process_gpu` label.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/histopathology --input ./samplesheet.csv --stain_reference example_tiff --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/histopathology -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/histopathology
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/histopathology releases page](https://github.com/nf-core/histopathology/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
