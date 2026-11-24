process LAZYSLIDE_CELLS {
    tag "${meta.id} (${model})"
    label 'process_medium'
    label 'process_gpu'

    conda {
        if (task.accelerator && (workflow.containerEngine || workflow.wave.enabled)) {
            return null
        }
        "${moduleDir}/${task.accelerator ? 'environment.cuda.yml' : 'environment.yml'}"
    }
    container "${task.accelerator
        ? ((task.arch ?: System.getProperty('os.arch')) in ['aarch64', 'arm64', 'linux/aarch64', 'linux/arm64', 'linux/arm64/v8']
            ? 'community.wave.seqera.io/library/pip_python_cellpose_dask_pruned:d21e545f93f33ca3'
            : 'community.wave.seqera.io/library/pip_python_cellpose_dask_pruned:922892ee7666738d')
        : workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/303a9fde68906168a9d87747014aad1ffcd8275d90de3449169d88db617e478f/data'
            : 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:f370c42ba94e2284'}"

    input:
    tuple val(meta), path(slide, stageAs: 'source.zarr'), path(tiles, stageAs: 'tiles.zarr')
    val model
    path checkpoint, stageAs: 'checkpoint/*', arity: '1'

    output:
    tuple val(meta), path("${prefix}.${model}.cells.zarr"), emit: cells
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    if (!(model in ['instanseg', 'cellpose', 'nulite', 'histoplus'])) {
        error("Unsupported segmentation model: ${model}")
    }
    template('lazyslide_cells.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir '${prefix}.${model}.cells.zarr'
    echo '{}' > versions.yml
    """
}
