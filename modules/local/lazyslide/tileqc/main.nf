process LAZYSLIDE_TILEQC {
    tag "${meta.id} (${method})"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d548d42355e0d2bf348d26eb274fc87e28bca784f7adb5ced1d9714adbc7c2a9/data'
        : 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:f73e7431e9863fc5'}"

    input:
    tuple val(meta), path(slide, stageAs: 'source.zarr'), path(tiles, stageAs: 'tiles.zarr')
    val method
    path checkpoint, stageAs: 'checkpoint.pt2'

    output:
    tuple val(meta), path("${prefix}.${method}.zarr"), emit: tiles
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    if (!(method in ['tileqc', 'focuslitenn', 'pathprofilerqc'])) {
        error("Unsupported tile QC method: ${method}")
    }
    if (method == 'tileqc' && checkpoint) {
        error('Pass an empty checkpoint input for tileqc.')
    }
    if (method != 'tileqc' && !checkpoint) {
        error("Supply a checkpoint file for ${method}.")
    }
    if (checkpoint instanceof Collection && checkpoint.size() > 1) {
        error('Supply one checkpoint file per tile QC task.')
    }
    template('lazyslide_tileqc.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    if (!(method in ['tileqc', 'focuslitenn', 'pathprofilerqc'])) {
        error("Unsupported tile QC method: ${method}")
    }
    """
    mkdir '${prefix}.${method}.zarr'
    echo '{}' > versions.yml
    """
}
