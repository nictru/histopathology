process LAZYSLIDE_TILETISSUES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d548d42355e0d2bf348d26eb274fc87e28bca784f7adb5ced1d9714adbc7c2a9/data'
        : 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:f73e7431e9863fc5'}"

    input:
    tuple val(meta), path(slide, stageAs: 'source.zarr'), path(tissue_mask, stageAs: 'tissue.zarr')
    val tile_px
    val mpp
    val overlap
    val background_fraction

    output:
    tuple val(meta), path("${prefix}.tiles.zarr"), emit: tiles
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    pixels = tile_px.toString().toInteger()
    resolution = mpp.toString().toBigDecimal()
    tile_overlap = overlap.toString().toBigDecimal()
    background = background_fraction.toString().toBigDecimal()
    if (pixels <= 0 || resolution <= 0 || tile_overlap < 0 || tile_overlap >= 1 || background < 0 || background > 1) {
        error('Expected positive tile size and MPP, overlap in [0, 1), and background fraction in [0, 1].')
    }
    template('lazyslide_tiletissues.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    """
    mkdir '${prefix}.tiles.zarr'
    echo '{}' > versions.yml
    """
}
