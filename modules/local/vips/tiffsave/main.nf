process VIPS_TIFFSAVE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1ea398788b2754d197335271055c3375c47965b1ab85dbc61df852c6b0eefde4/data'
        : 'community.wave.seqera.io/library/pip_python_openslide-bin_openslide-python_pruned:1754049372e890cd'}"

    input:
    tuple val(meta), path(slide, stageAs: 'source/*'), val(mpp)

    output:
    tuple val(meta), path("${prefix}.tif"), path("${prefix}.ingest.json"), emit: slide
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    template('vips_tiffsave.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    """
    touch '${prefix}.tif' '${prefix}.ingest.json' versions.yml
    """
}
