process WSIDATA_INGEST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8cd0517559b7bfc398a13c4aee2b347f91d6ba1753c3471ac068bbef6705acbb/data'
        : 'community.wave.seqera.io/library/bftools_openjdk_pip_python_pruned:e221d8d40eec43ef'}"

    input:
    tuple val(meta), path(slide, stageAs: 'source/*'), val(scene), val(mpp)

    output:
    tuple val(meta), path("${prefix}.zarr"), emit: slide
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    template('wsidata_ingest.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    """
    mkdir '${prefix}.zarr'
    echo '{}' > versions.yml
    """
}
