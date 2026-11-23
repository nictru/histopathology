process CZI_TOTIFF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3e2ae61a86eff04763d8258feff476a6f7d781e175db17815d361394f60790d4/data'
        : 'community.wave.seqera.io/library/pip_python_openslide-bin_openslide-python_pruned:e9a5fd377ca7bd3e'}"

    input:
    tuple val(meta), path(slide, stageAs: 'source/*'), val(mpp), val(scene)

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
    template('czi_totiff.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    """
    touch '${prefix}.tif' '${prefix}.ingest.json'
    echo '{}' > versions.yml
    """
}
