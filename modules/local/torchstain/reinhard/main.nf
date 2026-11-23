process TORCHSTAIN_REINHARD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/49a4f2fd6692f69102339bb139addd2dbfdd8ad4578d892b757670da5e5d515d/data'
        : 'community.wave.seqera.io/library/pip_python_pytorch-cpu_numpy_pruned:35643ac99012dd55'}"

    input:
    tuple val(meta), path(slide, stageAs: 'source/*')
    path fitted_reference, stageAs: 'reference.json'
    val modified
    val fit_reference

    output:
    tuple val(meta), path("${prefix}.reinhard_reference.json"), emit: reference, optional: true
    tuple val(meta), path("${prefix}.normalized.zarr"), path("${prefix}.normalization.json"), emit: slide, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    sample_id = meta.id
    prefix = task.ext.prefix ?: "${sample_id}"
    if (!(sample_id ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/) || !(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe sample ID or output prefix: ${sample_id}, ${prefix}")
    }
    normalization_method = modified ? 'reinhardmodified' : 'reinhard'
    template('torchstain_normalize.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    def command = fit_reference
        ? "touch '${prefix}.reinhard_reference.json'"
        : "mkdir '${prefix}.normalized.zarr' && touch '${prefix}.normalization.json'"
    """
    ${command}
    echo '{}' > versions.yml
    """
}
