process SQUIDPY_NORMALIZE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c9b8278f294041a770b0fbc68b12dc380f5224ea88ddb25f5645ec998c523277/data'
        : 'community.wave.seqera.io/library/mahotas_numpy_pip_python_pruned:019535063c082f89'}"

    input:
    tuple val(meta), path(slide, stageAs: 'source.zarr'), path(tissue_mask, stageAs: 'tissue.zarr')
    path fitted_reference, stageAs: 'reference.json'
    val method
    val fit_reference

    output:
    tuple val(meta), path("${prefix}.stain_reference.json"), emit: reference, optional: true
    tuple val(meta), path("${prefix}.normalized.zarr"), emit: slide, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    sample_id = meta.id
    prefix = task.ext.prefix ?: "${sample_id}"
    if (!(sample_id ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/) || !(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe sample ID or output prefix: ${sample_id}, ${prefix}")
    }
    if (!(method in ['macenko', 'reinhard', 'vahadane'])) {
        error("Unsupported stain normalization method: ${method}")
    }
    template('squidpy_normalize.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    def command = fit_reference
        ? "touch '${prefix}.stain_reference.json'"
        : "mkdir '${prefix}.normalized.zarr'"
    """
    ${command}
    echo '{}' > versions.yml
    """
}
