process LAZYSLIDE_FINDTISSUES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d548d42355e0d2bf348d26eb274fc87e28bca784f7adb5ced1d9714adbc7c2a9/data'
        : 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:f73e7431e9863fc5'}"

    input:
    tuple val(meta), path(slide, stageAs: 'source.zarr')

    output:
    tuple val(meta), path("${prefix}.tissue.zarr"), emit: mask
    path "${prefix}.tissue_mqc.json", emit: multiqc_files
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    template('lazyslide_findtissues.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    """
    mkdir '${prefix}.tissue.zarr'
    touch '${prefix}.tissue_mqc.json'
    echo '{}' > versions.yml
    """
}
