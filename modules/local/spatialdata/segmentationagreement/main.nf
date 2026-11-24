process SPATIALDATA_SEGMENTATIONAGREEMENT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container {
        def arm64 = (task.arch ?: System.getProperty('os.arch')) in ['aarch64', 'arm64', 'linux/aarch64', 'linux/arm64', 'linux/arm64/v8']
        if (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container) {
            return arm64
                ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/47717b7764e33e509c321074e5d3a06d841722a5fea34494510ee7bc7ef88c5c/data'
                : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2bacc7afedd39aabce34809c9dc500768099afb82dc345d800eb1400206e703c/data'
        }
        arm64
            ? 'community.wave.seqera.io/library/pip_python_matplotlib_numpy_pruned:d120b8172d318a56'
            : 'community.wave.seqera.io/library/pip_python_matplotlib_numpy_pruned:f47cf1913b2b3195'
    }

    input:
    tuple val(meta), path(segmentations, stageAs: 'segmentations/*', arity: '2..*')

    output:
    tuple val(meta), path("${prefix}.segmentation_agreement.csv"), emit: agreement
    path "${prefix}.segmentation_agreement_mqc.json", emit: multiqc_files
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    template('spatialdata_segmentationagreement.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch '${prefix}.segmentation_agreement.csv' '${prefix}.segmentation_agreement_mqc.json'
    echo '{}' > versions.yml
    """
}
