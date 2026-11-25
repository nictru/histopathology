process LAZYSLIDE_TISSUEREGIONS {
    tag "${meta.id} (${tile_model})"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container {
        def arm64 = (task.arch ?: System.getProperty('os.arch')) in ['aarch64', 'arm64', 'linux/aarch64', 'linux/arm64', 'linux/arm64/v8']
        if (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container) {
            return arm64
                ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/63aff325d57df75cbe3c7b570e28942a017ff89fa57dd8d69617b12d302518c0/data'
                : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3d3e9582ad18f2daeb8167c469e523a910d770398ae70954a8b46dc5d5094d8b/data'
        }
        return arm64
            ? 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:29a4d91eecf93c47'
            : 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:f7dc31a9e9082220'
    }

    input:
    tuple val(meta), path(domains, stageAs: 'domains.zarr'), path(slide, stageAs: 'slide.zarr'), path(tissue, stageAs: 'tissue.zarr')
    val tile_model

    output:
    tuple val(meta), path("${prefix}.${tile_model}.tissue_regions.zarr"), emit: regions
    path "${prefix}.${tile_model}.tissue_regions_mqc.json", emit: multiqc_files
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/) || !(tile_model ==~ /[a-z0-9][a-z0-9_.-]*/)) {
        error('Invalid tissue region prefix or tile encoder name.')
    }
    template('lazyslide_tissueregions.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir '${prefix}.${tile_model}.tissue_regions.zarr'
    echo '{}' > '${prefix}.${tile_model}.tissue_regions_mqc.json'
    echo '{}' > versions.yml
    """
}
