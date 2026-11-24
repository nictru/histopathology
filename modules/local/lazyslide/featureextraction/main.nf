process LAZYSLIDE_FEATUREEXTRACTION {
    tag "${meta.id} (${model})"
    label 'process_medium'
    label 'process_gpu'

    conda {
        if (task.accelerator && (workflow.containerEngine || workflow.wave.enabled)) {
            return null
        }
        "${moduleDir}/${task.accelerator ? 'environment.cuda.yml' : 'environment.yml'}"
    }
    container {
        def arm64 = (task.arch ?: System.getProperty('os.arch')) in ['aarch64', 'arm64', 'linux/aarch64', 'linux/arm64', 'linux/arm64/v8']
        if (task.accelerator) {
            return arm64
                ? 'community.wave.seqera.io/library/pip_python_cellpose_dask_pruned:d21e545f93f33ca3'
                : 'community.wave.seqera.io/library/pip_python_cellpose_dask_pruned:922892ee7666738d'
        }
        if (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container) {
            return arm64
                ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f09019fe360d60af6d41470156a66bd53888ba80cc0e25a5a4a770bec8ae8af/data'
                : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/303a9fde68906168a9d87747014aad1ffcd8275d90de3449169d88db617e478f/data'
        }
        arm64
            ? 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:a3e8c447ed0a7cd1'
            : 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:f370c42ba94e2284'
    }

    input:
    tuple val(meta), path(slide, stageAs: 'source.zarr'), path(tiles, stageAs: 'tiles.zarr'), val(model), path(checkpoint, stageAs: 'checkpoint/*', arity: '1')
    val batch_size

    output:
    tuple val(meta), path("${prefix}.${model}.embeddings.zarr"), emit: embeddings
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    if (!(model in ['resnet50', 'uni', 'uni2', 'gigapath-flash'])) {
        error("Unsupported embedding model: ${model}")
    }
    batch = batch_size.toString().toInteger()
    if (batch < 1) {
        error('Expected a positive embedding batch size.')
    }
    template('lazyslide_featureextraction.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir '${prefix}.${model}.embeddings.zarr'
    echo '{}' > versions.yml
    """
}
