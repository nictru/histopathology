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
                ? 'community.wave.seqera.io/library/pip_python_cellpose_dask_pruned:5bd948674751a4d4'
                : 'community.wave.seqera.io/library/pip_python_cellpose_dask_pruned:fcf74c2864653e72'
        }
        if (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container) {
            return arm64
                ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b1c5716b422f5a34775844100883be3c9d53fd9a8550a2f6390aa1e65f705310/data'
                : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f1a8f4bd7afdab5809b500c1cd8630437984e749f6770bfd9d61559c5b8efca/data'
        }
        arm64
            ? 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:64c989016dbcad18'
            : 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:3a9bdbc81931273f'
    }

    input:
    tuple val(meta), path(slide, stageAs: 'source.zarr'), path(tiles, stageAs: 'tiles.zarr'), val(model)
    val batch_size

    output:
    tuple val(meta), path("${prefix}.${model}.embeddings.zarr"), val(model), emit: embeddings
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    if (!(model ==~ /[a-z0-9][a-z0-9_.-]*/)) {
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
