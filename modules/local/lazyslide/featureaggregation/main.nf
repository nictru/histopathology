process LAZYSLIDE_FEATUREAGGREGATION {
    tag "${meta.id} (${tile_model}, ${method})"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container {
        def arm64 = (task.arch ?: System.getProperty('os.arch')) in ['aarch64', 'arm64', 'linux/aarch64', 'linux/arm64', 'linux/arm64/v8']
        if (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container) {
            return arm64
                ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/04/04d1e609a9cb2fe4c781e2dd782a6423c4332de81d3221cdaafd804283aa6ce8/data'
                : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/ca20b17069148713b0fb27c87afa1f5dfb2ded1cfb54c5da7e3f4d45e7b9b29c/data'
        }
        return arm64
            ? 'community.wave.seqera.io/library/ecos_pip_python_pytorch-cpu_pruned:134a18d954020251'
            : 'community.wave.seqera.io/library/ecos_pip_python_pytorch-cpu_pruned:2d223d1a21bfcdb1'
    }

    input:
    tuple val(meta), path(embeddings, stageAs: 'embeddings.zarr'), path(slide, stageAs: 'source.zarr'), val(tile_model), val(method)

    output:
    tuple val(meta), path("${prefix}.${tile_model}.${method}.slide_embedding.h5ad"), emit: embeddings
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/) || !([tile_model, method].every { name -> name ==~ /[a-z0-9][a-z0-9_.-]*/ })) {
        error('Invalid slide embedding prefix, tile model, or aggregation method.')
    }
    template('lazyslide_featureaggregation.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch '${prefix}.${tile_model}.${method}.slide_embedding.h5ad'
    echo '{}' > versions.yml
    """
}
