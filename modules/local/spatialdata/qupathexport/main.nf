process SPATIALDATA_QUPATHEXPORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container {
        def arm64 = (task.arch ?: System.getProperty('os.arch')) in ['aarch64', 'arm64', 'linux/aarch64', 'linux/arm64', 'linux/arm64/v8']
        if (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container) {
            return arm64
                ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f8a9e062fd5d029b951e763e93d44239cb1f990def3f85d667da843ad06dbbc3/data'
                : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a6c1d52d60aacb9aecded7b07ee39ce03459350c634f41fde0f4e234216fa247/data'
        }
        return arm64
            ? 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:cd8a22ff9c0c9af5'
            : 'community.wave.seqera.io/library/pip_python_pytorch-cpu_torchvision_pruned:b4bbaeffd6e10901'
    }

    input:
    tuple val(meta), val(kinds), path(stores, stageAs: 'stores/*', arity: '1..*')

    output:
    tuple val(meta), path("${prefix}"), emit: slide
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
        error("Unsafe output prefix: ${prefix}")
    }
    sources = groovy.json.JsonOutput.toJson(
        kinds.indices.collect { index ->
            [kind: kinds[index], path: stores[index].toString()]
        }
    )
    template('spatialdata_qupathexport.py')

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    mkdir '${prefix}'
    touch '${prefix}/original.ome.tif'
    echo '{}' > versions.yml
    """
}
