process QUPATH_PROJECT {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container {
        def arm64 = (task.arch ?: System.getProperty('os.arch')) in ['aarch64', 'arm64', 'linux/aarch64', 'linux/arm64', 'linux/arm64/v8']
        if (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container) {
            return arm64
                ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/71fbea6781753721b09cef48ad0103a4deb29e20062f800cd0dd11e53a6e2f0e/data'
                : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe0834ae193fceed7f2ebb8664f2c8124e8b43959e2d58c4029c697c6c1c6fa9/data'
        }
        return arm64
            ? 'community.wave.seqera.io/library/openjdk_xz:8e5ef81a5fd5a291'
            : 'community.wave.seqera.io/library/openjdk_xz:67cd529261f85418'
    }

    input:
    path slides, stageAs: 'slides/*', arity: '1..*'
    path archive, stageAs: 'QuPath-Linux.tar.xz'

    output:
    path 'qupath', emit: project
    path 'versions.yml', emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('qupath_project.sh')

    stub:
    """
    mkdir -p qupath/images qupath/data
    touch qupath/project.qpproj
    echo '{}' > versions.yml
    """
}
