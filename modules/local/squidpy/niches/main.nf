process SQUIDPY_NICHES {
    tag "${instances.name}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container {
        def arm64 = (task.arch ?: System.getProperty('os.arch')) in ['aarch64', 'arm64', 'linux/aarch64', 'linux/arm64', 'linux/arm64/v8']
        if (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container) {
            return arm64
                ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a8fa8610235ea4fa141a387087db283d0db6b0e32df9d6ec300078da247a1923/data'
                : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c28c1a6671b80bc112c236c0df9829a3f36855b8288bceab259512f256575821/data'
        }
        return arm64
            ? 'community.wave.seqera.io/library/centrosome_mahotas_numpy_pip_pruned:0b5a49b94cab9177'
            : 'community.wave.seqera.io/library/centrosome_mahotas_numpy_pip_pruned:82986b6673118e70'
    }

    input:
    tuple val(meta), path(instances), path(slide, stageAs: 'slide.zarr')
    val radius_um
    val resolution

    output:
    tuple val(meta), path("${prefix}.cell_niches.zarr"), emit: niches
    path "${prefix}.cell_niches.*_mqc.json", arity: '1..2', emit: multiqc_files
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: instances.name.replaceFirst(/(\.cells)?\.zarr$/, '')
    radius = radius_um.toString().toDouble()
    leiden_resolution = resolution.toString().toDouble()
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/) || !([radius, leiden_resolution].every { value -> Double.isFinite(value) && value > 0 })) {
        error('Spatial cell niches require a safe output prefix, positive radius, and positive resolution.')
    }
    template('squidpy_niches.py')

    stub:
    prefix = task.ext.prefix ?: instances.name.replaceFirst(/(\.cells)?\.zarr$/, '')
    """
    mkdir '${prefix}.cell_niches.zarr'
    echo '{}' > '${prefix}.cell_niches.overview_mqc.json'
    echo '{}' > '${prefix}.cell_niches.composition_mqc.json'
    echo '{}' > versions.yml
    """
}
