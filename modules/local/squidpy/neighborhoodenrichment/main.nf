process SQUIDPY_NEIGHBORHOODENRICHMENT {
    tag "${instances.name}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c9b8278f294041a770b0fbc68b12dc380f5224ea88ddb25f5645ec998c523277/data'
        : 'community.wave.seqera.io/library/mahotas_numpy_pip_python_pruned:019535063c082f89'}"

    input:
    tuple val(meta), path(instances), path(slide, stageAs: 'slide.zarr')
    val radius_um
    val permutations

    output:
    tuple val(meta), path("${prefix}.cell_neighborhoods.csv"), emit: neighborhoods
    path "${prefix}.cell_neighborhoods_mqc.json", optional: true, emit: multiqc_files
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: instances.name.replaceFirst(/(\.cells)?\.zarr$/, '')
    radius = radius_um.toString().toBigDecimal()
    n_perms = permutations.toString().toInteger()
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/) || radius <= 0 || n_perms < 2) {
        error('Expected a safe output prefix, a positive neighborhood radius, and at least two permutations.')
    }
    template('squidpy_neighborhoodenrichment.py')

    stub:
    prefix = task.ext.prefix ?: instances.name.replaceFirst(/(\.cells)?\.zarr$/, '')
    """
    touch '${prefix}.cell_neighborhoods.csv' '${prefix}.cell_neighborhoods_mqc.json'
    echo '{}' > versions.yml
    """
}
