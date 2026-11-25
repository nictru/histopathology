process LAZYSLIDE_SPATIALDOMAIN {
    tag "${meta.id} (${tile_model})"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'localhost:5000/histopathology-spatialdomain:0.13.0'

    input:
    tuple val(meta), path(embeddings, stageAs: 'embeddings.zarr'), val(tile_model), path(slide, stageAs: 'source.zarr')
    val resolution

    output:
    tuple val(meta), path("${prefix}.${tile_model}.tissue_domains.zarr"), emit: domains
    path "*_mqc.json", emit: multiqc_files
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/) || !(tile_model ==~ /[a-z0-9][a-z0-9_.-]*/)) {
        error('Invalid tissue domain prefix or tile encoder name.')
    }
    domain_resolution = resolution.toString().toDouble()
    if (!(Double.isFinite(domain_resolution) && domain_resolution > 0)) {
        error('Tissue domain resolution must be a finite positive number.')
    }
    template('lazyslide_spatialdomain.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir '${prefix}.${tile_model}.tissue_domains.zarr'
    echo '{}' > '${prefix}.${tile_model}.tissue_domains_mqc.json'
    echo '{}' > versions.yml
    """
}
