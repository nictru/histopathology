#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/histopathology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/histopathology
    Website: https://nf-co.re/histopathology
    Slack  : https://nfcore.slack.com/channels/histopathology
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HISTOPATHOLOGY  } from './workflows/histopathology'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_histopathology_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_histopathology_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_HISTOPATHOLOGY {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HISTOPATHOLOGY (
        samplesheet,
        params.stain_normalization,
        params.stain_reference?.toString(),
        params.tissue_segmentation_method,
        params.tile_px,
        params.tile_mpp == null ? '' : params.tile_mpp,
        params.tile_overlap,
        params.tile_background_fraction,
        params.grandqc_variant,
        params.segmentation_model,
        params.segmentation_mpp,
        params.cytoformer_organ?.toString(),
        params.embedding_model,
        params.embedding_batch_size,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    slides            = HISTOPATHOLOGY.out.slides            // channel: [ val(meta), path(slide) ]
    normalized_slides = HISTOPATHOLOGY.out.normalized_slides // channel: [ val(meta), path(slide) ]
    tiles             = HISTOPATHOLOGY.out.tiles             // channel: [ val(meta), path(tiles) ]
    tile_qc           = HISTOPATHOLOGY.out.tile_qc           // channel: [ val(meta), path(scored_tiles) ]
    slide_qc          = HISTOPATHOLOGY.out.slide_qc          // channel: [ val(meta), path(artifacts) ]
    cells             = HISTOPATHOLOGY.out.cells             // channel: [ val(meta), path(instances) ]
    cell_types        = HISTOPATHOLOGY.out.cell_types        // channel: [ val(meta), path(typed_instances) ]
    embeddings        = HISTOPATHOLOGY.out.embeddings        // channel: [ val(meta), path(embedding_store) ]
    multiqc_report    = HISTOPATHOLOGY.out.multiqc_report    // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        params.stain_reference?.toString(),
        params.segmentation_model,
        params.segmentation_mpp,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_HISTOPATHOLOGY (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_HISTOPATHOLOGY.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
