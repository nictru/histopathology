/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { WSIDATA_INGEST         } from '../modules/local/wsidata/ingest/main'
include { LAZYSLIDE_FINDTISSUES   } from '../modules/local/lazyslide/findtissues/main'
include { LAZYSLIDE_TILETISSUES   } from '../modules/local/lazyslide/tiletissues/main'
include { QC                     } from '../subworkflows/local/qc/main'
include { SEGMENTATION           } from '../subworkflows/local/segmentation/main'
include { STAIN_NORMALIZATION    } from '../subworkflows/local/stain_normalization/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_histopathology_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HISTOPATHOLOGY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    stain_normalization
    stain_reference
    tile_px
    tile_mpp
    tile_overlap
    tile_background_fraction
    focuslitenn_checkpoint
    pathprofilerqc_checkpoint
    grandqc_checkpoint
    grandqc_variant
    segmentation_model
    segmentation_mpp
    instanseg_checkpoint
    cellpose_checkpoint
    nulite_checkpoint
    histoplus_checkpoint
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_ingest = ch_samplesheet.map { meta, slide ->
        [meta, slide, meta.scene, meta.mpp]
    }
    WSIDATA_INGEST(ch_ingest)
    LAZYSLIDE_FINDTISSUES(WSIDATA_INGEST.out.slide)

    def ch_slides_with_masks = WSIDATA_INGEST.out.slide.join(LAZYSLIDE_FINDTISSUES.out.mask)

    LAZYSLIDE_TILETISSUES(ch_slides_with_masks, tile_px, tile_mpp, tile_overlap, tile_background_fraction)
    def ch_tiles = LAZYSLIDE_TILETISSUES.out.tiles.map { meta, tiles, _pixels -> [meta, tiles] }
    QC(
        ch_slides_with_masks,
        ch_tiles,
        focuslitenn_checkpoint,
        pathprofilerqc_checkpoint,
        grandqc_checkpoint,
        grandqc_variant,
    )

    def ch_normalized_slides = channel.empty()
    if (stain_reference) {
        STAIN_NORMALIZATION(ch_slides_with_masks, stain_normalization, stain_reference)
        ch_normalized_slides = STAIN_NORMALIZATION.out.slides
    }

    def segmentation_models = (!segmentation_model || segmentation_model == 'none')
        ? []
        : segmentation_model.split(',', -1).collect { it.trim() }.unique()
    def ch_cells = channel.empty()
    if (segmentation_models) {
        def unsupported_models = segmentation_models - ['instanseg', 'cellpose', 'nulite', 'histoplus']
        if (unsupported_models) {
            error("Unsupported segmentation models: ${unsupported_models.join(', ')}.")
        }
        SEGMENTATION(
            stain_reference ? ch_normalized_slides : WSIDATA_INGEST.out.slide,
            ch_slides_with_masks,
            ch_tiles,
            tile_px,
            tile_mpp,
            tile_overlap,
            tile_background_fraction,
            segmentation_models,
            segmentation_mpp,
            instanseg_checkpoint,
            cellpose_checkpoint,
            nulite_checkpoint,
            histoplus_checkpoint,
        )
        ch_cells = SEGMENTATION.out.cells
    }

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'histopathology_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'histopathology'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    slides            = WSIDATA_INGEST.out.slide                                  // channel: [ val(meta), path(slide) ]
    normalized_slides = ch_normalized_slides                                     // channel: [ val(meta), path(slide) ]
    tiles             = ch_tiles                                                  // channel: [ val(meta), path(tiles) ]
    tile_qc           = QC.out.tile_qc                                            // channel: [ val(meta), path(scored_tiles) ]
    slide_qc          = QC.out.slide_qc                                           // channel: [ val(meta), path(artifacts) ]
    cells             = ch_cells                                                  // channel: [ val(meta), path(instances) ]
    multiqc_report    = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions          = ch_versions                                              // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
