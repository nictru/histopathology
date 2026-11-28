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
include { CYTOFORMER_CLASSIFY    } from '../modules/local/cytoformer/classify/main'
include { SPATIALDATA_CELLTYPECOMPOSITION } from '../modules/local/spatialdata/celltypecomposition/main'
include { SQUIDPY_NEIGHBORHOODENRICHMENT } from '../modules/local/squidpy/neighborhoodenrichment/main'
include { SQUIDPY_NICHES         } from '../modules/local/squidpy/niches/main'
include { SPATIALDATA_CELLFEATURES } from '../modules/local/spatialdata/cellfeatures/main'
include { STAIN_NORMALIZATION    } from '../subworkflows/local/stain_normalization/main'
include { TILE_EMBEDDINGS         } from '../subworkflows/local/tile_embeddings/main'
include { LAZYSLIDE_FEATUREAGGREGATION } from '../modules/local/lazyslide/featureaggregation/main'
include { LAZYSLIDE_SPATIALDOMAIN } from '../modules/local/lazyslide/spatialdomain/main'
include { LAZYSLIDE_TISSUEREGIONS } from '../modules/local/lazyslide/tissueregions/main'
include { SPATIALDATA_QUPATHEXPORT } from '../modules/local/spatialdata/qupathexport/main'
include { QUPATH_PROJECT         } from '../modules/local/qupath/project/main'
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
    tissue_segmentation_method
    tile_px
    tile_mpp
    tile_overlap
    tile_background_fraction
    grandqc_variant
    segmentation_model
    segmentation_mpp
    cytoformer_organ
    cell_neighborhood_radius_um
    cell_neighborhood_permutations
    cell_niche_radius_um
    cell_niche_resolution
    cell_features
    cell_region_domain_model
    embedding_model
    embedding_batch_size
    slide_embedding_method
    tissue_domain_resolution
    export_qupath
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def segmentation_models = (!segmentation_model || segmentation_model == 'none')
        ? []
        : segmentation_model.split(',', -1).collect { model -> model.trim() }.unique()
    if (cytoformer_organ && !segmentation_models) {
        error('CytoFormer requires at least one --segmentation_model.')
    }
    if (cell_neighborhood_radius_um != null && !segmentation_models) {
        error('Cell-type neighborhood enrichment requires at least one --segmentation_model.')
    }
    def native_cell_typing_models = ['classpose', 'nulite', 'histoplus'].intersect(segmentation_models)
    if (cell_niche_radius_um != null && !cytoformer_organ && !native_cell_typing_models) {
        error('Spatial cell niches require Classpose, NuLite, HistoPLUS, or CytoFormer cell typing.')
    }
    if (cell_features && !segmentation_models) {
        error('Cell features require at least one --segmentation_model.')
    }

    def ch_ingest = ch_samplesheet.map { meta, slide ->
        [meta, slide, meta.scene, meta.mpp]
    }
    WSIDATA_INGEST(ch_ingest)
    def tissue_methods = tissue_segmentation_method.split(',', -1).collect { method -> method.trim() }.unique()
    def ch_tissue_inputs = WSIDATA_INGEST.out.slide.combine(channel.fromList(tissue_methods))
    LAZYSLIDE_FINDTISSUES(ch_tissue_inputs)

    def ch_primary_masks = LAZYSLIDE_FINDTISSUES.out.mask
        .filter { _meta, _mask, method -> method == tissue_methods[0] }
        .map { meta, mask, _method -> tuple(meta, mask) }
    def ch_slides_with_masks = WSIDATA_INGEST.out.slide.join(ch_primary_masks)

    LAZYSLIDE_TILETISSUES(ch_slides_with_masks, tile_px, tile_mpp, tile_overlap, tile_background_fraction)
    def ch_tiles = LAZYSLIDE_TILETISSUES.out.tiles.map { meta, tiles, _pixels -> [meta, tiles] }
    QC(
        ch_slides_with_masks,
        ch_tiles,
        grandqc_variant,
    )

    def ch_multiqc_files = LAZYSLIDE_FINDTISSUES.out.multiqc_files.mix(QC.out.multiqc_files)

    def ch_normalized_slides = channel.empty()
    if (stain_reference) {
        STAIN_NORMALIZATION(ch_slides_with_masks, stain_normalization, stain_reference)
        ch_normalized_slides = STAIN_NORMALIZATION.out.slides
        ch_multiqc_files = ch_multiqc_files.mix(STAIN_NORMALIZATION.out.multiqc_files)
    }
    def ch_analysis_slides = stain_reference ? ch_normalized_slides : WSIDATA_INGEST.out.slide

    def ch_cells = channel.empty()
    def ch_segmentation_agreement = channel.empty()
    if (segmentation_models) {
        def unsupported_models = segmentation_models - ['instanseg', 'cellpose', 'classpose', 'nulite', 'histoplus']
        if (unsupported_models) {
            error("Unsupported segmentation models: ${unsupported_models.join(', ')}.")
        }
        SEGMENTATION(
            ch_analysis_slides,
            ch_slides_with_masks,
            ch_tiles,
            tile_px,
            tile_mpp,
            tile_overlap,
            tile_background_fraction,
            segmentation_models,
            segmentation_mpp,
        )
        ch_cells = SEGMENTATION.out.cells
        ch_segmentation_agreement = SEGMENTATION.out.agreement
        ch_multiqc_files = ch_multiqc_files.mix(SEGMENTATION.out.multiqc_files)
    }

    def ch_cell_types = channel.empty()
    if (cytoformer_organ) {
        def ch_typing_inputs = ch_cells.combine(WSIDATA_INGEST.out.slide, by: 0)
        CYTOFORMER_CLASSIFY(ch_typing_inputs, cytoformer_organ)
        ch_cell_types = CYTOFORMER_CLASSIFY.out.cell_types
        ch_multiqc_files = ch_multiqc_files.mix(CYTOFORMER_CLASSIFY.out.multiqc_files)
    }

    def ch_classification_sources = ch_cells.mix(ch_cell_types)
    SPATIALDATA_CELLTYPECOMPOSITION(ch_classification_sources.combine(ch_slides_with_masks, by: 0))
    ch_multiqc_files = ch_multiqc_files.mix(SPATIALDATA_CELLTYPECOMPOSITION.out.multiqc_files)

    def ch_cell_neighborhoods = channel.empty()
    if (cell_neighborhood_radius_um != null) {
        SQUIDPY_NEIGHBORHOODENRICHMENT(
            ch_classification_sources.combine(WSIDATA_INGEST.out.slide, by: 0),
            cell_neighborhood_radius_um,
            cell_neighborhood_permutations,
        )
        ch_cell_neighborhoods = SQUIDPY_NEIGHBORHOODENRICHMENT.out.neighborhoods
        ch_multiqc_files = ch_multiqc_files.mix(SQUIDPY_NEIGHBORHOODENRICHMENT.out.multiqc_files)
    }

    def ch_cell_niches = channel.empty()
    if (cell_niche_radius_um != null) {
        def ch_native_cell_types = ch_cells.filter { _meta, store ->
            native_cell_typing_models.any { model -> store.name.endsWith(".${model}.cells.zarr") }
        }
        SQUIDPY_NICHES(
            ch_native_cell_types.mix(ch_cell_types).combine(WSIDATA_INGEST.out.slide, by: 0),
            cell_niche_radius_um,
            cell_niche_resolution,
        )
        ch_cell_niches = SQUIDPY_NICHES.out.niches
        ch_multiqc_files = ch_multiqc_files.mix(SQUIDPY_NICHES.out.multiqc_files)
    }

    def embedding_models = (!embedding_model || embedding_model == 'none')
        ? []
        : embedding_model.split(',', -1).collect { model -> model.trim() }.unique()
    def slide_methods = (!slide_embedding_method || slide_embedding_method == 'none')
        ? []
        : slide_embedding_method.split(',', -1).collect { method -> method.trim() }.unique()
    def pooling_methods = ['mean', 'median', 'max']
    def slide_pairs = [
        chief:     [tile_model: 'chief',               grid: [256, 1.0]],
        madeleine: [tile_model: 'conch-madeleine',      grid: [256, 1.0]],
        moozy:     [tile_model: 'lunit-dino-s8-moozy',  grid: [224, 0.5]],
        prism:     [tile_model: 'virchow',              grid: [224, 0.5]],
        titan:     [tile_model: 'titan',                grid: [512, 0.5]],
    ]
    def unsupported_methods = slide_methods - pooling_methods - slide_pairs.keySet()
    if (unsupported_methods) {
        error("Unsupported slide embedding methods: ${unsupported_methods.join(', ')}.")
    }
    def learned_methods = slide_methods - pooling_methods
    def required_grids = learned_methods.collectEntries { method ->
        def pair = slide_pairs[method]
        [(pair.tile_model): pair.grid]
    }
    embedding_models = (embedding_models + required_grids.keySet()).unique()
    if (slide_methods && !embedding_models) {
        error('Slide pooling requires --embedding_model or a learned --slide_embedding_method.')
    }
    if (tissue_domain_resolution != null && !embedding_models) {
        error('Tissue domain maps require --embedding_model or a learned --slide_embedding_method.')
    }
    if (cell_region_domain_model && (tissue_domain_resolution == null || !(cell_region_domain_model in embedding_models))) {
        error('--cell_region_domain_model requires --tissue_domain_resolution and a selected tile encoder.')
    }
    def ch_embeddings = channel.empty()
    def ch_slide_embeddings = channel.empty()
    def ch_tissue_domains = channel.empty()
    if (embedding_models) {
        def unsupported_models = embedding_models.findAll { model -> !(model ==~ /[a-z0-9][a-z0-9_.-]*/) || model == 'none' }
        if (unsupported_models) {
            error("Unsupported embedding models: ${unsupported_models.join(', ')}.")
        }
        TILE_EMBEDDINGS(
            ch_analysis_slides,
            ch_slides_with_masks,
            ch_tiles,
            tile_px,
            tile_mpp,
            tile_overlap,
            tile_background_fraction,
            embedding_models,
            embedding_batch_size,
            required_grids,
        )
        ch_embeddings = TILE_EMBEDDINGS.out.embeddings.map { meta, store, _model -> [meta, store] }
        if (tissue_domain_resolution != null) {
            LAZYSLIDE_SPATIALDOMAIN(
                TILE_EMBEDDINGS.out.embeddings.combine(WSIDATA_INGEST.out.slide, by: 0),
                tissue_domain_resolution,
            )
            ch_tissue_domains = LAZYSLIDE_SPATIALDOMAIN.out.domains
            ch_multiqc_files = ch_multiqc_files.mix(LAZYSLIDE_SPATIALDOMAIN.out.multiqc_files)
        }
        if (slide_methods) {
            def ch_slide_inputs = TILE_EMBEDDINGS.out.embeddings
                .combine(ch_analysis_slides, by: 0)
                .flatMap { meta, store, model, slide ->
                    slide_methods.findAll { method -> method in pooling_methods || slide_pairs[method].tile_model == model }
                        .collect { method -> [meta, store, slide, model, method] }
                }
            LAZYSLIDE_FEATUREAGGREGATION(ch_slide_inputs)
            ch_slide_embeddings = LAZYSLIDE_FEATUREAGGREGATION.out.embeddings
        }
    }

    def ch_tissue_regions = channel.empty()
    if (cell_region_domain_model) {
        def ch_region_domains = ch_tissue_domains.filter { _meta, store ->
            store.name.endsWith(".${cell_region_domain_model}.tissue_domains.zarr")
        }
        LAZYSLIDE_TISSUEREGIONS(ch_region_domains.join(ch_slides_with_masks), cell_region_domain_model)
        ch_tissue_regions = LAZYSLIDE_TISSUEREGIONS.out.regions
        ch_multiqc_files = ch_multiqc_files.mix(LAZYSLIDE_TISSUEREGIONS.out.multiqc_files)
    }

    def ch_cell_features = channel.empty()
    def ch_region_composition = channel.empty()
    if (cell_features || (cell_region_domain_model && segmentation_models)) {
        def ch_features_with_slide = ch_classification_sources.combine(WSIDATA_INGEST.out.slide, by: 0)
        def ch_feature_inputs = ch_features_with_slide.map { meta, instances, slide -> [meta, instances, slide, [], []] }
        if (cell_region_domain_model) {
            def ch_features_by_source = ch_features_with_slide.combine(ch_tissue_regions, by: 0)
                .map { meta, instances, slide, regions ->
                    [[meta.id, instances.name.replaceFirst(/(\.cells)?\.zarr$/, '')], meta, instances, slide, regions]
                }
            def ch_niches_by_source = ch_cell_niches.map { meta, store ->
                [[meta.id, store.name.replaceFirst(/\.cell_niches\.zarr$/, '')], store]
            }
            ch_feature_inputs = ch_features_by_source.join(ch_niches_by_source, remainder: true, failOnDuplicate: true)
                .map { _key, meta, instances, slide, regions, niches -> [meta, instances, slide, regions, niches ?: []] }
        }
        SPATIALDATA_CELLFEATURES(ch_feature_inputs)
        ch_cell_features = SPATIALDATA_CELLFEATURES.out.features
        ch_region_composition = SPATIALDATA_CELLFEATURES.out.composition
        ch_multiqc_files = ch_multiqc_files.mix(SPATIALDATA_CELLFEATURES.out.multiqc_files)
    }

    def ch_qupath_project = channel.empty()
    if (export_qupath) {
        def ch_export_stores = WSIDATA_INGEST.out.slide
            .map { meta, store -> [meta, 'original', store] }
            .mix(
                ch_normalized_slides.map { meta, store -> [meta, 'normalized', store] },
                LAZYSLIDE_FINDTISSUES.out.mask.map { meta, store, _method -> [meta, 'tissue', store] },
                ch_tiles.map { meta, store -> [meta, 'tiles', store] },
                QC.out.tile_qc.map { meta, store -> [meta, 'tile_qc', store] },
                QC.out.slide_qc.map { meta, store -> [meta, 'slide_qc', store] },
                ch_cells.map { meta, store -> [meta, 'cells', store] },
                ch_cell_types.map { meta, store -> [meta, 'cell_types', store] },
                ch_cell_niches.map { meta, store -> [meta, 'cell_niches', store] },
                ch_cell_features.map { meta, store -> [meta, 'cell_features', store] },
                ch_tissue_domains.map { meta, store -> [meta, 'domains', store] },
                ch_tissue_regions.map { meta, store -> [meta, 'regions', store] },
            )
            .groupTuple()
            .map { meta, kinds, stores ->
                def order = kinds.indices.toList().sort { i -> "${kinds[i]}:${stores[i].name}" }
                [meta, order.collect { i -> kinds[i] }, order.collect { i -> stores[i] }]
            }
        SPATIALDATA_QUPATHEXPORT(ch_export_stores)
        def ch_qupath_slides = SPATIALDATA_QUPATHEXPORT.out.slide
            .map { _meta, slide -> slide }
            .toSortedList { a, b -> a.name <=> b.name }
            .filter { slides -> slides.size() > 0 }
        def archive = workflow.stubRun
            ? file("${projectDir}/nextflow.config")
            : file('https://github.com/qupath/qupath/releases/download/v0.7.0/QuPath-v0.7.0-Linux.tar.xz')
        QUPATH_PROJECT(ch_qupath_slides, archive)
        ch_qupath_project = QUPATH_PROJECT.out.project
    }

    def ch_versions = channel.empty()

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
    cell_types        = ch_cell_types                                             // channel: [ val(meta), path(typed_instances) ]
    cell_type_composition = SPATIALDATA_CELLTYPECOMPOSITION.out.composition         // channel: [ val(meta), path(composition_csv) ]
    cell_neighborhoods = ch_cell_neighborhoods                                    // channel: [ val(meta), path(neighborhood_csv) ]
    cell_niches       = ch_cell_niches                                            // channel: [ val(meta), path(niche_store) ]
    cell_features     = ch_cell_features                                          // channel: [ val(meta), path(feature_store) ]
    region_composition = ch_region_composition                                   // channel: [ val(meta), path(composition_csv) ]
    segmentation_agreement = ch_segmentation_agreement                            // channel: [ val(meta), path(agreement_csv) ]
    embeddings        = ch_embeddings                                             // channel: [ val(meta), path(embedding_store) ]
    slide_embeddings  = ch_slide_embeddings                                       // channel: [ val(meta), path(slide_embedding) ]
    tissue_domains    = ch_tissue_domains                                         // channel: [ val(meta), path(domain_store) ]
    tissue_regions    = ch_tissue_regions                                         // channel: [ val(meta), path(region_store) ]
    qupath_project    = ch_qupath_project                                         // channel: path(project_directory)
    multiqc_report    = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions          = ch_versions                                              // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
