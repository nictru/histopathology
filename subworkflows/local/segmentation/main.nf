include { LAZYSLIDE_CELLS                                   } from '../../../modules/local/lazyslide/cells/main'
include { LAZYSLIDE_TILETISSUES as SEGMENTATION_TILETISSUES } from '../../../modules/local/lazyslide/tiletissues/main'
include { SPATIALDATA_SEGMENTATIONAGREEMENT                 } from '../../../modules/local/spatialdata/segmentationagreement/main'

workflow SEGMENTATION {
    take:
    ch_slides            // channel: [ val(meta), path(slide) ]
    ch_slides_with_masks // channel: [ val(meta), path(original), path(tissue_mask) ]
    ch_tiles             // channel: [ val(meta), path(tiles) ]
    tile_px              // value: regular tile width in pixels
    tile_mpp             // value: regular tile resolution
    tile_overlap         // value: regular tile overlap
    background_fraction  // value: maximum background fraction
    models               // value: list of segmentation model names
    segmentation_mpp     // value: segmentation tile resolution

    main:
    def tile_sizes = [
        instanseg: 512,
        cellpose: 256,
        classpose: 256,
        nulite: 256,
        histoplus: 224,
    ]
    def selected_models = tile_sizes.keySet().findAll { model -> model in models }
    def selected_sizes = selected_models.collect { model -> tile_sizes[model] }.unique()
    def segmentation_overlap = 0.125
    def regular_pixels = tile_px.toString().toInteger()
    def reuse_grid = selected_sizes.contains(regular_pixels) && tile_mpp && tile_mpp.toString().toBigDecimal() == segmentation_mpp.toString().toBigDecimal() && tile_overlap.toString().toBigDecimal() == segmentation_overlap

    def ch_segmentation_tiles = reuse_grid
        ? ch_tiles.map { meta, tiles -> [meta, tiles, regular_pixels] }
        : channel.empty()
    def sizes_to_build = reuse_grid ? selected_sizes - regular_pixels : selected_sizes
    if (sizes_to_build) {
        SEGMENTATION_TILETISSUES(ch_slides_with_masks, sizes_to_build, segmentation_mpp, segmentation_overlap, background_fraction)
        ch_segmentation_tiles = ch_segmentation_tiles.mix(SEGMENTATION_TILETISSUES.out.tiles)
    }

    def ch_models = channel.fromList(selected_models).map { model -> tuple(tile_sizes[model], model) }

    def ch_originals = ch_slides_with_masks.map { meta, original, _mask -> tuple(meta, original) }
    def ch_grids = ch_slides
        .join(ch_originals)
        .combine(ch_segmentation_tiles, by: 0)
        .map { meta, slide, original, tiles, pixels -> tuple(pixels, meta, slide, original, tiles) }
    def ch_inputs = ch_grids
        .combine(ch_models, by: 0)
        .map { _pixels, meta, slide, original, tiles, model -> tuple(meta, slide, original, tiles, model) }
    LAZYSLIDE_CELLS(ch_inputs)

    def ch_agreement = channel.empty()
    def ch_multiqc_files = LAZYSLIDE_CELLS.out.multiqc_files.mix(LAZYSLIDE_CELLS.out.cell_types_multiqc)
    if (selected_models.size() > 1) {
        def ch_segmentations = LAZYSLIDE_CELLS.out.cells
            .groupTuple(size: selected_models.size())
            .map { meta, stores -> tuple(meta, stores.sort { store -> store.name }) }
        SPATIALDATA_SEGMENTATIONAGREEMENT(ch_segmentations)
        ch_agreement = SPATIALDATA_SEGMENTATIONAGREEMENT.out.agreement
        ch_multiqc_files = ch_multiqc_files.mix(SPATIALDATA_SEGMENTATIONAGREEMENT.out.multiqc_files)
    }

    emit:
    cells         = LAZYSLIDE_CELLS.out.cells // channel: [ val(meta), path(instances) ]
    agreement     = ch_agreement              // channel: [ val(meta), path(agreement_csv) ]
    multiqc_files = ch_multiqc_files          // channel: [ path(multiqc_file) ]
}
