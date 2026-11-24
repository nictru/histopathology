include { LAZYSLIDE_CELLS                                        } from '../../../modules/local/lazyslide/cells/main'
include { LAZYSLIDE_TILETISSUES as SEGMENTATION_TILETISSUES      } from '../../../modules/local/lazyslide/tiletissues/main'
include { HUGGINGFACE_DOWNLOAD as HUGGINGFACE_DOWNLOAD_HISTOPLUS } from '../../../modules/local/huggingface/download/main'

workflow SEGMENTATION {
    take:
    ch_slides            // channel: [ val(meta), path(slide) ]
    ch_slides_with_masks // channel: [ val(meta), path(slide), path(tissue_mask) ]
    ch_tiles             // channel: [ val(meta), path(tiles) ]
    tile_px              // value: regular tile width in pixels
    tile_mpp             // value: regular tile resolution
    tile_overlap         // value: regular tile overlap
    background_fraction  // value: maximum background fraction
    models               // value: list of segmentation model names
    segmentation_mpp     // value: segmentation tile resolution
    instanseg_checkpoint // value: InstanSeg checkpoint path or URL
    cellpose_checkpoint  // value: Cellpose checkpoint path or URL
    nulite_checkpoint    // value: NuLite checkpoint path or URL
    histoplus_checkpoint // value: HistoPLUS checkpoint path or URL

    main:
    def sources = [
        instanseg: [pixels: 512, checkpoint: instanseg_checkpoint],
        cellpose: [pixels: 256, checkpoint: cellpose_checkpoint],
        nulite: [pixels: 256, checkpoint: nulite_checkpoint],
        histoplus: [pixels: 224, checkpoint: histoplus_checkpoint],
    ]
    def selected_models = sources.keySet().findAll { model -> model in models }
    def selected_sizes = selected_models.collect { model -> sources[model].pixels }.unique()
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

    def download_histoplus = 'histoplus' in selected_models && !histoplus_checkpoint
    def provided_models = download_histoplus ? selected_models - 'histoplus' : selected_models
    def ch_checkpoints = channel.fromList(provided_models)
        .map { model ->
            tuple(sources[model].pixels, model, file(sources[model].checkpoint, checkIfExists: true))
        }
    if (download_histoplus) {
        HUGGINGFACE_DOWNLOAD_HISTOPLUS(
            channel.value(
                [
                    [id: 'histoplus'],
                    'model',
                    'Owkin-Bioptimus/histoplus',
                    'cde2eee81af9e39b03802fc33d4f284733b5ee5e',
                    'histoplus_cellvit_segmentor_20x.pt',
                ]
            )
        )
        def ch_downloaded = HUGGINGFACE_DOWNLOAD_HISTOPLUS.out.resource.map { meta, checkpoint -> tuple(sources[meta.id].pixels, meta.id, checkpoint) }
        ch_checkpoints = ch_checkpoints.mix(ch_downloaded)
    }

    def ch_grids = ch_slides
        .combine(ch_segmentation_tiles, by: 0)
        .map { meta, slide, tiles, pixels -> tuple(pixels, meta, slide, tiles) }
    def ch_inputs = ch_grids
        .combine(ch_checkpoints, by: 0)
        .map { _pixels, meta, slide, tiles, model, checkpoint -> tuple(meta, slide, tiles, model, checkpoint) }
    LAZYSLIDE_CELLS(ch_inputs)

    emit:
    cells = LAZYSLIDE_CELLS.out.cells // channel: [ val(meta), path(instances) ]
}
