include { LAZYSLIDE_CELLS as LAZYSLIDE_CELLS_INSTANSEG      } from '../../../modules/local/lazyslide/cells/main'
include { LAZYSLIDE_CELLS as LAZYSLIDE_CELLS_CELLPOSE       } from '../../../modules/local/lazyslide/cells/main'
include { LAZYSLIDE_CELLS as LAZYSLIDE_CELLS_NULITE         } from '../../../modules/local/lazyslide/cells/main'
include { LAZYSLIDE_CELLS as LAZYSLIDE_CELLS_HISTOPLUS      } from '../../../modules/local/lazyslide/cells/main'
include { LAZYSLIDE_TILETISSUES as SEGMENTATION_TILETISSUES } from '../../../modules/local/lazyslide/tiletissues/main'
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
    def model_pixels = [instanseg: 512, cellpose: 256, nulite: 256, histoplus: 224]
    def selected_sizes = models.collect { m -> model_pixels[m] }.unique()
    def segmentation_overlap = 0.125
    def regular_pixels = tile_px.toString().toInteger()
    def reuse_grid = selected_sizes.contains(regular_pixels) && tile_mpp &&
        tile_mpp.toString().toBigDecimal() == segmentation_mpp.toString().toBigDecimal() &&
        tile_overlap.toString().toBigDecimal() == segmentation_overlap

    def ch_segmentation_tiles = reuse_grid
        ? ch_tiles.map { meta, tiles -> [meta, tiles, regular_pixels] }
        : channel.empty()
    def sizes_to_build = reuse_grid ? selected_sizes - regular_pixels : selected_sizes
    if (sizes_to_build) {
        SEGMENTATION_TILETISSUES(ch_slides_with_masks, sizes_to_build, segmentation_mpp, segmentation_overlap, background_fraction)
        ch_segmentation_tiles = ch_segmentation_tiles.mix(SEGMENTATION_TILETISSUES.out.tiles)
    }

    def grids = ch_slides.combine(ch_segmentation_tiles, by: 0)
        .branch { meta, slide, tiles, pixels ->
            tiles_512: pixels == 512
                return [meta, slide, tiles]
            tiles_256: pixels == 256
                return [meta, slide, tiles]
            tiles_224: pixels == 224
                return [meta, slide, tiles]
        }

    def ch_cells = channel.empty()
    if ('instanseg' in models) {
        LAZYSLIDE_CELLS_INSTANSEG(grids.tiles_512, 'instanseg', file(instanseg_checkpoint, checkIfExists: true))
        ch_cells = ch_cells.mix(LAZYSLIDE_CELLS_INSTANSEG.out.cells)
    }
    if ('cellpose' in models) {
        LAZYSLIDE_CELLS_CELLPOSE(grids.tiles_256, 'cellpose', file(cellpose_checkpoint, checkIfExists: true))
        ch_cells = ch_cells.mix(LAZYSLIDE_CELLS_CELLPOSE.out.cells)
    }
    if ('nulite' in models) {
        LAZYSLIDE_CELLS_NULITE(grids.tiles_256, 'nulite', file(nulite_checkpoint, checkIfExists: true))
        ch_cells = ch_cells.mix(LAZYSLIDE_CELLS_NULITE.out.cells)
    }
    if ('histoplus' in models) {
        def checkpoint
        if (histoplus_checkpoint) {
            checkpoint = file(histoplus_checkpoint, checkIfExists: true)
        }
        else {
            HUGGINGFACE_DOWNLOAD_HISTOPLUS(channel.value([
                [id: 'histoplus'],
                'model',
                'Owkin-Bioptimus/histoplus',
                'cde2eee81af9e39b03802fc33d4f284733b5ee5e',
                'histoplus_cellvit_segmentor_20x.pt',
            ]))
            checkpoint = HUGGINGFACE_DOWNLOAD_HISTOPLUS.out.resource.map { _meta, resource -> resource }
        }
        LAZYSLIDE_CELLS_HISTOPLUS(grids.tiles_224, 'histoplus', checkpoint)
        ch_cells = ch_cells.mix(LAZYSLIDE_CELLS_HISTOPLUS.out.cells)
    }

    emit:
    cells = ch_cells // channel: [ val(meta), path(instances) ]
}
