include { LAZYSLIDE_FEATUREEXTRACTION                             } from '../../../modules/local/lazyslide/featureextraction/main'
include { LAZYSLIDE_TILETISSUES as EMBEDDING_TILETISSUES          } from '../../../modules/local/lazyslide/tiletissues/main'

workflow TILE_EMBEDDINGS {
    take:
    ch_slides                 // channel: [ val(meta), path(slide) ]
    ch_slides_with_masks      // channel: [ val(meta), path(slide), path(tissue_mask) ]
    ch_tiles                  // channel: [ val(meta), path(tiles) ]
    tile_px                   // value: regular tile width in pixels
    tile_mpp                  // value: regular tile resolution
    tile_overlap              // value: regular tile overlap
    background_fraction       // value: maximum background fraction
    models                    // value: list of embedding model names
    batch_size                // value: inference batch size

    main:
    def model_tile_sizes = [
        conceptclip: 384,
        'conch_v1.5': 448,
        medsiglip: 448,
        musk: 384,
        titan: 448,
    ]
    def embedding_sizes = models.collect { model -> model_tile_sizes.get(model, 256) }.unique()
    def embedding_mpp = 0.5
    def regular_pixels = tile_px.toString().toInteger()
    def reuse_grid = embedding_sizes.contains(regular_pixels) && tile_mpp && tile_mpp.toString().toBigDecimal() == embedding_mpp && tile_overlap.toString().toBigDecimal() == 0
    def ch_embedding_tiles = reuse_grid
        ? ch_tiles.map { meta, tiles -> [meta, tiles, regular_pixels] }
        : channel.empty()
    def sizes_to_build = reuse_grid ? embedding_sizes - regular_pixels : embedding_sizes
    if (sizes_to_build) {
        EMBEDDING_TILETISSUES(ch_slides_with_masks, sizes_to_build, embedding_mpp, 0, background_fraction)
        ch_embedding_tiles = ch_embedding_tiles.mix(EMBEDDING_TILETISSUES.out.tiles)
    }

    def ch_models = channel.fromList(models).map { model -> tuple(model_tile_sizes.get(model, 256), model) }
    def ch_inputs = ch_slides.combine(ch_embedding_tiles, by: 0)
        .map { meta, slide, tiles, pixels -> tuple(pixels, meta, slide, tiles) }
        .combine(ch_models, by: 0)
        .map { _pixels, meta, slide, tiles, model -> tuple(meta, slide, tiles, model) }
    LAZYSLIDE_FEATUREEXTRACTION(ch_inputs, batch_size)

    emit:
    embeddings = LAZYSLIDE_FEATUREEXTRACTION.out.embeddings // channel: [ val(meta), path(embedding_store) ]
}
