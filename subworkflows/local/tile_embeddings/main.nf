include { LAZYSLIDE_FEATUREEXTRACTION                        } from '../../../modules/local/lazyslide/featureextraction/main'
include { LAZYSLIDE_TILETISSUES as EMBEDDING_TILETISSUES     } from '../../../modules/local/lazyslide/tiletissues/main'
include { LAZYSLIDE_TILETISSUES as EMBEDDING_TILETISSUES_10X } from '../../../modules/local/lazyslide/tiletissues/main'

workflow TILE_EMBEDDINGS {
    take:
    ch_slides            // channel: [ val(meta), path(slide) ]
    ch_slides_with_masks // channel: [ val(meta), path(slide), path(tissue_mask) ]
    ch_tiles             // channel: [ val(meta), path(tiles) ]
    tile_px              // value: regular tile width in pixels
    tile_mpp             // value: regular tile resolution
    tile_overlap         // value: regular tile overlap
    background_fraction  // value: maximum background fraction
    models               // value: list of embedding model names
    batch_size           // value: inference batch size
    grid_overrides       // value: model name to [tile size, MPP] for paired slide encoders

    main:
    def model_tile_sizes = [
        conceptclip: 384,
        'conch_v1.5': 448,
        medsiglip: 448,
        musk: 384,
        titan: 448,
    ]
    def model_grids = models.collectEntries { model ->
        [model, grid_overrides.get(model, [model_tile_sizes.get(model, 256), 0.5])]
    }
    def embedding_grids = model_grids.values().toList().unique()
    def regular_grid = tile_mpp && tile_overlap.toString().toBigDecimal() == 0
        ? embedding_grids.find { grid -> grid[0] == tile_px.toString().toInteger() && grid[1] == tile_mpp.toString().toBigDecimal() }
        : null
    def grids_to_build = embedding_grids - [regular_grid]
    def ch_embedding_tiles = regular_grid
        ? ch_tiles.map { meta, tiles -> [meta, tiles, regular_grid] }
        : channel.empty()
    def sizes_20x = grids_to_build.findAll { grid -> grid[1] == 0.5 }.collect { grid -> grid[0] }
    if (sizes_20x) {
        EMBEDDING_TILETISSUES(ch_slides_with_masks, sizes_20x, 0.5, 0, background_fraction)
        ch_embedding_tiles = ch_embedding_tiles.mix(EMBEDDING_TILETISSUES.out.tiles.map { meta, tiles, pixels -> [meta, tiles, [pixels, 0.5]] })
    }
    def sizes_10x = grids_to_build.findAll { grid -> grid[1] == 1.0 }.collect { grid -> grid[0] }
    if (sizes_10x) {
        EMBEDDING_TILETISSUES_10X(ch_slides_with_masks, sizes_10x, 1.0, 0, background_fraction)
        ch_embedding_tiles = ch_embedding_tiles.mix(EMBEDDING_TILETISSUES_10X.out.tiles.map { meta, tiles, pixels -> [meta, tiles, [pixels, 1.0]] })
    }

    def ch_models = channel.fromList(models).map { model -> tuple(model_grids[model], model) }
    def ch_inputs = ch_slides
        .combine(ch_embedding_tiles, by: 0)
        .map { meta, slide, tiles, grid -> tuple(grid, meta, slide, tiles) }
        .combine(ch_models, by: 0)
        .map { _grid, meta, slide, tiles, model -> tuple(meta, slide, tiles, model) }
    LAZYSLIDE_FEATUREEXTRACTION(ch_inputs, batch_size)

    emit:
    embeddings = LAZYSLIDE_FEATUREEXTRACTION.out.embeddings // channel: [ val(meta), path(embedding_store), val(model) ]
}
