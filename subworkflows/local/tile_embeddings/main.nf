include { LAZYSLIDE_FEATUREEXTRACTION                             } from '../../../modules/local/lazyslide/featureextraction/main'
include { LAZYSLIDE_TILETISSUES as EMBEDDING_TILETISSUES          } from '../../../modules/local/lazyslide/tiletissues/main'
include { HUGGINGFACE_DOWNLOAD as HUGGINGFACE_DOWNLOAD_EMBEDDINGS } from '../../../modules/local/huggingface/download/main'

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
    resnet50_checkpoint       // value: ResNet50 checkpoint path or URL
    uni_checkpoint            // value: optional UNI checkpoint path or URL
    uni2_checkpoint           // value: optional UNI2 checkpoint path or URL
    gigapath_flash_checkpoint // value: optional GigaPath-Flash checkpoint path or URL

    main:
    def embedding_tile_px = 256
    def embedding_mpp = 0.5
    def reuse_grid = tile_px.toString().toInteger() == embedding_tile_px && tile_mpp && tile_mpp.toString().toBigDecimal() == embedding_mpp && tile_overlap.toString().toBigDecimal() == 0
    def ch_embedding_tiles = ch_tiles
    if (!reuse_grid) {
        EMBEDDING_TILETISSUES(ch_slides_with_masks, embedding_tile_px, embedding_mpp, 0, background_fraction)
        ch_embedding_tiles = EMBEDDING_TILETISSUES.out.tiles.map { meta, tiles, _pixels -> [meta, tiles] }
    }

    def sources = [
        resnet50: [checkpoint: resnet50_checkpoint],
        uni: [checkpoint: uni_checkpoint, repository: 'MahmoodLab/UNI', revision: 'b55a5ec6cade1a39edfe6534189a9b8ca7a022f0'],
        uni2: [checkpoint: uni2_checkpoint, repository: 'MahmoodLab/UNI2-h', revision: 'd517a8dd47902dd7c308b3c36f63bce47e7b9a43'],
        'gigapath-flash': [checkpoint: gigapath_flash_checkpoint, repository: 'prov-gigapath/prov-gigapath-flash', revision: '685a3c816fb7bfd4fec697d7f1ed7da57f2e8e86'],
    ]
    def requests = channel.fromList(sources.keySet().findAll { model -> model in models })
        .map { model -> tuple(model, sources[model]) }
        .branch { model, source ->
            provided: source.checkpoint || !source.repository
            return tuple(model, file(source.checkpoint, checkIfExists: true))
            download: true
            return tuple([id: model], 'model', source.repository, source.revision, 'pytorch_model.bin')
        }

    HUGGINGFACE_DOWNLOAD_EMBEDDINGS(requests.download)
    def ch_downloaded = HUGGINGFACE_DOWNLOAD_EMBEDDINGS.out.resource.map { meta, checkpoint -> tuple(meta.id, checkpoint) }
    def ch_checkpoints = requests.provided.mix(ch_downloaded)
    LAZYSLIDE_FEATUREEXTRACTION(ch_slides.join(ch_embedding_tiles).combine(ch_checkpoints), batch_size)

    emit:
    embeddings = LAZYSLIDE_FEATUREEXTRACTION.out.embeddings // channel: [ val(meta), path(embedding_store) ]
}
