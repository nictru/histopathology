include { LAZYSLIDE_FEATUREEXTRACTION as LAZYSLIDE_FEATUREEXTRACTION_RESNET50       } from '../../../modules/local/lazyslide/featureextraction/main'
include { LAZYSLIDE_FEATUREEXTRACTION as LAZYSLIDE_FEATUREEXTRACTION_UNI            } from '../../../modules/local/lazyslide/featureextraction/main'
include { LAZYSLIDE_FEATUREEXTRACTION as LAZYSLIDE_FEATUREEXTRACTION_UNI2           } from '../../../modules/local/lazyslide/featureextraction/main'
include { LAZYSLIDE_FEATUREEXTRACTION as LAZYSLIDE_FEATUREEXTRACTION_GIGAPATH_FLASH } from '../../../modules/local/lazyslide/featureextraction/main'
include { LAZYSLIDE_TILETISSUES as EMBEDDING_TILETISSUES                            } from '../../../modules/local/lazyslide/tiletissues/main'
include { HUGGINGFACE_DOWNLOAD as HUGGINGFACE_DOWNLOAD_UNI                          } from '../../../modules/local/huggingface/download/main'
include { HUGGINGFACE_DOWNLOAD as HUGGINGFACE_DOWNLOAD_UNI2                         } from '../../../modules/local/huggingface/download/main'
include { HUGGINGFACE_DOWNLOAD as HUGGINGFACE_DOWNLOAD_GIGAPATH_FLASH               } from '../../../modules/local/huggingface/download/main'

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

    def ch_inputs = ch_slides.join(ch_embedding_tiles)
    def ch_embeddings = channel.empty()
    if ('resnet50' in models) {
        LAZYSLIDE_FEATUREEXTRACTION_RESNET50(ch_inputs, 'resnet50', file(resnet50_checkpoint, checkIfExists: true), batch_size)
        ch_embeddings = ch_embeddings.mix(LAZYSLIDE_FEATUREEXTRACTION_RESNET50.out.embeddings)
    }
    if ('uni' in models) {
        def checkpoint
        if (uni_checkpoint) {
            checkpoint = file(uni_checkpoint, checkIfExists: true)
        }
        else {
            HUGGINGFACE_DOWNLOAD_UNI('model', 'MahmoodLab/UNI', 'b55a5ec6cade1a39edfe6534189a9b8ca7a022f0', 'pytorch_model.bin')
            checkpoint = HUGGINGFACE_DOWNLOAD_UNI.out.resource
        }
        LAZYSLIDE_FEATUREEXTRACTION_UNI(ch_inputs, 'uni', checkpoint, batch_size)
        ch_embeddings = ch_embeddings.mix(LAZYSLIDE_FEATUREEXTRACTION_UNI.out.embeddings)
    }
    if ('uni2' in models) {
        def checkpoint
        if (uni2_checkpoint) {
            checkpoint = file(uni2_checkpoint, checkIfExists: true)
        }
        else {
            HUGGINGFACE_DOWNLOAD_UNI2('model', 'MahmoodLab/UNI2-h', 'd517a8dd47902dd7c308b3c36f63bce47e7b9a43', 'pytorch_model.bin')
            checkpoint = HUGGINGFACE_DOWNLOAD_UNI2.out.resource
        }
        LAZYSLIDE_FEATUREEXTRACTION_UNI2(ch_inputs, 'uni2', checkpoint, batch_size)
        ch_embeddings = ch_embeddings.mix(LAZYSLIDE_FEATUREEXTRACTION_UNI2.out.embeddings)
    }
    if ('gigapath-flash' in models) {
        def checkpoint
        if (gigapath_flash_checkpoint) {
            checkpoint = file(gigapath_flash_checkpoint, checkIfExists: true)
        }
        else {
            HUGGINGFACE_DOWNLOAD_GIGAPATH_FLASH('model', 'prov-gigapath/prov-gigapath-flash', '685a3c816fb7bfd4fec697d7f1ed7da57f2e8e86', 'pytorch_model.bin')
            checkpoint = HUGGINGFACE_DOWNLOAD_GIGAPATH_FLASH.out.resource
        }
        LAZYSLIDE_FEATUREEXTRACTION_GIGAPATH_FLASH(ch_inputs, 'gigapath-flash', checkpoint, batch_size)
        ch_embeddings = ch_embeddings.mix(LAZYSLIDE_FEATUREEXTRACTION_GIGAPATH_FLASH.out.embeddings)
    }

    emit:
    embeddings = ch_embeddings // channel: [ val(meta), path(embedding_store) ]
}
