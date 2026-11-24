include { LAZYSLIDE_TILEQC                                    } from '../../../modules/local/lazyslide/tileqc/main'
include { LAZYSLIDE_TILEQC as LAZYSLIDE_TILEQC_FOCUSLITENN    } from '../../../modules/local/lazyslide/tileqc/main'
include { LAZYSLIDE_TILEQC as LAZYSLIDE_TILEQC_PATHPROFILERQC } from '../../../modules/local/lazyslide/tileqc/main'
include { LAZYSLIDE_TILETISSUES as GRANDQC_TILETISSUES        } from '../../../modules/local/lazyslide/tiletissues/main'
include { LAZYSLIDE_GRANDQC                                   } from '../../../modules/local/lazyslide/grandqc/main'

workflow QC {
    take:
    ch_slides_with_masks      // channel: [ val(meta), path(slide), path(tissue_mask) ]
    ch_tiles                  // channel: [ val(meta), path(tiles) ]
    grandqc_variant           // value: optional GrandQC model variant

    main:
    def resolved_grandqc_variant = grandqc_variant ?: '7x'

    def ch_slides = ch_slides_with_masks.map { meta, slide, _mask -> [meta, slide] }
    def ch_slides_with_tiles = ch_slides.join(ch_tiles)
    LAZYSLIDE_TILEQC(ch_slides_with_tiles, 'tileqc')
    LAZYSLIDE_TILEQC_FOCUSLITENN(ch_slides_with_tiles, 'focuslitenn')
    LAZYSLIDE_TILEQC_PATHPROFILERQC(ch_slides_with_tiles, 'pathprofilerqc')
    def ch_tile_qc = LAZYSLIDE_TILEQC.out.tiles.mix(
        LAZYSLIDE_TILEQC_FOCUSLITENN.out.tiles,
        LAZYSLIDE_TILEQC_PATHPROFILERQC.out.tiles,
    )

    def grandqc_mpp = ['5x': 2, '7x': 1.5, '10x': 1][resolved_grandqc_variant]
    GRANDQC_TILETISSUES(ch_slides_with_masks, 512, grandqc_mpp, 0.5, 1)
    def ch_grandqc_tiles = GRANDQC_TILETISSUES.out.tiles.map { meta, tiles, _pixels -> [meta, tiles] }
    LAZYSLIDE_GRANDQC(ch_slides.join(ch_grandqc_tiles), resolved_grandqc_variant)

    def ch_multiqc_files = LAZYSLIDE_TILEQC.out.multiqc_files.mix(
        LAZYSLIDE_TILEQC_FOCUSLITENN.out.multiqc_files,
        LAZYSLIDE_TILEQC_PATHPROFILERQC.out.multiqc_files,
        LAZYSLIDE_GRANDQC.out.multiqc_files,
    )

    emit:
    tile_qc  = ch_tile_qc                      // channel: [ val(meta), path(scored_tiles) ]
    slide_qc = LAZYSLIDE_GRANDQC.out.artifacts // channel: [ val(meta), path(artifacts) ]
    multiqc_files = ch_multiqc_files // channel: [ path(multiqc_file) ]
}
