include { LAZYSLIDE_TILEQC                                    } from '../../../modules/local/lazyslide/tileqc/main'
include { LAZYSLIDE_TILEQC as LAZYSLIDE_TILEQC_FOCUSLITENN    } from '../../../modules/local/lazyslide/tileqc/main'
include { LAZYSLIDE_TILEQC as LAZYSLIDE_TILEQC_PATHPROFILERQC } from '../../../modules/local/lazyslide/tileqc/main'
include { LAZYSLIDE_TILETISSUES as GRANDQC_TILETISSUES        } from '../../../modules/local/lazyslide/tiletissues/main'
include { LAZYSLIDE_GRANDQC                                   } from '../../../modules/local/lazyslide/grandqc/main'

workflow QC {
    take:
    ch_slides_with_masks      // channel: [ val(meta), path(slide), path(tissue_mask) ]
    ch_tiles                  // channel: [ val(meta), path(tiles) ]
    focuslitenn_checkpoint    // value: FocusLiteNN checkpoint path or URL
    pathprofilerqc_checkpoint // value: PathProfilerQC checkpoint path or URL
    grandqc_checkpoint        // value: GrandQC checkpoint path or URL
    grandqc_variant           // value: optional GrandQC model variant

    main:
    def grandqc_checkpoint_file = file(grandqc_checkpoint, checkIfExists: true)
    def checkpoint_variant = [
        'GrandQC_MPP2_exported.pt2': '5x',
        'GrandQC_MPP15_exported.pt2': '7x',
        'GrandQC_MPP1_exported.pt2': '10x',
    ][grandqc_checkpoint_file.name]
    def resolved_grandqc_variant = grandqc_variant ?: checkpoint_variant
    if (!resolved_grandqc_variant) {
        error("Set --grandqc_variant for the checkpoint ${grandqc_checkpoint_file.name}.")
    }
    if (checkpoint_variant && checkpoint_variant != resolved_grandqc_variant) {
        error("${grandqc_checkpoint_file.name} is a ${checkpoint_variant} checkpoint. Supply a checkpoint for --grandqc_variant ${resolved_grandqc_variant}.")
    }

    def ch_slides = ch_slides_with_masks.map { meta, slide, _mask -> [meta, slide] }
    def ch_slides_with_tiles = ch_slides.join(ch_tiles)
    LAZYSLIDE_TILEQC(ch_slides_with_tiles, 'tileqc', [])
    LAZYSLIDE_TILEQC_FOCUSLITENN(ch_slides_with_tiles, 'focuslitenn', file(focuslitenn_checkpoint, checkIfExists: true))
    LAZYSLIDE_TILEQC_PATHPROFILERQC(ch_slides_with_tiles, 'pathprofilerqc', file(pathprofilerqc_checkpoint, checkIfExists: true))
    def ch_tile_qc = LAZYSLIDE_TILEQC.out.tiles.mix(
        LAZYSLIDE_TILEQC_FOCUSLITENN.out.tiles,
        LAZYSLIDE_TILEQC_PATHPROFILERQC.out.tiles,
    )

    def grandqc_mpp = ['5x': 2, '7x': 1.5, '10x': 1][resolved_grandqc_variant]
    GRANDQC_TILETISSUES(ch_slides_with_masks, 512, grandqc_mpp, 0.5, 1)
    def ch_grandqc_tiles = GRANDQC_TILETISSUES.out.tiles.map { meta, tiles, _pixels -> [meta, tiles] }
    LAZYSLIDE_GRANDQC(ch_slides.join(ch_grandqc_tiles), grandqc_checkpoint_file, resolved_grandqc_variant)

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
