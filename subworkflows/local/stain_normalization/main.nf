include { SQUIDPY_NORMALIZE as SQUIDPY_NORMALIZE_FIT } from '../../../modules/local/squidpy/normalize/main'
include { SQUIDPY_NORMALIZE                          } from '../../../modules/local/squidpy/normalize/main'

workflow STAIN_NORMALIZATION {
    take:
    ch_slides    // channel: [ val(meta), path(slide), path(tissue_mask) ]
    method       // value: stain normalization method
    reference_id // value: reference sample ID

    main:
    if (!reference_id) {
        error("Set --stain_reference to a sample ID from the samplesheet.")
    }
    def ch_reference_slide = ch_slides
        .filter { meta, _slide, _tissue_mask -> meta.id == reference_id }
        .ifEmpty { error("Reference sample '${reference_id}' did not produce an ingested slide and tissue mask.") }
        .first()

    SQUIDPY_NORMALIZE_FIT(ch_reference_slide, [], method, true)
    def fitted_reference = SQUIDPY_NORMALIZE_FIT.out.reference.map { _meta, reference -> reference }
    SQUIDPY_NORMALIZE(ch_slides, fitted_reference, method, false)

    emit:
    slides = SQUIDPY_NORMALIZE.out.slide // channel: [ val(meta), path(slide) ]
    multiqc_files = SQUIDPY_NORMALIZE.out.multiqc_files // channel: [ path(multiqc_file) ]
}
