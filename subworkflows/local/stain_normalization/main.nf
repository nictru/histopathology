include { TORCHSTAIN_MACENKO                             } from '../../../modules/local/torchstain/macenko/main'
include { TORCHSTAIN_REINHARD                            } from '../../../modules/local/torchstain/reinhard/main'
include { TORCHSTAIN_REINHARD as TORCHSTAIN_REINHARD_FIT } from '../../../modules/local/torchstain/reinhard/main'

workflow STAIN_NORMALIZATION {
    take:
    ch_slides    // channel: [ val(meta), path(slide) ]
    method       // value: stain normalization method
    reference_id // value: sample ID for either Reinhard method

    main:
    def ch_reference = channel.empty()
    if (method == 'macenko') {
        TORCHSTAIN_MACENKO(ch_slides)
        ch_normalized = TORCHSTAIN_MACENKO.out.slide
    }
    else if (method in ['reinhard', 'reinhardmodified']) {
        if (!reference_id) {
            error("Set --reinhard_reference to a sample ID when using --stain_normalization ${method}.")
        }
        def modified = method == 'reinhardmodified'
        def ch_reference_slide = ch_slides
            .filter { meta, _slide -> meta.id == reference_id }
            .ifEmpty { error("Reference sample '${reference_id}' was not found in the samplesheet.") }
            .first()

        TORCHSTAIN_REINHARD_FIT(ch_reference_slide, [], modified, true)
        ch_reference = TORCHSTAIN_REINHARD_FIT.out.reference
        def fitted_reference = ch_reference.map { _meta, reference -> reference }.first()
        TORCHSTAIN_REINHARD(ch_slides, fitted_reference, modified, false)
        ch_normalized = TORCHSTAIN_REINHARD.out.slide
    }
    else {
        error("Unsupported stain normalization method: ${method}")
    }

    emit:
    slides    = ch_normalized // channel: [ val(meta), path(slide), path(report) ]
    reference = ch_reference  // channel: [ val(meta), path(reference) ]
}
