include { CZI_TOTIFF    } from '../../../modules/local/czi/totiff/main'
include { VIPS_TIFFSAVE } from '../../../modules/local/vips/tiffsave/main'

workflow SLIDE_CONVERT {
    take:
    ch_slides // channel: [ val(meta), path(slide), val(scene), val(mpp) ]

    main:
    def by_format = ch_slides.branch { meta, slide, scene, mpp ->
        def extension = slide.extension.toLowerCase(java.util.Locale.ROOT)
        czi: extension == 'czi'
        return [meta, slide, mpp == null ? '' : mpp, scene == null ? '' : scene]
        tiff_jpeg: extension in ['tif', 'tiff', 'jpg', 'jpeg']
        return [meta, slide, mpp == null ? '' : mpp]
        unsupported: true
        error("Unsupported slide format for ${meta.id}: ${slide.name}")
    }

    CZI_TOTIFF(by_format.czi)
    VIPS_TIFFSAVE(by_format.tiff_jpeg)

    emit:
    slides = CZI_TOTIFF.out.slide.mix(VIPS_TIFFSAVE.out.slide) // channel: [ val(meta), path(slide), path(report) ]
}
