/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Functional Annotation with rGREAT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RGREAT } from '../../../modules/local/rgreat/main'

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_super_enhancers  // channel: [ meta, super_enhancers ]
    genome              // val: genome name

    main:
    ch_versions = Channel.empty()

    //
    // Run rGREAT functional annotation
    //
    RGREAT (
        ch_super_enhancers,
        genome
    )
    ch_versions = ch_versions.mix(RGREAT.out.versions.first())

    emit:
    go_results = RGREAT.out.go_results  // channel: [ meta, tsv ]
    go_volcano = RGREAT.out.go_volcano  // channel: [ meta, pdf ]
    versions   = ch_versions             // channel: [ path ]
}
