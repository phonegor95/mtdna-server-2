nextflow.enable.dsl = 2
include { CALL_MUTECT2 } from '../../workflows/call_mutect2'
include { FILTER_VARIANTS } from '../../modules/local/filter_variants'
workflow {
    CALL_MUTECT2(Channel.fromPath(params.test_bam),
        Channel.value(file(params.test_reference)),
        Channel.value([file(params.test_fai), file(params.test_dict)]),
        Channel.value('chrM'), 'mutect2_fusion')
    FILTER_VARIANTS(CALL_MUTECT2.out.mutect2_ch)
}
