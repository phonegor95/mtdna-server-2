include { MUTECT2 } from '../modules/local/mutect2'
include { PARABRICKS_MUTECT2 } from '../modules/local/parabricks_mutect2'

workflow CALL_MUTECT2 {
    take:
    bams
    reference
    fasta_indexes
    contig
    method

    main:
    if (!(params.mutect2_backend in ['gatk', 'parabricks'])) {
        error "mutect2_backend must be gatk or parabricks"
    }
    if (params.mutect2_backend == 'parabricks') {
        PARABRICKS_MUTECT2(bams, reference, fasta_indexes, contig)
        prepared = PARABRICKS_MUTECT2.out.raw_calls
    } else {
        prepared = bams.map { bam -> tuple(bam, []) }
    }
    MUTECT2(prepared, reference, fasta_indexes, contig, method)

    emit:
    mutect2_ch = MUTECT2.out.mutect2_ch
}
