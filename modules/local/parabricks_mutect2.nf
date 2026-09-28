process PARABRICKS_MUTECT2 {
    tag "${bam_file.baseName}"
    container params.parabricks_container
    cpus params.parabricks_cpus
    memory params.parabricks_memory
    maxForks 1
    containerOptions { workflow.containerEngine == 'docker' ? "--gpus ${params.parabricks_gpus}" : '--nv' }

    input:
    path bam_file
    path reference
    path fasta_index_files
    val detected_contig

    output:
    tuple path(bam_file), path('raw.vcf.gz*'), emit: raw_calls

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def runner = groovy.json.JsonOutput.toJson(params.parabricks_runner)
    """
    # Runtime signature: ${params.parabricks_signature}
    ${quote(params.parabricks_python)} ${quote(moduleDir + '/../../bin/parabricks_mutect.py')} \
        --bam ${quote(bam_file)} --reference ${quote(reference)} \
        --contig ${quote(detected_contig)} --base-quality ${params.baseQ} \
        --gpus ${params.parabricks_gpus} --threads ${task.cpus} \
        --samtools ${quote(params.parabricks_samtools)} --runner-json ${quote(runner)}
    """
}
