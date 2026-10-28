process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'link'

    input:
    tuple val(sample_id), val(type), path(r1), path(r2)

    output:
    tuple val(sample_id), val(type), path("${sample_id}_trim_R1.fastq.gz"), path("${sample_id}_trim_R2.fastq.gz"), emit: reads
    path("${sample_id}.fastp.json"), emit: json
    path("${sample_id}.fastp.html"), emit: html

    script:
    """
    fastp -i ${r1} -I ${r2} \
        -o ${sample_id}_trim_R1.fastq.gz -O ${sample_id}_trim_R2.fastq.gz \
        --json ${sample_id}.fastp.json --html ${sample_id}.fastp.html \
        --thread ${task.cpus}
    """
}
