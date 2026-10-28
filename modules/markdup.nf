process MARKDUP {
    tag "$sample_id"
    publishDir "${params.outdir}/bam/markdup", mode: 'link'

    input:
    tuple val(sample_id), val(type), path(bam)

    output:
    tuple val(sample_id), val(type), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bai"), emit: bam
    path("${sample_id}.markdup.metrics.txt"), emit: metrics

    script:
    """
    gatk --java-options "-Xmx${(task.memory.toMega() * 0.8) as int}m" MarkDuplicates \
        -I ${bam} -O ${sample_id}.markdup.bam -M ${sample_id}.markdup.metrics.txt \
        --CREATE_INDEX true --TMP_DIR .
    """
}
