process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'link'

    input:
    tuple val(sample_id), val(type), path(r1), path(r2)

    output:
    path("*.zip"), emit: zip
    path("*.html"), emit: html

    script:
    """
    fastqc -t ${task.cpus} ${r1} ${r2}
    """
}
