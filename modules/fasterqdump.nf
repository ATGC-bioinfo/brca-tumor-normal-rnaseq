process FASTERQDUMP {
    tag "$srr"
    publishDir "${params.outdir}/fastq", mode: 'link'

    input:
    tuple val(srr), val(sample_id), val(type), path(sra_file)

    output:
    tuple val(sample_id), val(type), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: reads

    script:
    /* Use the local .sra from PREFETCH so nothing is downloaded twice */
    """
    fasterq-dump --split-files --threads ${task.cpus} --temp . -O . ./${sra_file}
    pigz -p ${task.cpus} -c ${srr}_1.fastq > ${sample_id}_R1.fastq.gz
    pigz -p ${task.cpus} -c ${srr}_2.fastq > ${sample_id}_R2.fastq.gz
    rm ${srr}_1.fastq ${srr}_2.fastq
    """
}
