/*
 * RNA-seq alignment (GATK RNA best practice, with HISAT2 as the splice-aware aligner):
 * HISAT2 -> (MarkDuplicates, in markdup.nf) -> SplitNCigarReads
 * HISAT2 instead of STAR: STAR's human index build needs 30+ GB RAM; HISAT2's needs ~8 GB,
 * so the pipeline runs on a 32 GB workstation.
 */

process HISAT2_INDEX {
    tag "${fasta.name}"
    storeDir "${params.outdir}/reference"

    input:
    path(fasta)
    path(gtf)

    output:
    path("hisat2_index"), emit: index
    path("splicesites.txt"), emit: splicesites

    script:
    /* Plain genome index (no --ss/--exon: those make hisat2-build need ~160 GB RAM for human).
       Known splice sites from the GTF are given at alignment time instead, which costs no extra RAM. */
    """
    mkdir hisat2_index
    hisat2-build -p ${task.cpus} ${fasta} hisat2_index/genome
    hisat2_extract_splice_sites.py ${gtf} > splicesites.txt
    """
}

process HISAT2_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/bam/hisat2", mode: 'link', pattern: "*.{bam,summary.txt}"

    input:
    tuple val(sample_id), val(type), path(r1), path(r2)
    path(hisat2_index)
    path(splicesites)

    output:
    tuple val(sample_id), val(type), path("${sample_id}.sorted.bam"), emit: bam
    path("${sample_id}.hisat2.summary.txt"), emit: log

    script:
    /* HISAT2 gives uniquely mapped reads MAPQ 60 (what GATK expects), so no MAPQ fix is needed.
       Read group (SM = sample_id) is required by Mutect2 to tell tumor from normal. */
    """
    hisat2 -p ${task.cpus} -x ${hisat2_index}/genome \
        -1 ${r1} -2 ${r2} \
        --known-splicesite-infile ${splicesites} --no-unal \
        --rg-id ${sample_id} --rg SM:${sample_id} --rg PL:ILLUMINA --rg LB:${sample_id} \
        --summary-file ${sample_id}.hisat2.summary.txt \
    | samtools sort -@ 4 -m 1G -o ${sample_id}.sorted.bam -
    """
}

process SPLITNCIGAR {
    tag "$sample_id"
    publishDir "${params.outdir}/bam/split", mode: 'link'

    input:
    tuple val(sample_id), val(type), path(bam), path(bai)
    path(genome_fasta)
    path(genome_fai)
    path(genome_dict)

    output:
    tuple val(sample_id), val(type), path("${sample_id}.split.bam"), path("${sample_id}.split.bai"), emit: bam

    script:
    /* Splits reads that span introns (N in CIGAR) into exon pieces and hard-clips
       overhangs into introns, so Mutect2 does not see intron gaps as deletions */
    """
    gatk --java-options "-Xmx${(task.memory.toMega() * 0.8) as int}m" SplitNCigarReads \
        -R ${genome_fasta} -I ${bam} -O ${sample_id}.split.bam --tmp-dir .
    """
}
