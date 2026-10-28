/*
 * Build reference indexes. main.nf only calls these when the index
 * is not already next to the input file.
 * storeDir keeps them in results/reference: later runs skip the build
 * and reuse the same files, so downstream steps stay cached too.
 */

process SAMTOOLS_FAIDX {
    tag "${fasta.name}"
    storeDir "${params.outdir}/reference"

    input:
    path(fasta)

    output:
    path("${fasta.name}.fai")

    script:
    """
    samtools faidx ${fasta}
    """
}

process GATK_DICT {
    tag "${fasta.name}"
    storeDir "${params.outdir}/reference"

    input:
    path(fasta)

    output:
    path("${fasta.baseName}.dict")

    script:
    """
    gatk CreateSequenceDictionary -R ${fasta} -O ${fasta.baseName}.dict
    """
}

process INDEX_VCF {
    tag "${vcf.name}"
    storeDir "${params.outdir}/reference"

    input:
    path(vcf)

    output:
    path("${vcf.name}.{tbi,idx}")

    script:
    /* .vcf.gz -> .tbi, plain .vcf -> .idx (VCF must be bgzip-compressed, like the GATK files) */
    """
    gatk IndexFeatureFile -I ${vcf}
    """
}
