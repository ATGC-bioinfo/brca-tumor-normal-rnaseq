/*
 * Annotation of the final (pass) calls:
 * SnpEff (gene / effect / protein change) -> SnpSift + ClinVar -> MAF -> maftools report
 */

process PREP_CLINVAR {
    tag "${vcf.name}"
    storeDir "${params.outdir}/reference/clinvar"

    input:
    path(vcf)

    output:
    tuple path("clinvar.chr.vcf.gz"), path("clinvar.chr.vcf.gz.tbi")

    script:
    /* NCBI ClinVar names contigs 1..22, X, Y, MT; our reference uses chr1.., chrM */
    """
    zcat ${vcf} | awk 'BEGIN{OFS="\\t"}
        /^##contig=<ID=/ { sub(/ID=/, "ID=chr"); sub(/ID=chrMT/, "ID=chrM"); print; next }
        /^#/ { print; next }
        { \$1 = (\$1 == "MT") ? "chrM" : "chr" \$1; print }' | bgzip -c > clinvar.chr.vcf.gz
    tabix -p vcf clinvar.chr.vcf.gz
    """
}

process SNPEFF {
    tag "$pair_id"
    publishDir "${params.outdir}/annotation/vcf", mode: 'link'

    input:
    tuple val(pair_id), path(vcf), path(tbi)
    path(snpeff_cache)
    tuple path(clinvar), path(clinvar_tbi)

    output:
    tuple val(pair_id), path("${pair_id}.snpeff.clinvar.vcf.gz"), path("${pair_id}.snpeff.clinvar.vcf.gz.tbi"), emit: vcf
    path("${pair_id}.snpeff.html"), emit: html
    path("${pair_id}.snpeff.genes.txt"), emit: genes

    script:
    /* ClinVar: CLNSIG = germline significance, ONC = oncogenicity, SCI = somatic clinical impact (AMP/ASCO tier) */
    """
    snpEff -Xmx${(task.memory.toMega() * 0.8) as int}m -dataDir \$PWD/${snpeff_cache} \
        -stats ${pair_id}.snpeff.html ${params.snpeff_db} ${vcf} \
    | SnpSift annotate -info CLNSIG,CLNDN,CLNREVSTAT,CLNVC,ALLELEID,ONC,ONCDN,SCI,SCIDN ${clinvar} - \
    | bgzip -c > ${pair_id}.snpeff.clinvar.vcf.gz
    tabix -p vcf ${pair_id}.snpeff.clinvar.vcf.gz
    """
}

process VCF_TO_MAF {
    tag "$pair_id"
    publishDir "${params.outdir}/annotation/maf", mode: 'link'

    input:
    tuple val(pair_id), val(tumor_id), val(normal_id), path(vcf), path(tbi)

    output:
    path("${pair_id}.maf"), emit: maf

    script:
    """
    snpeff_to_maf.py ${vcf} \
        --tumor ${tumor_id} --normal ${normal_id} --barcode ${pair_id} > ${pair_id}.maf
    """
}

process MAFTOOLS {
    tag "cohort"
    publishDir "${params.outdir}/annotation/maftools", mode: 'link'

    input:
    path(mafs)

    output:
    path("*.{csv,pdf,png}")

    script:
    """
    maftools_report.R . ${mafs}
    """
}
