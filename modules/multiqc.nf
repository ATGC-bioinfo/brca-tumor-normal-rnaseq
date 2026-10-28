process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'link'

    input:
    path(fastqc_zips)
    path(fastp_jsons)

    output:
    path("multiqc_report.html")
    path("multiqc_report_data")

    script:
    """
    multiqc . -n multiqc_report.html --no-version-check
    """
}
