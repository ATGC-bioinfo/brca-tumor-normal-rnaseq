process PLOTTING {
    tag "plots"
    publishDir "${params.outdir}/plots", mode: 'link'

    input:
    path(filtered_vcfs)   // all Mutect2 calls with FILTER labels
    path(pass_vcfs)       // final somatic calls
    path(gtf)             // gene strands for the RNA-editing check

    output:
    path("*.{png,pdf}")

    script:
    """
    generate_plots.py \
        --filtered ${filtered_vcfs} --pass ${pass_vcfs} --gtf ${gtf} --outdir .
    """
}
