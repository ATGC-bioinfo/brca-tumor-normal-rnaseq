process MUTECT2 {
    tag "${pair_id}:${contig}"

    input:
    tuple val(pair_id), val(tumor_id), path(tumor_bam), path(tumor_bai), val(normal_id), path(normal_bam), path(normal_bai), val(contig)
    path(genome_fasta)
    path(genome_fai)
    path(genome_dict)
    path(gnomad)    // [vcf, tbi] or [] when not given

    output:
    tuple val(pair_id), val(tumor_id), val(normal_id), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), path("${prefix}.f1r2.tar.gz"), path("${prefix}.vcf.gz.stats"), emit: raw

    script:
    /* One job per chromosome; safe file name even for odd contig names */
    prefix = "${pair_id}.${contig.replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
    def gnomad_vcf = gnomad.find { f -> !f.name.endsWith('.tbi') && !f.name.endsWith('.idx') }
    def gnomad_arg = gnomad_vcf ? "--germline-resource ${gnomad_vcf}" : ""
    """
    gatk --java-options "-Xmx${(task.memory.toMega() * 0.8) as int}m" Mutect2 \
        -R ${genome_fasta} \
        -L "${contig}" \
        -I ${tumor_bam} -tumor ${tumor_id} \
        -I ${normal_bam} -normal ${normal_id} \
        ${gnomad_arg} \
        --native-pair-hmm-threads ${task.cpus} \
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \
        -O ${prefix}.vcf.gz
    """
}

process MERGE_MUTECT2 {
    tag "${pair_id}"
    publishDir "${params.outdir}/mutect2/raw", mode: 'link', pattern: "${pair_id}.unfiltered.*"

    input:
    tuple val(pair_id), val(tumor_id), val(normal_id), path(vcfs), path(tbis), path(f1r2s), path(stats)

    output:
    tuple val(pair_id), val(tumor_id), val(normal_id), path("${pair_id}.unfiltered.vcf.gz"), path("${pair_id}.unfiltered.vcf.gz.tbi"), path(f1r2s), path("${pair_id}.unfiltered.vcf.gz.stats"), emit: raw

    script:
    """
    gatk MergeVcfs ${vcfs.collect { v -> "-I ${v}" }.join(' ')} -O ${pair_id}.unfiltered.vcf.gz
    gatk MergeMutectStats ${stats.collect { st -> "--stats ${st}" }.join(' ')} -O ${pair_id}.unfiltered.vcf.gz.stats
    """
}

process LEARNORIENTATIONMODEL {
    tag "$pair_id"
    publishDir "${params.outdir}/mutect2/orientation", mode: 'link'

    input:
    tuple val(pair_id), val(tumor_id), val(normal_id), path(vcf), path(tbi), path(f1r2s), path(stats)

    output:
    tuple val(pair_id), path("${pair_id}.read-orientation-model.tar.gz"), emit: artifact

    script:
    """
    gatk LearnReadOrientationModel ${[f1r2s].flatten().collect { f -> "-I ${f}" }.join(' ')} -O ${pair_id}.read-orientation-model.tar.gz
    """
}

process FILTERMUTECTCALLS {
    tag "$pair_id"
    publishDir "${params.outdir}/mutect2/filtered", mode: 'link'

    input:
    tuple val(pair_id), val(tumor_id), val(normal_id), path(vcf), path(tbi), path(f1r2s), path(stats), path(artifact_model)
    path(genome_fasta)
    path(genome_fai)
    path(genome_dict)

    output:
    tuple val(pair_id), path("${pair_id}.somatic.filtered.vcf.gz"), path("${pair_id}.somatic.filtered.vcf.gz.tbi"), emit: filtered

    script:
    /* No contamination table: GATK CalculateContamination gives impossible values on RNA-seq */
    """
    gatk FilterMutectCalls \
        -R ${genome_fasta} \
        -V ${vcf} \
        --orientation-bias-artifact-priors ${artifact_model} \
        --stats ${stats} \
        -O ${pair_id}.somatic.filtered.vcf.gz
    """
}

process PASS_FILTER {
    tag "$pair_id"
    publishDir "${params.outdir}/mutect2/pass", mode: 'link'

    input:
    tuple val(pair_id), val(tumor_id), val(normal_id), path(vcf), path(tbi)
    path(editing_sites)    // REDIportal TABLE1_hg38 or [] when not given

    output:
    tuple val(pair_id), path("${pair_id}.somatic.pass.vcf.gz"), path("${pair_id}.somatic.pass.vcf.gz.tbi"), emit: pass

    script:
    /* Only PASS calls that meet the cutoffs in nextflow.config (params.min_* / max_*) */
    def editing_arg = editing_sites ? "--editing_sites ${editing_sites}" : ""
    """
    pass_filter.py ${vcf} \
        --tumor ${tumor_id} --normal ${normal_id} \
        --min_tumor_dp ${params.min_tumor_dp} --min_normal_dp ${params.min_normal_dp} \
        --min_tumor_alt ${params.min_tumor_alt} --min_tumor_af ${params.min_tumor_af} \
        --max_normal_af ${params.max_normal_af} --max_normal_alt ${params.max_normal_alt} \
        --max_pop_af ${params.max_pop_af} --max_indel_len ${params.max_indel_len} \
        --cluster_dist ${params.cluster_dist} ${editing_arg} \
        | bgzip -c > ${pair_id}.somatic.pass.vcf.gz
    tabix -p vcf ${pair_id}.somatic.pass.vcf.gz
    """
}
