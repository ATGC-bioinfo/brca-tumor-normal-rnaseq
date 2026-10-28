/*
 * Somatic SNV/indel calling from tumor-normal RNA-seq (breast cancer).
 * FASTQ -> fastp -> HISAT2 -> MarkDuplicates -> SplitNCigarReads -> Mutect2 -> filters
 *       -> SnpEff + ClinVar -> MAF -> maftools report + plots
 * See README.md for every step.
 */

include { PREFETCH } from './modules/prefetch'
include { FASTERQDUMP } from './modules/fasterqdump'
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { MULTIQC } from './modules/multiqc'
include { SAMTOOLS_FAIDX; GATK_DICT; INDEX_VCF } from './modules/reference'
include { HISAT2_INDEX; HISAT2_ALIGN; SPLITNCIGAR } from './modules/hisat2'
include { MARKDUP } from './modules/markdup'
include { MUTECT2; MERGE_MUTECT2; LEARNORIENTATIONMODEL; FILTERMUTECTCALLS; PASS_FILTER } from './modules/mutect2'
include { PREP_CLINVAR; SNPEFF; VCF_TO_MAF; MAFTOOLS } from './modules/annotate'
include { PLOTTING } from './modules/plotting'

/* An index file next to its input (e.g. genome.fa.fai), or null if it is not there */
def existing(path) {
    def f = file(path)
    return f.exists() ? f : null
}

workflow {
    ['genome_fasta', 'gtf', 'snpeff_cache', 'clinvar'].each { p ->
        if (!params[p]) error "Missing required parameter --${p} (see README.md, 'How to run')"
    }
    def fasta = file(params.genome_fasta, checkIfExists: true)
    if (fasta.name.endsWith('.gz')) error "--genome_fasta must be uncompressed (gunzip it first)"
    def gtf = file(params.gtf, checkIfExists: true)

    // ---------------- Reference indexes: reused if next to the FASTA, else built once into results/reference ----------------
    fasta_ch = channel.value(fasta)
    def fai = existing("${fasta}.fai")
    def dict = existing("${fasta.parent}/${fasta.baseName}.dict")
    fai_ch  = fai  ? channel.value(fai)  : SAMTOOLS_FAIDX(fasta_ch)
    dict_ch = dict ? channel.value(dict) : GATK_DICT(fasta_ch)
    hisat2 = HISAT2_INDEX(fasta_ch, gtf)

    // ---------------- Samples ----------------
    rows_ch = channel.fromPath(params.samplesheet, checkIfExists: true).splitCsv(header: true)
    if (params.download_from_sra) {
        fastq_ch = FASTERQDUMP(PREFETCH(rows_ch.map { r -> tuple(r.srr, r.sample_id, r.type) }).sra_dir).reads
    } else {
        fastq_ch = rows_ch.map { r -> tuple(r.sample_id, r.type, file(r.fastq1.trim(), checkIfExists: true), file(r.fastq2.trim(), checkIfExists: true)) }
    }

    // ---------------- QC + trimming ----------------
    fastqc_ch = FASTQC(fastq_ch)
    trimmed_ch = FASTP(fastq_ch)
    MULTIQC(fastqc_ch.zip.collect(), trimmed_ch.json.collect())

    // ---------------- Alignment: HISAT2 -> MarkDuplicates -> SplitNCigarReads ----------------
    bam_aln = HISAT2_ALIGN(trimmed_ch.reads, hisat2.index, hisat2.splicesites)
    markdup_bam = MARKDUP(bam_aln.bam).bam    // sample_id, type, bam, bai
    analysis_bam = SPLITNCIGAR(markdup_bam, fasta_ch, fai_ch, dict_ch).bam

    // ---------------- Pair tumor + normal by patient_id ----------------
    patient_ch = rows_ch.map { r -> tuple(r.sample_id, r.patient_id) }
    bam_pid = analysis_bam.join(patient_ch)    // sid, type, bam, bai, pid
    tumors  = bam_pid.filter { row -> row[1] == 'tumor' }.map { sid, _type, bam, bai, pid -> tuple(pid, sid, bam, bai) }
    normals = bam_pid.filter { row -> row[1] == 'normal' }.map { sid, _type, bam, bai, pid -> tuple(pid, sid, bam, bai) }
    pairs = tumors.join(normals)    // pid, tumor_id, tbam, tbai, normal_id, nbam, nbai
    pair_names = pairs.map { pid, tid, _tbam, _tbai, nid, _nbam, _nbai -> tuple(pid, tid, nid) }

    // ---------------- Mutect2 (one job per chromosome, merged per patient) ----------------
    gnomad_ch = []    // optional germline resource; its index is built if missing
    if (params.gnomad_af) {
        def gnomad = file(params.gnomad_af, checkIfExists: true)
        def gidx = existing("${gnomad}.tbi")
        gnomad_ch = gidx ? channel.value([gnomad, gidx]) : INDEX_VCF(gnomad).map { idx -> [gnomad, idx] }
    }
    contigs = fai_ch.splitCsv(sep: '\t').map { row -> row[0] }
    mutect_parts = MUTECT2(pairs.combine(contigs), fasta_ch, fai_ch, dict_ch, gnomad_ch)
    mutect_raw = MERGE_MUTECT2(mutect_parts.raw.groupTuple(by: [0, 1, 2])).raw
    orient = LEARNORIENTATIONMODEL(mutect_raw).artifact
    filtered = FILTERMUTECTCALLS(mutect_raw.join(orient), fasta_ch, fai_ch, dict_ch).filtered

    // ---------------- Final call set: PASS + project cutoffs - known RNA editing sites ----------------
    editing_ch = params.editing_sites ? channel.value(file(params.editing_sites, checkIfExists: true)) : channel.value([])
    passed = PASS_FILTER(pair_names.join(filtered), editing_ch).pass

    // ---------------- Annotation: SnpEff + ClinVar -> MAF -> maftools ----------------
    clinvar_ch = PREP_CLINVAR(file(params.clinvar, checkIfExists: true))
    annotated = SNPEFF(passed, file(params.snpeff_cache, checkIfExists: true), clinvar_ch).vcf
    MAFTOOLS(VCF_TO_MAF(pair_names.join(annotated)).maf.collect())

    // ---------------- Plots ----------------
    PLOTTING(filtered.map { row -> row[1] }.collect(), passed.map { row -> row[1] }.collect(), gtf)
}
