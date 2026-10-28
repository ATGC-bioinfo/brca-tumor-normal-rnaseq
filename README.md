# Breast Cancer Tumor–Normal Somatic Variant Calling (GATK Mutect2)

A Nextflow (DSL2) pipeline for **RNA-seq** of matched tumor/normal breast cancer samples. It downloads
reads from SRA (or uses local FASTQ), runs read QC and trimming, aligns with the splice-aware aligner **HISAT2** following the
GATK RNA-seq best practice (with HISAT2 in place of STAR), calls somatic SNVs/indels with **GATK Mutect2** in tumor–normal mode, filters
them (GATK filters + project cutoffs + known RNA editing sites), annotates them with **SnpEff** and
**ClinVar**, converts them to **MAF**, and produces a **maftools** cohort report (oncoplot, mutation
spectrum, COSMIC signatures). All tools and databases are free.

> **About the data:** GEO GSE351089, RNA-seq (poly(A)/rRNA-depleted total RNA), Illumina HiSeq 1500,
> 2×~148 bp, ~6–8 Gb per sample. Somatic calling from RNA only sees **expressed** genes, depth reflects
> expression level, and allele-specific expression can shift VAFs. See [Known limitations](#known-limitations).

---

## Contents

1. [Workflow at a glance](#workflow-at-a-glance)
2. [Samples](#samples)
3. [Requirements](#requirements)
4. [How to run](#how-to-run)
5. [Step-by-step: what each process does](#step-by-step-what-each-process-does)
6. [Variant filtering criteria](#variant-filtering-criteria)
7. [Parameters](#parameters)
8. [Outputs](#outputs)
9. [Resources, caching and re-running](#resources-caching-and-re-running)
10. [Known limitations](#known-limitations)
11. [Code layout](#code-layout)
12. [Citations](#citations)

---

## Workflow at a glance

```
                       samplesheet.csv
                             │
          ┌──────────────────┴──────────────────┐
   (download_from_sra=true)             (download_from_sra=false)
   PREFETCH ─► FASTERQDUMP                 local fastq1/fastq2
          └──────────────────┬──────────────────┘
                             │  raw paired FASTQ per sample
                ┌────────────┼────────────┐
              FASTQC       FASTP (trim)    │
                └─────┬──────┘             │
                   MULTIQC                 │
                             ▼
     Reference: SAMTOOLS_FAIDX, GATK_DICT, HISAT2_INDEX, INDEX_VCF (gnomAD)
                (built once only if missing; stored in results/reference)
                             ▼
     HISAT2_ALIGN ────────► MARKDUP ─► SPLITNCIGAR      (per sample)
                             ▼
           pair tumor + normal BAMs by patient_id
                             ▼
     MUTECT2 (one job per contig per pair)  ─►  MERGE_MUTECT2
                             ▼
                  LEARNORIENTATIONMODEL
                             ▼
                    FILTERMUTECTCALLS   (labels every call: PASS or reason)
                             ▼
                       PASS_FILTER      (PASS + cutoffs − RNA editing, long indels, clusters)
                        │         │
                        │      PLOTTING
                        ▼
     PREP_CLINVAR ─► SNPEFF (SnpEff + SnpSift ClinVar) ─► VCF_TO_MAF ─► MAFTOOLS (cohort report)
```

---

## Samples

`samplesheet.csv` has one row per sequencing run:

| Column | Meaning |
| --- | --- |
| `sample_id` | Short name used in all file names and as the BAM read-group `SM` (e.g. `T1`, `N1`) |
| `type` | `tumor` or `normal` — exactly these lowercase words |
| `fastq1`, `fastq2` | Local FASTQ paths; only used when `--download_from_sra false` |
| `srr` | SRA run accession; used when `--download_from_sra true` (default) |
| `patient_id` | Links tumor and normal of the same patient. Each patient needs exactly one tumor and one normal |

Current samples:

| Patient | Normal | Tumor |
| --- | --- | --- |
| P1 | N1 (SRR41219225) | T1 (SRR41219228) |
| P2 | N2 (SRR41219226) | T2 (SRR41219229) |
| P3 | N3 (SRR41219227) | T3 (SRR41219230) |

---

## Requirements

- **Nextflow** ≥ 24.04 (tested with 25.10.3) and **conda** (Miniconda / Miniforge).
- Every tool runs from a pinned conda environment in `envs/`; Nextflow creates them automatically on the first
  run (~15 min), so nothing needs to be installed or activated by hand:

| Env file | Tools | Used by |
| --- | --- | --- |
| `envs/qc.yml` | SRA Toolkit 3.2.1, pigz, FastQC 0.12.1, fastp 0.24.0, MultiQC 1.33 | download + QC steps |
| `envs/hisat2.yml` | HISAT2 2.2.1, samtools 1.21 | `HISAT2_INDEX`, `HISAT2_ALIGN` |
| `envs/somatic.yml` | GATK 4.6.2.0, samtools/htslib 1.24, SnpEff + SnpSift 5.4.0c, Python 3.10 + matplotlib | everything else |
| `envs/maftools.yml` | R 4.5, maftools 2.26.0, BSgenome.Hsapiens.UCSC.hg38, nnls | `MAFTOOLS` (kept separate so R cannot clash with GATK) |

- **Hardware:** a 32 GB RAM workstation is enough (HISAT2 needs ~8 GB to build its index, ~5 GB to align).
  Disk: ~25 GB for references + ~150 GB for a full run of the 6 samples.

### Reference files (all free)

Download once into `refs/` inside this folder (the pipeline's default location; `--refs /other/folder`
points it elsewhere). File names must match the ones below. The FASTA, GTF, ClinVar and REDIportal files must all use
GRCh38 / hg38 coordinates; ClinVar's `1, 2, …, MT` contig names are converted to `chr` names by the pipeline.

```bash
mkdir -p refs && cd refs

# 1. Genome: UCSC hg38, primary chromosomes only (chr1-22, X, Y, M; no ALT/random contigs)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip hg38.fa.gz
samtools faidx hg38.fa
samtools faidx hg38.fa $(printf 'chr%s ' {1..22} X Y M) > genome.fa && rm hg38.fa hg38.fa.fai

# 2. Gene annotation: GENCODE v43, same chromosomes as the FASTA
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz \
  | gunzip | awk '/^#/ || $1 ~ /^chr([0-9]+|X|Y|M)$/' > gencode.v43.main_chr.gtf

# 3. gnomAD allele frequencies (Mutect2 germline resource), GATK bundle
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz{,.tbi}

# 4. ClinVar GRCh38
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz{,.tbi}

# 5. REDIportal v3 known RNA editing sites (~1.5 GB)
wget https://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v3.txt.gz

# 6. SnpEff database (MANE transcripts), using the pipeline's own env
conda env create -f ../envs/somatic.yml
conda run -n bc-somatic snpEff download -dataDir $PWD/snpeff GRCh38.mane.1.5.ensembl
```

samtools in step 1 comes from the same env (`conda activate bc-somatic`), or any other samtools install.

---

## How to run

With the reference files in `refs/` (see [Reference files](#reference-files-all-free)):

```bash
git clone https://github.com/ATGC-bioinfo/brca-tumor-normal-rnaseq.git && cd brca-tumor-normal-rnaseq
nextflow run main.nf -profile standard -resume
```

That is the same as writing out every default location:

```bash
nextflow run main.nf -profile standard \
  --samplesheet   samplesheet.csv \
  --genome_fasta  refs/genome.fa \
  --gtf           refs/gencode.v43.main_chr.gtf \
  --gnomad_af     refs/af-only-gnomad.hg38.vcf.gz \
  --editing_sites refs/TABLE1_hg38_v3.txt.gz \
  --snpeff_cache  refs/snpeff \
  --clinvar       refs/clinvar.vcf.gz \
  --outdir        results \
  -resume
```

Files kept elsewhere: give that one parameter (e.g. `--genome_fasta /data/hg38/genome.fa`), or `--refs /data/refs`
if all of them sit in one folder with the names above.

- By default the reads are downloaded from SRA using the `srr` column of `samplesheet.csv`. If the FASTQs
  are already on disk, put their paths in `fastq1`/`fastq2` (relative to the folder you run from, e.g.
  `data/fastq/`) and add `--download_from_sra false`.
- The first run builds the HISAT2 index (~1–1.5 h) into `results/reference/`; later runs reuse it.
- `-resume` reuses every step that already finished (see [caching](#resources-caching-and-re-running)).
- Run in the background: `nextflow run ... -resume -bg > run.log 2>&1`.
- On a SLURM cluster use `-profile slurm`.

---

## Step-by-step: what each process does

Each step lists **tool**, **what it does / why**, **key settings**, and **output** (relative to `results/`).
Code is in `modules/*.nf`; the wiring between steps is in `main.nf`.

### 0. Reference preparation — `modules/reference.nf`

`main.nf` first looks for each index **next to the FASTA**. Only if it is missing does it run the
build process. Built indexes are kept with `storeDir` in `results/reference/`, so they are made once
and reused by every later run.

| Process | Tool | Produces | Needed by |
| --- | --- | --- | --- |
| `SAMTOOLS_FAIDX` | `samtools faidx` | `genome.fa.fai` (contig names + lengths) | GATK tools; also the list of contigs Mutect2 is split over |
| `GATK_DICT` | `gatk CreateSequenceDictionary` | `genome.dict` | All GATK tools |
| `HISAT2_INDEX` (`modules/hisat2.nf`) | `hisat2-build <fasta>` + `hisat2_extract_splice_sites.py <gtf>` | `hisat2_index/`, `splicesites.txt` (~1–1.5 h, ~8 GB RAM) | `HISAT2_ALIGN` |
| `INDEX_VCF` | `gatk IndexFeatureFile` | `.tbi` for the gnomAD VCF (only if missing) | `MUTECT2` |

### 1. Download from SRA — `PREFETCH` (`modules/prefetch.nf`)

- **Tool:** `prefetch <SRR>` (SRA Toolkit).
- **Why:** downloads the compressed `.sra` archive once, so extraction does not re-download.
- **Output:** kept in `work/` only (the `.sra` archive is not needed after extraction).
- Skipped when `--download_from_sra false`.

### 2. Convert to FASTQ — `FASTERQDUMP` (`modules/fasterqdump.nf`)

- **Tool:** `fasterq-dump --split-files` on the local `.sra`, then `pigz` to gzip.
- **Why:** produces paired-end R1/R2 FASTQ; files are renamed from the SRR to the `sample_id`.
- **Output:** `fastq/<sample>_R1.fastq.gz`, `fastq/<sample>_R2.fastq.gz`

### 3. Raw read QC — `FASTQC` (`modules/fastqc.nf`)

- **Tool:** `fastqc` on the **raw** (untrimmed) reads.
- **Why:** per-base quality, GC content, adapter content, duplication — to spot bad libraries early.
  (High duplication is expected for RNA-seq because highly expressed genes are sequenced many times.)
- **Output:** `fastqc/*_fastqc.html`, `*_fastqc.zip`

### 4. Trimming and filtering — `FASTP` (`modules/fastp.nf`)

- **Tool:** `fastp` with **default settings**, paired-end mode.
- **What the defaults do:** adapter detection by read-pair overlap and trimming; discard reads where
  > 40% of bases are below Q15; discard reads shorter than 15 bp after trimming; discard reads with > 5 N bases.
  No fixed-length or sliding-window trimming is applied.
- **Output:** `fastp/<sample>_trim_R{1,2}.fastq.gz` (used for alignment), `<sample>.fastp.json/html` (QC report).

### 5. QC summary — `MULTIQC` (`modules/multiqc.nf`)

- **Tool:** `multiqc` over all FastQC zips and fastp JSONs.
- **Output:** `multiqc/multiqc_report.html` — one page comparing all samples.

### 6. Spliced alignment — `HISAT2_ALIGN` (`modules/hisat2.nf`)

- **Tool:** `hisat2 -x hisat2_index/genome -1 R1 -2 R2 --known-splicesite-infile splicesites.txt --rg-id <s> --rg SM:<s> --rg PL:ILLUMINA --rg LB:<s> | samtools sort`
- **Why a splice-aware aligner:** RNA reads often span an exon–exon junction (the intron is spliced out). A DNA aligner
  (BWA) soft-clips or mis-places these reads, which creates false variants near splice sites.
- **Why HISAT2, not STAR:** STAR is the aligner in the GATK RNA-seq best practice, but building its human index
  needs 30+ GB RAM and was repeatedly killed (out of memory) on the 32 GB machine this pipeline was built on. HISAT2 is also splice-aware,
  builds its index in ~8 GB and aligns in ~5 GB. Its output goes through the same MarkDuplicates → SplitNCigarReads steps.
- **Known splice sites:** extracted from the GENCODE GTF and given at alignment time, so HISAT2 knows annotated
  junctions without the very memory-hungry `--ss/--exon` index build; it still finds novel junctions itself.
- **MAPQ:** HISAT2 gives uniquely mapped reads MAPQ 60, which is what GATK expects (no conversion needed).
- **Read group** (`SM` = `sample_id`) is how Mutect2 tells tumor from normal later.
- **Resources:** 12 CPUs, 12 GB.
- **Output:** `bam/hisat2/<sample>.sorted.bam`, `<sample>.hisat2.summary.txt` (check the overall alignment rate ≥ 85 %).

### 7a. Duplicate marking — `MARKDUP` (`modules/markdup.nf`)

- **Tool:** `gatk MarkDuplicates --CREATE_INDEX true`.
- **Why:** flags PCR/optical duplicates (reads with identical start positions and orientation) so they
  are not counted as independent evidence for a variant. Reads are **marked, not removed**; Mutect2
  ignores marked duplicates by default. (Duplication is naturally high in RNA-seq for highly expressed genes.)
- **Output:** `bam/markdup/<sample>.markdup.bam` + `.bai`, `<sample>.markdup.metrics.txt`.

### 7b. Split spliced reads — `SPLITNCIGAR` (`modules/hisat2.nf`)

- **Tool:** `gatk SplitNCigarReads`.
- **Why:** a spliced read has an `N` (skipped intron) in its CIGAR. This step splits it into one piece per
  exon and hard-clips any overhang into the intron, so Mutect2 does not treat intron gaps as deletions
  or call variants in overhangs. Required by the GATK RNA-seq best practice.
- **Output:** `bam/split/<sample>.split.bam` + `.bai`  ← **BAMs used for variant calling**.

**Not performed:** base quality score recalibration (BQSR) — no known-sites resource is used.

### 8. Tumor–normal pairing — `main.nf`

Marked BAMs are joined to `patient_id` from the samplesheet, split into tumor and normal by `type`,
and joined per patient into `(patient, tumor_id, tumor_bam, normal_id, normal_bam)`.

### 9. Somatic calling — `MUTECT2` (`modules/mutect2.nf`)

- **Tool:** `gatk Mutect2 -I tumor.bam -tumor <T> -I normal.bam -normal <N> -L <contig> --germline-resource gnomAD --f1r2-tar-gz …`
- **Why tumor–normal mode:** a variant seen in the tumor but also in the matched normal is germline
  (or an artifact) and gets down-weighted; only tumor-specific events are reported as somatic.
- **Germline resource (gnomAD):** gives the population allele frequency of each site (`POPAF` in INFO).
  Mutect2 uses it as a prior: common population variants are very likely germline.
- **Parallelisation:** one job **per contig per patient** (contigs read from the `.fai`), with 2 CPUs each.
  This is purely for speed; the results are identical to a single whole-genome run.
- **F1R2 counts** (`--f1r2-tar-gz`): per-site read-orientation counts used in step 11.
- **No intervals/BED and no panel of normals** are used — Mutect2 scans every contig.
- **Output (per contig, intermediate in `work/`):** VCF, `.stats`, `.f1r2.tar.gz`.

### 10. Merge per patient — `MERGE_MUTECT2`

- **Tools:** `gatk MergeVcfs` (per-contig VCFs → one VCF) and `gatk MergeMutectStats` (per-contig stats → one).
  The merged stats file is required by FilterMutectCalls to compute the false-discovery-rate threshold.
- **Output:** `mutect2/raw/<patient>.unfiltered.vcf.gz` (+ `.tbi`, `.stats`) — **every candidate call, no filters applied**.

### 11. Read-orientation model — `LEARNORIENTATIONMODEL`

- **Tool:** `gatk LearnReadOrientationModel` on all F1R2 files of a patient.
- **Why:** learns strand/orientation-specific artifacts (e.g. FFPE deamination C>T, oxidative 8-oxoG G>T)
  that appear on only one read orientation. Used as a prior in filtering.
- **Output:** `mutect2/orientation/<patient>.read-orientation-model.tar.gz`

### 12. Contamination estimate — not used

GATK's `GetPileupSummaries` + `CalculateContamination` were tested and removed: on this RNA-seq data they
estimated 60–81 % contamination for every tumor (impossible), because allele-specific expression skews the
allele counts at common SNPs that the method relies on. FilterMutectCalls then flagged 6,000–10,000 calls per
patient as `contamination`. FilterMutectCalls therefore assumes zero cross-sample contamination.

### 13. Mutect2 filtering — `FILTERMUTECTCALLS`

- **Tool:** `gatk FilterMutectCalls -V unfiltered.vcf.gz --stats … --orientation-bias-artifact-priors …`
- **What it does:** computes, for each call, the probability that it is an error, germline, an artifact,
  etc., and picks a threshold that balances sensitivity and precision (optimises F-score by default).
  **It does not remove variants**: every record stays in the file and the `FILTER` column is set to
  `PASS` or to one or more reasons, e.g. `weak_evidence`, `germline`, `strand_bias`, `orientation`,
  `clustered_events`, `multiallelic`, `slippage`, `base_qual`, `map_qual`, `normal_artifact`,
  `haplotype`, `fragment`, `position`.
- GATK defaults are used — no custom thresholds at this step.
- **Output:** `mutect2/filtered/<patient>.somatic.filtered.vcf.gz` (+ `.tbi`)

### 14. Project cutoffs — `PASS_FILTER` (`bin/pass_filter.py`)

- **Tool:** a small Python script (standard library only), then `bgzip` + `tabix`.
- **What it does:** keeps only records with `FILTER == PASS` **and** that meet all cutoffs in
  [Variant filtering criteria](#variant-filtering-criteria), then removes **known RNA editing sites**
  (see below). Tumor and normal columns are found by sample name.
- **RNA editing:** ADAR enzymes change A to inosine in RNA, which the sequencer reads as G (A>G, or T>C on
  the minus strand). The DNA is unchanged, but because editing levels differ between tumor and normal, Mutect2
  reports these as "somatic". Any variant matching a REDIportal site (same chrom, position, ref and edited
  base) is dropped. On this data it removed 44–89 calls per patient.
- **Output:** `mutect2/pass/<patient>.somatic.pass.vcf.gz` (+ `.tbi`) — **the final high-confidence call set**.
  The exact cutoffs used are written into the VCF header as a `##pass_filter=` line; the counts
  (passed cutoffs / removed as editing / removed as clustered / written) are printed in the task log (`.command.err`).

### 15. Interpretation plots — `PLOTTING` (`bin/generate_plots.py`)

- Reads all Mutect2 calls (filtered VCFs), the final pass VCFs, and the GTF (gene strands).
- **Output (`plots/`, each as `.png` + `.pdf`):**

| Plot | Question it answers | How to read it |
| --- | --- | --- |
| `filter_tumor_vs_normal_vaf` | Do the filters keep only tumor-specific, well-supported calls? | Every Mutect2 call, tumor VAF vs normal VAF. Red (final) points should sit at normal VAF ≈ 0. |
| `mrna_strand_spectrum` | Are the calls still contaminated by RNA editing? | Substitutions written on the gene's mRNA strand. ADAR editing makes mRNA **A>G** much taller than **T>C**; DNA mutations show no such A>G excess. |
| `vaf_distribution` | Clonal or subclonal? | Tumor VAF histogram per patient. A peak near purity/2 = mutations in all tumor cells; a low-VAF pile-up = subclones or residual noise. |

### 16. ClinVar preparation — `PREP_CLINVAR` (`modules/annotate.nf`)

- NCBI ClinVar names contigs `1`, `2`, …, `MT`; the reference uses `chr1`, …, `chrM`. This one-time step
  renames them (awk + bgzip + tabix) and stores the result in `results/reference/clinvar/`.

### 17. Functional annotation — `SNPEFF` (`modules/annotate.nf`)

- **Tools:** `snpEff GRCh38.mane.1.5.ensembl` then `SnpSift annotate` with ClinVar.
- **SnpEff** adds an `ANN` field: gene, effect (missense, stop_gained, frameshift, splice, UTR…), impact
  (HIGH / MODERATE / LOW / MODIFIER), transcript, HGVS c. and p. notation. The **MANE** database uses one
  agreed reference transcript per protein-coding gene (NCBI + Ensembl), so each variant gets one clear answer.
- **SnpSift + ClinVar** adds (matched on chrom, position, ref and alt):
  - `CLNSIG` — germline clinical significance (Pathogenic, Likely_pathogenic, Benign, Uncertain…)
  - `ONC` — oncogenicity (Oncogenic, Likely_oncogenic…)
  - `SCI` — somatic clinical impact, AMP/ASCO/CAP tier (e.g. `Tier_I_-_Strong`)
  - `CLNDN`, `ONCDN`, `SCIDN` — disease names; `CLNREVSTAT` — review status; `CLNVC` — variant type; `ALLELEID`
- **Output:** `annotation/vcf/<patient>.snpeff.clinvar.vcf.gz` (+ `.tbi`), `<patient>.snpeff.html` (effect summary),
  `<patient>.snpeff.genes.txt` (counts per gene).
- Mitochondrial (chrM) variants get no gene annotation: the MANE database has no mitochondrial transcripts.

### 18. VCF → MAF — `VCF_TO_MAF` (`bin/snpeff_to_maf.py`)

- Converts each annotated VCF into a **MAF** (Mutation Annotation Format, the TCGA standard table that
  maftools and cBioPortal read). One row per variant; key columns: `Hugo_Symbol`, `Variant_Classification`,
  `HGVSp_Short` (e.g. `p.H1047R`), `t_depth`, `t_alt_count`, `t_vaf`, `ClinVar_*`. `Tumor_Sample_Barcode` = patient ID.
- SnpEff effects are mapped to MAF classes (e.g. `missense_variant` → `Missense_Mutation`,
  `frameshift_variant` → `Frame_Shift_Del/Ins`, `3_prime_UTR_variant` → `3'UTR`); the most severe effect is used.
  Indel coordinates follow the MAF convention (shared anchor base removed, `-` for an empty allele).
- Self-check: `python3 bin/snpeff_to_maf.py --demo`.
- **Output:** `annotation/maf/<patient>.maf`

### 19. Cohort report — `MAFTOOLS` (`bin/maftools_report.R`, env `envs/maftools.yml`)

All patient MAFs are merged, then (output in `annotation/maftools/`, PDF + PNG for every plot):

| Output | What it shows |
| --- | --- |
| `maf_summary` | Variant classes, types, SNV classes, variants per sample, top genes |
| `oncoplot_top25` | Waterfall plot of the 25 most frequently mutated genes (non-silent mutations) |
| `oncoplot_breast_drivers` | Same, limited to known breast cancer drivers (PIK3CA, TP53, GATA3, CDH1, MAP3K1, KMT2C, PTEN, AKT1, ERBB2, ESR1, …); only drawn if ≥ 2 are hit |
| `lollipop_<GENE>` | Position of each mutation on the protein with Pfam domains, one plot per mutated breast driver gene |
| `mutation_spectrum` + `.csv` | Fraction of the six substitution classes (C>T, C>A, …) per sample |
| `signature_96_profiles` | 96-trinucleotide-context mutation profile of each sample |
| `signature_contributions` + `.csv` | Each sample fitted (non-negative least squares) as a mix of COSMIC SBS v3.4 signatures reported in breast cancer (SBS1, 2, 3, 5, 8, 13, 17a/b, 18, 37, 39, 40a/b/c, 41); `reconstruction_cosine` shows how well the mix explains the sample (> 0.9 good) |
| `signature_cosine_similarity.csv` | Similarity of each sample to all 86 COSMIC signatures |
| `apobec_enrichment.csv` | APOBEC (SBS2/SBS13) enrichment test per sample |
| `breast_driver_mutations.csv` | All mutations (incl. silent/UTR) in breast cancer driver genes |
| `clinvar_significant.csv` | Mutations that are Pathogenic/Likely pathogenic, (Likely) oncogenic, or somatic Tier I/II in ClinVar |
| `gene_summary.csv`, `sample_summary.csv` | Counts per gene / per sample by variant class |

Signature analysis excludes chrM. With only 3 samples, signatures are **fitted to known COSMIC
signatures** rather than discovered de novo (de novo extraction needs many samples).

---

## Variant filtering criteria

Applied in `PASS_FILTER`, on top of GATK's own filters. Values live in `nextflow.config` (`params`) and
can be changed on the command line (e.g. `--min_tumor_af 0.10`). A variant is kept only if **all** hold:

| Criterion | VCF field | Default | Rationale |
| --- | --- | --- | --- |
| GATK filter | `FILTER` | `PASS` | Passed all Mutect2 statistical filters |
| Tumor depth | tumor `FORMAT/DP` | ≥ 20 | Enough reads for a stable VAF estimate |
| Normal depth | normal `FORMAT/DP` | ≥ 20 | Enough reads to be confident the variant is absent from normal. RNA: low normal depth means the gene is barely expressed in normal tissue, so a germline variant can look somatic (raised from 10 after the first run) |
| Tumor alt reads | tumor `FORMAT/AD` (alt) | ≥ 5 | 2–3 supporting reads are often sequencing/alignment noise |
| Tumor VAF | tumor `FORMAT/AF` | ≥ 0.05 | Below 5% at this depth calls are mostly artifacts |
| Normal VAF | normal `FORMAT/AF` | ≤ 0.02 **or** normal alt reads ≤ 1 | Ensures it is somatic, not germline / tumor-in-normal |
| Population AF | `INFO/POPAF` (−log10 AF) | AF < 0.001 | Common in gnomAD → likely germline |
| RNA editing | REDIportal v3 | not a known site | Known A-to-I editing sites are RNA changes, not DNA mutations |
| Indel length | `REF`/`ALT` length difference | ≤ 20 bp | Longer "deletions" in RNA are usually reads across an unannotated splice junction |
| Clustered calls | distance to the next call | > 10 bp | 2+ calls within 10 bp are typically misalignment or editing clusters; all calls in the cluster are dropped |

Details:
- For multi-allelic sites the strongest alt allele is used for each check (in practice these are almost
  always already filtered as `multiallelic`).
- A missing value (`.`) in any required field fails the variant.
- **RNA-seq caveat:** depth reflects gene expression, not sequencing depth, so lowly expressed genes are
  under-represented, and tumor VAF can be skewed by allele-specific expression.
- Novel editing sites not in REDIportal can still remain. Check: editing is on the transcript strand (A>G on `+` genes,
  T>C on `−` genes). In the first run the remaining A>G/T>C calls were mostly on the *opposite* strand, so they are
  not mainly editing.

---

## Parameters

Set in `nextflow.config`; override with `--name value`.

| Parameter | Default | Description |
| --- | --- | --- |
| `--refs` | `refs/` | Folder holding all reference files below (default locations) |
| `--genome_fasta` | `refs/genome.fa` | Uncompressed reference FASTA |
| `--gtf` | `refs/gencode.v43.main_chr.gtf` | GENCODE GTF; known splice sites for HISAT2 |
| `--samplesheet` | `samplesheet.csv` | Sample table (see [Samples](#samples)) |
| `--gnomad_af` | `refs/af-only-gnomad.hg38.vcf.gz` | gnomAD AF-only VCF; Mutect2 germline resource (`null` to skip) |
| `--download_from_sra` | `true` | `true` = fetch reads by `srr`; `false` = use local `fastq1`/`fastq2` |
| `--outdir` | `results` | Output folder |
| `--min_tumor_dp` … `--max_pop_af`, `--max_indel_len` (20), `--cluster_dist` (10, 0 = off) | see table above | Cutoffs for `PASS_FILTER` |
| `--editing_sites` | `refs/TABLE1_hg38_v3.txt.gz` | REDIportal known RNA editing sites (`null` to skip) |
| `--snpeff_cache` | `refs/snpeff` | Folder holding the SnpEff database |
| `--snpeff_db` | `GRCh38.mane.1.5.ensembl` | SnpEff database name |
| `--clinvar` | `refs/clinvar.vcf.gz` | NCBI ClinVar GRCh38 `clinvar.vcf.gz` (with `.tbi`) |
| `--max_cpus` | all CPUs | Upper limit for any single job and for the whole run |
| `--max_memory` | 90% of RAM | Same, for memory |
| `--max_time` | `72 h` | Upper limit per job |

Profiles: `standard` (this computer), `slurm` (cluster). Both use the conda envs in `envs/`.

---

## Outputs

```
results/
├── fastq/               # downloaded FASTQs (only with --download_from_sra true)
├── reference/            # built once: .fai, .dict, hisat2_index/, splicesites.txt, gnomAD .tbi, clinvar/ (chr-renamed)
├── fastqc/               # FastQC on raw reads
├── fastp/                # trimmed FASTQ + fastp json/html
├── multiqc/              # multiqc_report.html
├── bam/
│   ├── hisat2/           # <sample>.sorted.bam, <sample>.hisat2.summary.txt (mapping stats)
│   ├── markdup/          # <sample>.markdup.bam/.bai + duplication metrics
│   └── split/            # <sample>.split.bam/.bai   ← BAMs used for calling
├── mutect2/
│   ├── raw/              # <patient>.unfiltered.vcf.gz (+ .tbi, .stats) — all candidate calls
│   ├── orientation/      # read-orientation model
│   ├── filtered/         # <patient>.somatic.filtered.vcf.gz — ALL calls, FILTER = PASS or reason
│   └── pass/             # <patient>.somatic.pass.vcf.gz — PASS + project cutoffs  ← final call set
├── annotation/
│   ├── vcf/              # <patient>.snpeff.clinvar.vcf.gz, snpeff.html, snpeff.genes.txt
│   ├── maf/              # <patient>.maf
│   └── maftools/         # oncoplots, spectrum, signatures, driver + ClinVar tables
└── plots/                # filter_tumor_vs_normal_vaf, mrna_strand_spectrum, vaf_distribution (.png/.pdf)
```

Outputs in `results/` are **hard links** to the files in `work/` (`publishDir mode: 'link'`): they are normal
files, take no extra disk space, and stay intact if `work/` is deleted.

Useful one-liners:

```bash
# count final variants per patient
for f in results/mutect2/pass/*.vcf.gz; do echo "$f $(zgrep -vc '^#' $f)"; done
# ClinVar-significant / oncogenic hits
column -s, -t < results/annotation/maftools/clinvar_significant.csv | less -S
# why were calls rejected? (FILTER reason counts)
zgrep -v '^#' results/mutect2/filtered/P1.somatic.filtered.vcf.gz | cut -f7 | sort | uniq -c | sort -rn
```

---

## Resources, caching and re-running

- **CPU / RAM** are auto-detected (`--max_cpus`, `--max_memory`). Per-step requests are in the
  `process` block of `nextflow.config`: HISAT2 index 8 CPUs / 12 GB, HISAT2 align 12 CPUs / 12 GB, fastp 8 CPUs,
  MarkDuplicates 2 CPUs / 10 GB, SplitNCigarReads 2 CPUs / 8 GB, Mutect2 2 CPUs / 4 GB per contig job,
  SnpEff 8 GB, most other steps 1 CPU.
- **Retries:** a job killed for memory/time (exit 137, 143, 104, 134, 139) is retried up to 2 times with
  doubled memory; other failures let running jobs finish, then stop the run.
- **`-resume`** reuses finished tasks from the **most recent run in this folder**. Any `nextflow run`
  (even `-preview`) becomes the most recent run, so resume a specific earlier run with its session ID:
  `nextflow log` lists runs; then `-resume <session-id>`.
- `work/` holds all intermediate files; delete it only when you no longer need to resume.
- **Cache mode is `lenient`** (`nextflow.config`): tasks are matched by input path + size, not timestamps,
  so copying or touching input files does not make finished steps run again.

---

## Known limitations

1. **RNA-seq only sees expressed genes**, and depth/VAF depend on expression and allele-specific expression.
   Mutations in silent genes (or on a silenced allele) are invisible.
2. **RNA editing**: known sites are removed, but novel editing sites can remain.
3. **No panel of normals** — recurrent technical artifacts are not removed. GATK's public
   `1000g_pon.hg38.vcf.gz` could be added via `--panel-of-normals`.
4. **No contamination estimate:** GATK CalculateContamination gives impossible values on RNA-seq (see step 12).
5. **No BQSR.**
6. **chrM variants are not annotated** (no mitochondrial transcripts in the MANE SnpEff database) and are left
   out of signature analysis.
7. **3 patients**: cohort statistics (driver significance, de novo signatures) are not meaningful; interpret per patient.
8. **No COSMIC / OncoKB lookup** (licensing); ClinVar `ONC`/`SCI` fields cover well-known cancer hotspots.

---

## Code layout

```
main.nf               # workflow: wiring of all steps, index detection, tumor/normal pairing
nextflow.config       # parameters (filter cutoffs, reference paths), CPU/memory per step, conda envs, profiles
samplesheet.csv       # samples: id, tumor/normal, FASTQ paths or SRA run, patient
refs/                 # reference files (not in git; download commands in README)
envs/                 # pinned conda environments (qc, hisat2, somatic, maftools)
modules/
  prefetch.nf  fasterqdump.nf  fastqc.nf  fastp.nf  multiqc.nf  markdup.nf
  reference.nf        # SAMTOOLS_FAIDX, GATK_DICT, INDEX_VCF
  hisat2.nf           # HISAT2_INDEX, HISAT2_ALIGN, SPLITNCIGAR
  mutect2.nf          # MUTECT2, MERGE_MUTECT2, LEARNORIENTATIONMODEL, FILTERMUTECTCALLS, PASS_FILTER
  annotate.nf         # PREP_CLINVAR, SNPEFF, VCF_TO_MAF, MAFTOOLS
  plotting.nf         # PLOTTING
bin/                  # helper scripts; Nextflow puts bin/ on the PATH of every task
  pass_filter.py      # cutoffs + RNA editing / indel / cluster filters (step 14)
  snpeff_to_maf.py    # SnpEff/ClinVar VCF -> MAF (step 18)
  maftools_report.R   # cohort plots and tables (step 19)
  generate_plots.py   # interpretation plots (step 15)
```

The Python scripts have a built-in self-check: `python3 bin/<script>.py --demo`.

---

## Citations

- Di Tommaso P. et al. *Nextflow enables reproducible computational workflows.* Nat Biotechnol 2017.
- Benjamin D. et al. *Calling somatic SNVs and indels with Mutect2.* bioRxiv 2019. doi:10.1101/861054
- Van der Auwera G. & O'Connor B. *Genomics in the Cloud* (GATK best practices). O'Reilly 2020.
- Kim D. et al. *Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype.* Nature Biotechnology 2019.
- Cingolani P. et al. *A program for annotating and predicting the effects of SNPs, SnpEff.* Fly 2012.
- Cingolani P. et al. *Using Drosophila melanogaster as a model for genotoxic chemical mutational studies (SnpSift).* Front Genet 2012.
- Landrum M. et al. *ClinVar: improvements to accessing data.* Nucleic Acids Res 2020.
- Morales J. et al. *A joint NCBI and EMBL-EBI transcript set for clinical genomics (MANE).* Nature 2022.
- Frankish A. et al. *GENCODE 2023.* Nucleic Acids Res 2023.
- D'Addabbo P. et al. *REDIportal: toward an integrated view of the A-to-I editing.* Nucleic Acids Res 2025.
- Mayakonda A. et al. *Maftools: efficient and comprehensive analysis of somatic variants in cancer.* Genome Res 2018.
- Alexandrov L. et al. *The repertoire of mutational signatures in human cancer (COSMIC SBS).* Nature 2020.
- Danecek P. et al. *Twelve years of SAMtools and BCFtools.* GigaScience 2021.
- Chen S. et al. *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics 2018.
- Andrews S. *FastQC.* Babraham Bioinformatics, 2010.
- Ewels P. et al. *MultiQC.* Bioinformatics 2016.
- Karczewski K. et al. *The mutational constraint spectrum… (gnomAD).* Nature 2020.
- Leinonen R. et al. *The Sequence Read Archive.* Nucleic Acids Res 2011.
