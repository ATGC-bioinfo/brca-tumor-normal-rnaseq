#!/usr/bin/env Rscript
# Cohort plots and tables from per-patient MAFs (made by snpeff_to_maf.py).
# Usage: Rscript maftools_report.R <outdir> <maf1> <maf2> ...
suppressPackageStartupMessages({
  library(maftools)
  library(nnls)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
mafs <- args[-1]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(outdir, f)

# Well-known breast cancer driver genes (TCGA / Pereira et al. 2016)
breast_drivers <- c("PIK3CA", "TP53", "GATA3", "CDH1", "MAP3K1", "KMT2C", "PTEN", "AKT1",
                    "ERBB2", "ESR1", "NF1", "RB1", "CBFB", "RUNX1", "TBX3", "MAP2K4",
                    "ARID1A", "FOXA1", "SF3B1", "KMT2D", "NCOR1", "BRCA1", "BRCA2")

# COSMIC SBS signatures reported in breast cancer (PCAWG Breast-AdenoCA); used for fitting.
# ponytail: fixed list; for a large cohort use de-novo extraction (extractSignatures) instead.
breast_sigs <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS13", "SBS17a", "SBS17b",
                 "SBS18", "SBS37", "SBS39", "SBS40a", "SBS40b", "SBS40c", "SBS41")

maf <- merge_mafs(mafs = mafs, verbose = FALSE)

# ---- Tables ----
write.csv(getSampleSummary(maf), out("sample_summary.csv"), row.names = FALSE)
write.csv(getGeneSummary(maf), out("gene_summary.csv"), row.names = FALSE)
all_muts <- rbind(maf@data, maf@maf.silent, fill = TRUE)
write.csv(all_muts[Hugo_Symbol %in% breast_drivers], out("breast_driver_mutations.csv"), row.names = FALSE)
# ClinVar hits: (likely) pathogenic germline, (likely) oncogenic, or somatic Tier I/II
write.csv(all_muts[grepl("athogenic", ClinVar_CLNSIG) | grepl("ncogenic", ClinVar_ONC) | grepl("Tier_I", ClinVar_SCI)],
          out("clinvar_significant.csv"), row.names = FALSE)

pdf_png <- function(name, w, h, expr) {
  pdf(out(paste0(name, ".pdf")), width = w, height = h); expr(); dev.off()
  png(out(paste0(name, ".png")), width = w, height = h, units = "in", res = 300); expr(); dev.off()
}

# ---- Plots ----
pdf_png("maf_summary", 10, 7, function() plotmafSummary(maf, addStat = "median", titvRaw = FALSE))
pdf_png("oncoplot_top25", 9, 8, function() oncoplot(maf, top = 25, showTumorSampleBarcodes = TRUE))
hits <- intersect(breast_drivers, maf@data$Hugo_Symbol)
if (length(hits) >= 2) {  # oncoplot needs at least 2 genes
  pdf_png("oncoplot_breast_drivers", 9, 6, function()
    oncoplot(maf, genes = hits, showTumorSampleBarcodes = TRUE, keepGeneOrder = FALSE))
}
# Lollipop: where on the protein each mutated breast driver is hit (domains from maftools' bundled Pfam data)
for (g in intersect(breast_drivers, maf@data$Hugo_Symbol)) {
  pdf_png(paste0("lollipop_", g), 9, 4, function()
    lollipopPlot(maf, gene = g, AACol = "HGVSp_Short", labelPos = "all", showMutationRate = TRUE))
}
titv_res <- titv(maf, plot = FALSE, useSyn = TRUE)
write.csv(titv_res$fraction.contribution, out("mutation_spectrum_fraction.csv"), row.names = FALSE)
pdf_png("mutation_spectrum", 9, 6, function() plotTiTv(titv_res))

# ---- Mutational signatures (96 trinucleotide contexts) ----
tnm <- trinucleotideMatrix(maf, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", useSyn = TRUE,
                           ignoreChr = "chrM")  # mito: no SnpEff annotation, heteroplasmy noise
if (!is.null(tnm$APOBEC_scores)) write.csv(tnm$APOBEC_scores, out("apobec_enrichment.csv"), row.names = FALSE)
counts <- t(tnm$nmf_matrix)                     # 96 contexts x samples
db <- readRDS(system.file("extdata", "SBS_v34_signatures.RDs", package = "maftools"))
sigs <- db$db[rownames(counts), intersect(breast_sigs, colnames(db$db))]

# 96-context profile per sample
profile <- sweep(counts, 2, pmax(colSums(counts), 1), "/")
pdf_png("signature_96_profiles", 10, 2.2 * ncol(profile), function()
  plotSignatures(nmfRes = list(signatures = profile), sig_db = "SBS_v34", show_title = TRUE,
                 title_size = 0.8, yaxisLim = NA))

# Cosine similarity of each sample to every COSMIC signature, and best match
cmp <- compareSignatures(nmfRes = list(signatures = profile), sig_db = "SBS_v34", verbose = FALSE)
write.csv(cmp$cosine_similarities, out("signature_cosine_similarity.csv"))

# Fit each sample as a non-negative mix of the breast cancer signatures (NNLS)
fit <- sapply(colnames(counts), function(s) {
  w <- nnls(as.matrix(sigs), counts[, s])$x
  setNames(w / max(sum(w), 1e-9), colnames(sigs))
})
recon <- sapply(colnames(counts), function(s) {
  r <- as.matrix(sigs) %*% (fit[, s] * sum(counts[, s]))
  sum(r * counts[, s]) / sqrt(sum(r^2) * sum(counts[, s]^2))
})
write.csv(rbind(round(fit, 3), reconstruction_cosine = round(recon, 3), mutations = colSums(counts)),
          out("signature_contributions.csv"))
pdf_png("signature_contributions", 8, 5, function() {
  par(mar = c(5, 4, 2, 8), xpd = TRUE)
  cols <- hcl.colors(nrow(fit), "Set 3")
  barplot(fit, col = cols, ylab = "Fraction of mutations", border = NA,
          main = "COSMIC SBS signature contributions (breast cancer set)")
  legend("topright", inset = c(-0.25, 0), legend = rownames(fit), fill = cols, cex = 0.7, bty = "n")
})

cat("maftools report written to", outdir, "\n")
