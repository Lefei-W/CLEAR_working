# Install CLEAR as described in README.md, then edit these input paths.

# ============================================================
# 1. Load environment
# ============================================================

suppressPackageStartupMessages({
  library(CLEAR)
  library(SummarizedExperiment)
})


# ============================================================
# 2. Inputs (edit these)
# ============================================================
# Two inputs snATAC matrix in RSE and GWAS in GR
se_path    <- "ArchR_matrix/breast_330017peak.rds"
trait_path <- "traits/BCAC_FM_GR.rds"

# Custom names for folders
se_name    <- "breast_inhouse"
trait_name <- "BCAC_5375"

# Should be standard at least hg38 
gtf_path   <- NULL  # Optional local GTF for TSS and nearest-gene annotation.


# ============================================================
# Output structure
# ------------------------------------------------------------
# All outputs land under <wd>/<se_name>_x_<trait_name>/ with
# sub-folders:
#   results/         tables, summary files
#   plots/           PDFs
#   data/            generic working data (rarely used)
#   processed_data/  cached processed SE (raw + L2 assays)
#   report/          standalone HTML explorer
# ============================================================
run_tag <- paste0(se_name, "_x_", trait_name)
out_dir <- file.path(getwd(), run_tag)

results_dir   <- file.path(out_dir, "results")
plots_dir     <- file.path(out_dir, "plots")
data_dir      <- file.path(out_dir, "data")
processed_dir <- file.path(out_dir, "processed_data")
for (d in c(out_dir, results_dir, plots_dir, data_dir, processed_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

processed_se_path <- file.path(processed_dir, paste0(se_name, "_processed.rds"))

message("Output folder: ", out_dir)



# ============================================================
# 3. Load + cache SE; compute L2 metrics
# ------------------------------------------------------------
# This tutorial processes the input explicitly. The master runner additionally
# validates and reuses caches keyed by normalisation settings.
# ============================================================

se <- readRDS(se_path)
se <- prepareSE(se)
se <- addSEMetrics(se)
# Opt-in alternative, starting from the unprocessed input:
# se <- addSEMetrics(prepareSE(readRDS(se_path)), per_cell = "yes", depth = "colsum")
se_all <- flagAccessiblePeaks(se)  # No filter by default.
# Opt-in alternative:
# se_all <- flagAccessiblePeaks(se, min_max_quantile = 0.1)
se <- filterAccessiblePeaks(se_all)
saveRDS(se, processed_se_path)

mat_l2 <- assay(se, "raw_l2")
message("SE: ", nrow(se), " peaks x ", ncol(se), " cell types")


# ============================================================
# 4. Inspect the dendrogram and decide k_lineage
# ------------------------------------------------------------
# plotLineageDendrogram now draws dashed cut-height lines for
# k = 3..6 (default `k_preview` arg) on top of the chosen-k
# rectangles. Open the PDF, pick the k that best separates the
# biology, then continue.
# ============================================================

# Pick k after inspecting the plot:
k_lineage  <- 4L
lineage_obj <- addLineageMap(mat_l2, k_lineage = k_lineage)
plotLineageDendrogram(lineage_obj, se_name, plots_dir, k_preview = 3:6)


# ============================================================
# 5. Load trait + GWAS overlap
# ============================================================
snps     <- readRDS(trait_path)
overlap  <- addGWASOverlap(se, snps, se_all = se_all)
se_gwas  <- overlap$se_gwas
gwas_mat <- overlap$gwas_mat

if (nrow(gwas_mat) > 0) {
  n_sig <- length(unique(na.omit(unlist(
    strsplit(as.character(rowData(se_gwas)$signal), ",")))))
  saveRDS(se_gwas, file.path(processed_dir,
                             paste0(se_name, "_", trait_name, "_",
                                    nrow(gwas_mat), "peaks_", n_sig, "signals.rds")))
}

# ============================================================
# 6. Correlation heatmaps + per-cell-type densities
# ============================================================
plotAssayCorHeatmaps(se, se_name, plots_dir)
plotAssayDensities(se, se_name, plots_dir,
                   lineage_palette = lineage_obj$palette)
plotLineageCorHeatmap(mat_l2, lineage_obj$cell_lineage, se_name, plots_dir,
                      lineage_colours = lineage_obj$lineage_colours)
plotLineageDensities(mat_l2, lineage_obj$cell_lineage, se_name, plots_dir,
                     lineage_colours = lineage_obj$lineage_colours)

if (nrow(gwas_mat) > 1) {
  tag <- paste0(se_name, "_", trait_name, "_GWAS")
  plotAssayCorHeatmaps(se_gwas, tag, plots_dir)
  plotAssayDensities(se_gwas, tag, plots_dir,
                     lineage_palette = lineage_obj$palette)
  plotLineageCorHeatmap(gwas_mat, lineage_obj$cell_lineage, tag, plots_dir,
                        lineage_colours = lineage_obj$lineage_colours)
  plotLineageDensities(gwas_mat, lineage_obj$cell_lineage, tag, plots_dir,
                       lineage_colours = lineage_obj$lineage_colours)
}

# 7. One combined peak table. Legacy separate-table APIs remain available.
input_summary <- getInputSummary(
  se_name, trait_name, se_path, trait_path, se, snps, overlap, results_dir,
  cell_lineage = lineage_obj$cell_lineage)
peak_summary <- addSpecPeakTable(
  gwas_mat, snps, se_gwas, results_dir,
  cell_lineage = lineage_obj$cell_lineage, gtf_path = gtf_path)

# 8. Global enrichment on the same retained peak universe.
global_cell <- addGlobalEnrichment(
  se, snps, mat_l2, results_dir, plots_dir, se_name, trait_name,
  level = "cell", lineage_palette = lineage_obj$palette)
global_lineage <- addGlobalEnrichment(
  se, snps, mat_l2, results_dir, plots_dir, se_name, trait_name,
  level = "lineage", cell_lineage = lineage_obj$cell_lineage,
  lineage_colours = lineage_obj$lineage_colours)

# 9. Open this HTML file locally; it has no server or network dependency.
report_file <- writeSpecHtmlReport(
  se_name, trait_name, mat_l2, gwas_mat, se_gwas, snps,
  peak_cell = peak_summary, lineage_obj = lineage_obj,
  global_enrichment_cell = global_cell,
  global_enrichment_lineage = global_lineage,
  input_summary = input_summary,
  output_file = file.path(out_dir, "report", "index.html"),
  mapping_file = file.path(results_dir, "signal_variant_peak_map.txt"),
  signal_file = file.path(results_dir, "signal_summary.txt"))
