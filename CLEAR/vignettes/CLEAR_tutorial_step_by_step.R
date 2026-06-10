setwd('/Users/lefeiwang/Desktop/CLEAR_test/')

install.packages(c("remotes", "BiocManager"),
                 repos = "https://cloud.r-project.org")

remotes::install_github(
  repo         = "Lefei-W/CLEAR_working",
  subdir       = "CLEAR",
  dependencies = TRUE,     # CRAN + Bioconductor (auto via biocViews:)
  upgrade      = "never"
)

# ============================================================
# 1. Load environment
# ============================================================

suppressPackageStartupMessages({
  library(CLEAR)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(S4Vectors)
  library(IRanges)
  library(rtracklayer)
  library(ComplexHeatmap)
  library(dplyr)
  library(ggplot2)
  library(GenomeInfoDb)
  library(corpcor)
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
gtf_path   <- "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"


# ============================================================
# Output structure
# ------------------------------------------------------------
# All outputs land under <wd>/<se_name>_x_<trait_name>/ with
# four sub-folders:
#   results/         tables, summary files
#   plots/           PDFs
#   data/            generic working data (rarely used)
#   processed_data/  cached processed SE (raw + L2 assays)
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
# Processed SE (raw + raw_l2 + raw_l2_rank assays) is written
# to processed_data/<se_name>_processed.rds; subsequent runs
# reuse the cache instead of re-normalising.
# ============================================================

CLEAR::prepareSE
CLEAR::addSEMetrics

loadCLEARInputs(se_path, trait_path)


se <- readRDS(se_path)
se <- prepareSE(se)
se <- addSEMetrics(se)
saveRDS(se, processed_se_path)


# Safety net: refresh cache if assays are incomplete
required_assays <- c("raw", "raw_l2", "raw_l2_rank")
if (!all(required_assays %in% assayNames(se))) {
  message("Cached SE missing required assays; refreshing")
  se <- readRDS(se_path)
  se <- prepareSE(se)
  se <- addSEMetrics(se)
  saveRDS(se, processed_se_path)
}

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
overlap  <- addGWASOverlap(se, snps)
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
