options(timeout = 600)
options(download.file.method = "libcurl")

install.packages(
  c("remotes", "BiocManager"),
  repos = "https://cloud.r-project.org"
)

remotes::install_github(
  repo = "Lefei-W/CLEAR_personal",
  subdir = "CLEAR",
  dependencies = TRUE,
  upgrade = "never"
)

BiocManager::install(
  c(
    "SummarizedExperiment",
    "GenomicRanges",
    "S4Vectors",
    "IRanges",
    "ComplexHeatmap",
    "rtracklayer",
    "GenomeInfoDb"
  ),
  ask = FALSE,
  update = FALSE
)
# ============================================================
# 00_load_CLEAR_environment.R
# ============================================================

## Optional: force local user library first
local_lib <- file.path(
  Sys.getenv("USERPROFILE"),
  "R",
  "win-library",
  paste0(R.version$major, ".", R.version$minor)
)

if (dir.exists(local_lib)) {
  .libPaths(unique(c(local_lib, .libPaths())))
}

message("R version:")
message(R.version.string)

message("Library paths:")
print(.libPaths())

## Check CLEAR is installed
if (!requireNamespace("CLEAR", quietly = TRUE)) {
  stop("CLEAR is not installed in the current .libPaths(). Check .libPaths() or reinstall CLEAR.")
}

## Load packages
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
})

message("CLEAR environment loaded successfully.")


se_name    <- "breast_inhouse"
se_path    <- "ArchR_matrix/breast_330017peak.rds"          # SummarizedExperiment, peaks x cell types
trait_name <- "BCAC_5375"
trait_path <- "traits/BCAC_FM_GR.rds"        # GRanges
gtf_path   <- "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
k_lineage  <- 4L                                # set after inspecting the dendrogram
processed_dir <- "processed_data"

results_dir <- "results"
plots_dir   <- "plots"
data_dir    <- "data"
for (d in c(results_dir, plots_dir, data_dir, processed_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# Step A -- Load + cache the SE, then plot the dendrogram only
# ============================================================
inspectLineages <- function() {
  message("\n--- Loading SE / preparing metrics ---")
  if (file.exists(se_path) && any(c("raw", "raw_l2") %in%
                                  assayNames(readRDS(se_path)))) {
    se <<- readRDS(se_path)
  } else {
    se_raw <- readRDS(se_path)
    se_raw <- prepareSE(se_raw)
    se_raw <- addSEMetrics(se_raw)
    saveRDS(se_raw, se_path)
    se <<- se_raw
  }
  
  mat_l2 <<- assay(se, "raw_l2")
  
  # k=2 is fine for the dendrogram itself; only the tree is needed here
  message("--- Plotting lineage dendrogram (k_lineage placeholder = 2) ---")
  lineage_obj_preview <<- addLineageMap(mat_l2, k_lineage = 2L)
  plotLineageDendrogram(lineage_obj_preview, se_name, plots_dir)
  
  message("\nDendrogram written to: ",
          normalizePath(file.path(plots_dir,
                                  paste0(se_name, "_lineage_dendrogram.pdf"))))
  message("Inspect it, then set k_lineage and call runCLEAR().")
  invisible(NULL)
}

# ============================================================
# Step B -- Full variant-level CLEAR analysis
# ============================================================
runCLEAR <- function() {
  
  message("\n========================================")
  message("CLEAR (variant-level): ", se_name, " x ", trait_name)
  message("  k_lineage = ", k_lineage)
  message("========================================\n")
  
  # --- Step 1: SE + trait ---
  if (!exists("se", envir = .GlobalEnv) ||
      !all(c("raw", "raw_l2", "raw_l2_rank") %in% assayNames(se))) {
    message("--- Step 1: load + prepare SE ---")
    se_raw <- readRDS(se_path)
    se_raw <- prepareSE(se_raw)
    se_raw <- addSEMetrics(se_raw)
    saveRDS(se_raw, se_path)
    se <<- se_raw
  }
  snps   <<- readRDS(trait_path)
  mat_l2 <<- assay(se, "raw_l2")
  message("  SE: ", nrow(se), " peaks x ", ncol(se), " cell types")
  message("  SNPs: ", length(snps))
  
  # --- Lineage map (user-defined k) ---
  lineage_obj <<- addLineageMap(mat_l2, k_lineage = k_lineage)
  plotLineageDendrogram(lineage_obj, se_name, plots_dir)
  
  # --- Step 3: GWAS overlap ---
  message("\n--- Step 3: GWAS overlap ---")
  overlap  <<- addGWASOverlap(se, snps)
  se_gwas  <<- overlap$se_gwas
  gwas_mat <<- overlap$gwas_mat
  message("  GWAS peaks: ", nrow(gwas_mat))
  if (nrow(gwas_mat) > 0) {
    n_sig <- length(unique(na.omit(unlist(
      strsplit(as.character(rowData(se_gwas)$signal), ",")))))
    saveRDS(se_gwas, file.path(results_dir,
                               paste0(se_name, "_", trait_name, "_",
                                      nrow(gwas_mat), "peaks_", n_sig, "signals.rds")))
  }
  
  # --- Step 4: heatmaps + densities ---
  message("\n--- Step 4: Correlation heatmaps and densities ---")
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
  
  # --- Step 5: input summary ---
  message("\n--- Step 5: Input summary ---")
  getInputSummary(se_name, trait_name, se_path, trait_path,
                  se, snps, overlap, results_dir,
                  cell_lineage = lineage_obj$cell_lineage)
  
  # --- Step 6: variant-peak specificity ---
  message("\n--- Step 6: L2 specificity summary ---")
  high_thresh <- 1 / sqrt(2)
  mid_thresh  <- 1 / sqrt(ncol(mat_l2))
  
  specificity_summary_l2 <<- addSpecificity(
    gwas_mat, snps, high_thresh, mid_thresh, results_dir,
    cell_lineage = lineage_obj$cell_lineage, se_gwas = se_gwas,
    gtf_path = gtf_path
  )
  lineage_specificity_summary_l2 <<- addLineageSpecificity(
    gwas_mat, snps, results_dir,
    cell_lineage = lineage_obj$cell_lineage,
    se_gwas = se_gwas, gtf_path = gtf_path
  )
  specificity_comparison <<- getSpecificityComparison(
    specificity_summary_l2, lineage_specificity_summary_l2, results_dir
  )
  
  plotSpecificityBar(specificity_summary_l2, se_name, trait_name, plots_dir)
  plotSpecificityDensity(mat_l2, gwas_mat, se_name, trait_name, plots_dir)
  plotSpecificityWeighted(specificity_summary_l2, se_name, trait_name, plots_dir,
                          lineage_palette = lineage_obj$palette)
  plotWeightVsTss(specificity_summary_l2, se_name, trait_name, plots_dir)
  
  plotCelltypeStackbars(se_gwas, value_metric = "raw_l2",
                        cumvar_threshold = 0.8,
                        outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_cum80.pdf")),
                        lineage_palette = lineage_obj$palette)
  plotLineageStackbars(se_gwas, value_metric = "raw_l2",
                       outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_lineage_stackbars.pdf")),
                       cell_lineage = lineage_obj$cell_lineage,
                       lineage_colours = lineage_obj$lineage_colours)
  plotCelltypeStackbarsAll(se_gwas, value_metric = "raw_l2",
                           outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_all.pdf")),
                           lineage_palette = lineage_obj$palette)
  plotLineageStackbarsAll(se_gwas, value_metric = "raw_l2",
                          outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_lineage_stackbars_all.pdf")),
                          cell_lineage = lineage_obj$cell_lineage,
                          lineage_colours = lineage_obj$lineage_colours)
  plotCelltypeStackbarsByLocus(se_gwas, value_metric = "raw_l2",
                               outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_by_locus.pdf")),
                               lineage_palette = lineage_obj$palette)
  plotLineageStackbarsByLocus(se_gwas, value_metric = "raw_l2",
                              outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_lineage_stackbars_by_locus.pdf")),
                              cell_lineage = lineage_obj$cell_lineage,
                              lineage_colours = lineage_obj$lineage_colours)
  
  # --- Step 7: dominance ---
  message("\n--- Step 7: Dominance plots ---")
  plotDominance(mat_l2, gwas_mat, se_name, trait_name, plots_dir,
                lineage_palette = lineage_obj$palette)
  plotDominanceLineage(mat_l2, gwas_mat, lineage_obj$cell_lineage,
                       se_name, trait_name, plots_dir,
                       lineage_colours = lineage_obj$lineage_colours)
  
  # --- Step 9: global enrichment (cell + lineage) ---
  message("\n--- Step 9: Global rank enrichment (cell) ---")
  addGlobalEnrichment(se, snps, mat_l2, results_dir, plots_dir,
                      se_name, trait_name, level = "cell",
                      lineage_palette = lineage_obj$palette)
  
  message("\n--- Step 9b: Global rank enrichment (lineage) ---")
  addGlobalEnrichment(se, snps, mat_l2, results_dir, plots_dir,
                      se_name, trait_name, level = "lineage",
                      cell_lineage = lineage_obj$cell_lineage,
                      lineage_colours = lineage_obj$lineage_colours)
  
  message("\n========================================")
  message("Done. Results: ", normalizePath(results_dir))
  message("       Plots:  ", normalizePath(plots_dir))
  message("========================================\n")
  invisible(NULL)
}

message("Ready. Edit INPUTS, then:")
message("  1) inspectLineages()   # look at the dendrogram in plots/")
message("  2) set k_lineage <- <n>")
message("  3) runCLEAR()")

getwd()
local_work <- file.path(Sys.getenv("USERPROFILE"), "R_install_tmp")
dir.create(local_work, recursive = TRUE, showWarnings = FALSE)
setwd(local_work)
getwd()
setwd('/Users/lefeiW/Desktop/CLEAR_test/')
inspectLineages()
k_lineage <- 4
runCLEAR()

