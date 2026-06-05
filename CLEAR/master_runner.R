#!/usr/bin/env Rscript
# ============================================================
# CLEAR/master_runner.R
# Single snATAC SE x GWAS trait pair  -> runs the full CLEAR
# pipeline using the package API (post-refactor names).
#
# Usage:
#   CLI:         Rscript master_runner.R <se_name> <se_path> <trait_name> <trait_path> [gtf_path] [n_perm] [k_lineage] [processed_dir]
#   Interactive: edit USER INPUTS below, then call runCLEAR()
# ============================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(corpcor)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(gridExtra)
  library(ggsci)
  library(rtracklayer)
})

# ---------------------------------------------------------------
# Load the CLEAR package
# ---------------------------------------------------------------
# Prefer an installed package; fall back to devtools::load_all() for
# in-development use from the source directory.
if (requireNamespace("CLEAR", quietly = TRUE)) {
  library(CLEAR)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  pkg_dir <- if (basename(getwd()) == "CLEAR") getwd() else file.path(getwd(), "CLEAR")
  if (!dir.exists(pkg_dir)) pkg_dir <- "."
  devtools::load_all(pkg_dir, quiet = TRUE)
} else {
  stop("Install the CLEAR package or install devtools to use load_all().")
}

# ===========================================================
# INPUTS
# ===========================================================
args <- commandArgs(TRUE)

if (length(args) >= 4) {
  se_name       <- args[1]
  se_path       <- args[2]
  trait_name    <- args[3]
  trait_path    <- args[4]
  gtf_path      <- if (length(args) >= 5) args[5] else "/working/lab_jonathb/lefeiW/projects/CLEAR_data/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz"
  n_perm        <- if (length(args) >= 6) as.integer(args[6]) else 1000L
  # k_lineage: number of lineage groups to cut the dendrogram into (user choice)
  k_lineage     <- if (length(args) >= 7) as.integer(args[7]) else 4L
  processed_dir <- if (length(args) >= 8) args[8] else file.path(dirname(dirname(getwd())), "processed_data")
} else {
  # --- Interactive / RStudio defaults: edit these ---
  se_name       <- "breast_full_330k"
  se_path       <- "/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/breast_330017peak.rds"
  trait_name    <- "BCAC_FM"
  trait_path    <- "/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/data/BCAC_FM_GR.rds"
  gtf_path      <- "/working/lab_jonathb/lefeiW/projects/CLEAR_data/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz"
  n_perm        <- 100L
  # k_lineage: number of lineage groups to cut the dendrogram into
  k_lineage     <- 4L
  processed_dir <- file.path(dirname(dirname(getwd())), "processed_data")
}

# ===========================================================
# runCLEAR()
# Writes results to the global environment so se, gwas_mat,
# locus_results, etc. are inspectable afterwards (matches the
# legacy runner's interactive workflow).
# ===========================================================
runCLEAR <- function() {

  message("\n========================================")
  message("CLEAR analysis (locus): ", se_name, " x ", trait_name)
  message("  SE path:       ", se_path)
  message("  Trait path:    ", trait_path)
  message("  GTF path:      ", gtf_path)
  message("  n_perm:        ", n_perm)
  message("  k_lineage:     ", k_lineage)
  message("  processed_dir: ", processed_dir)
  message("========================================\n")

  # Output directories
  results_dir <<- "results"
  plots_dir   <<- "plots"
  data_dir    <<- "data"
  for (d in c(results_dir, plots_dir, data_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

  # ------------------------------------------------------------
  # Step 1 -- Load trait + (cached) processed SE
  # ------------------------------------------------------------
  message("\n--- Step 1: Load trait and processed/raw SE ---")
  snps <<- readRDS(trait_path)
  message("  SNPs: ", length(snps))

  processed_se_path <- se_path
  if (file.exists(processed_se_path)) {
    message("  Found processed SE cache: ", processed_se_path)
    se <<- readRDS(processed_se_path)
  } else {
    message("  No processed SE cache; loading and preparing raw SE")
    se <<- readRDS(se_path)
    se <<- prepareSE(se)
    se <<- addSEMetrics(se)
    saveRDS(se, processed_se_path)
    message("  Saved processed SE cache: ", processed_se_path)
  }

  required_assays <- c("raw", "raw_l2", "raw_l2_rank")
  if (!all(required_assays %in% assayNames(se))) {
    message("  Cached SE missing expected assays; refreshing cache")
    se <<- readRDS(se_path)
    se <<- prepareSE(se)
    se <<- addSEMetrics(se)
    saveRDS(se, processed_se_path)
  }

  message("  SE: ", nrow(se), " peaks x ", ncol(se), " cell types")

  message("  Building lineage map and plotting dendrogram")
  lineage_obj <<- addLineageMap(assay(se, "raw_l2"), k_lineage = k_lineage)
  plotLineageDendrogram(lineage_obj, se_name, plots_dir)

  # ------------------------------------------------------------
  # Step 2 -- L2 matrix
  # ------------------------------------------------------------
  message("\n--- Step 2: Extract L2 matrix ---")
  mat_l2 <<- assay(se, "raw_l2")

  # ------------------------------------------------------------
  # Step 3 -- GWAS overlap
  # ------------------------------------------------------------
  message("\n--- Step 3: GWAS overlap ---")
  overlap  <<- addGWASOverlap(se, snps)
  se_gwas  <<- overlap$se_gwas
  gwas_mat <<- overlap$gwas_mat
  message("  GWAS peaks: ", nrow(gwas_mat))
  if (nrow(gwas_mat) > 0) {
    n_sig_overlap <- length(unique(na.omit(unlist(strsplit(as.character(rowData(se_gwas)$signal), ",")))))
    gwas_se_name  <- paste0(se_name, "_", trait_name, "_",
                            nrow(gwas_mat), "peaks_", n_sig_overlap, "signals.rds")
    saveRDS(se_gwas, file.path(results_dir, gwas_se_name))
  }

  # ------------------------------------------------------------
  # Step 4 -- Correlation heatmaps + densities
  # ------------------------------------------------------------
  message("\n--- Step 4: Correlation heatmaps and densities ---")
  plotAssayCorHeatmaps(se, se_name, plots_dir)
  plotAssayDensities(se, se_name, plots_dir, lineage_palette = lineage_obj$palette)
  plotLineageCorHeatmap(mat_l2, lineage_obj$cell_lineage, se_name, plots_dir,
                        lineage_colours = lineage_obj$lineage_colours)
  plotLineageDensities(mat_l2, lineage_obj$cell_lineage, se_name, plots_dir,
                       lineage_colours = lineage_obj$lineage_colours)
  if (nrow(gwas_mat) > 1) {
    plotAssayCorHeatmaps(se_gwas, paste0(se_name, "_", trait_name, "_GWAS"), plots_dir)
    plotAssayDensities(se_gwas, paste0(se_name, "_", trait_name, "_GWAS"), plots_dir,
                       lineage_palette = lineage_obj$palette)
    plotLineageCorHeatmap(gwas_mat, lineage_obj$cell_lineage,
                          paste0(se_name, "_", trait_name, "_GWAS"), plots_dir,
                          lineage_colours = lineage_obj$lineage_colours)
    plotLineageDensities(gwas_mat, lineage_obj$cell_lineage,
                         paste0(se_name, "_", trait_name, "_GWAS"), plots_dir,
                         lineage_colours = lineage_obj$lineage_colours)
  }

  # ------------------------------------------------------------
  # Step 5 -- Input summary
  # ------------------------------------------------------------
  message("\n--- Step 5: Input summary ---")
  getInputSummary(se_name, trait_name, se_path, trait_path,
                  se, snps, overlap, results_dir,
                  cell_lineage = lineage_obj$cell_lineage)

  # ------------------------------------------------------------
  # Step 6 -- L2 specificity summary
  # ------------------------------------------------------------
  message("\n--- Step 6: L2 specificity summary ---")
  high_thresh <- 1 / sqrt(2)               # one cell type / lineage > 50% squared mass
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

  plotCelltypeStackbars(
    se_gwas, value_metric = "raw_l2", cumvar_threshold = 0.8,
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_cum80.pdf")),
    lineage_palette = lineage_obj$palette
  )
  plotLineageStackbars(
    se_gwas, value_metric = "raw_l2",
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_lineage_stackbars.pdf")),
    cell_lineage = lineage_obj$cell_lineage,
    lineage_colours = lineage_obj$lineage_colours
  )
  plotCelltypeStackbarsAll(
    se_gwas, value_metric = "raw_l2",
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_all.pdf")),
    lineage_palette = lineage_obj$palette
  )
  plotLineageStackbarsAll(
    se_gwas, value_metric = "raw_l2",
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_lineage_stackbars_all.pdf")),
    cell_lineage = lineage_obj$cell_lineage,
    lineage_colours = lineage_obj$lineage_colours
  )
  plotCelltypeStackbarsByLocus(
    se_gwas, value_metric = "raw_l2",
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_by_locus.pdf")),
    lineage_palette = lineage_obj$palette
  )
  plotCelltypeStackbarsByLocus(
    se_gwas, value_metric = "raw_l2", multi_peak_only = TRUE,
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_by_locus_multipeak.pdf")),
    lineage_palette = lineage_obj$palette
  )
  plotLineageStackbarsByLocus(
    se_gwas, value_metric = "raw_l2",
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_lineage_stackbars_by_locus.pdf")),
    cell_lineage = lineage_obj$cell_lineage,
    lineage_colours = lineage_obj$lineage_colours
  )
  plotLineageStackbarsByLocus(
    se_gwas, value_metric = "raw_l2", multi_peak_only = TRUE,
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_lineage_stackbars_by_locus_multipeak.pdf")),
    cell_lineage = lineage_obj$cell_lineage,
    lineage_colours = lineage_obj$lineage_colours
  )

  # ------------------------------------------------------------
  # Step 7 -- Dominance plots
  # ------------------------------------------------------------
  message("\n--- Step 7: Dominance plots ---")
  plotDominance(mat_l2, gwas_mat, se_name, trait_name, plots_dir,
                lineage_palette = lineage_obj$palette)
  plotDominanceLineage(mat_l2, gwas_mat, lineage_obj$cell_lineage,
                       se_name, trait_name, plots_dir,
                       lineage_colours = lineage_obj$lineage_colours)

  # ------------------------------------------------------------
  # Steps 8 / 8b / 8c -- Locus- and signal-level coherence (DISABLED)
  # ------------------------------------------------------------
  # The locus-level coherence pathway (addLocusCoherence,
  # addDominantPeakCoherence, addSignalSummary, plotCoherence,
  # plotSignalComparison) has been moved to CLEAR/R/old_02-06/ while
  # we focus on variant-level analysis. Re-source those files and
  # uncomment the block below to re-enable.
  #
  # message("\n--- Step 8: Locus-level permutation concordance ---")
  # locus_results <<- addLocusCoherence(se, se_gwas, mat_l2, gwas_mat,
  #                                     gtf_path, lineage_obj, n_perm, results_dir)
  # if (nrow(locus_results)) plotCoherence(locus_results, se_name, trait_name, plots_dir)
  #
  # message("\n--- Step 8b: Dominant-peaks-only concordance ---")
  # dominant_results <<- addDominantPeakCoherence(
  #   se = se, se_gwas = se_gwas, mat_l2 = mat_l2, gwas_mat = gwas_mat,
  #   gtf_path = gtf_path, lineage_obj = lineage_obj,
  #   specificity_summary_l2 = specificity_summary_l2,
  #   n_perm = n_perm, results_dir = results_dir
  # )
  #
  # message("\n--- Step 8c: Signal-level summary ---")
  # signal_summary <<- addSignalSummary(
  #   specificity_summary_l2 = specificity_summary_l2,
  #   locus_results          = if (nrow(locus_results)) locus_results else NULL,
  #   dominant_results       = dominant_results,
  #   results_dir            = results_dir
  # )
  # if (nrow(signal_summary)) {
  #   plotSignalComparison(
  #     signal_summary            = signal_summary,
  #     locus_results_dom_lineage = dominant_results$lineage,
  #     se_name = se_name, trait_name = trait_name, plots_dir = plots_dir
  #   )
  # }

  # ------------------------------------------------------------
  # Step 9 -- Global rank enrichment (cell + lineage)
  # ------------------------------------------------------------
  message("\n--- Step 9: Global rank enrichment (cell) ---")
  addGlobalEnrichment(se, snps, mat_l2, results_dir, plots_dir, se_name, trait_name,
                      level = "cell",
                      lineage_palette = lineage_obj$palette)

  message("\n--- Step 9b: Global rank enrichment (lineage) ---")
  addGlobalEnrichment(se, snps, mat_l2, results_dir, plots_dir, se_name, trait_name,
                      level = "lineage",
                      cell_lineage = lineage_obj$cell_lineage,
                      lineage_colours = lineage_obj$lineage_colours)

  message("\n========================================")
  message("CLEAR analysis complete: ", se_name, " x ", trait_name)
  message("Results: ", normalizePath(results_dir))
  message("Plots:   ", normalizePath(plots_dir))
  message("========================================\n")
  invisible(NULL)
}

# Auto-run when called from CLI; print instructions interactively
if (length(commandArgs(TRUE)) >= 4) {
  runCLEAR()
} else {
  message("Ready. Edit the USER INPUTS block, then call runCLEAR()")
}
