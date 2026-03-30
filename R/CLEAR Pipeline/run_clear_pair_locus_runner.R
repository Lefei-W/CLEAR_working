#!/usr/bin/env Rscript
# ============================================================
# run_clear_pair_locus.R  (master runner)
# One snATAC SE x GWAS trait pair — sources Functions/ then
# runs the 9-step CLEAR pipeline.
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

# ---------------------------
# Source function modules
# ---------------------------
script_dir <- dirname(path = "/working/lab_jonathb/lefeiW/projects/CLEAR_2026/Functions/")
functions_dir <- file.path(script_dir, "Functions")

for (f in list.files(functions_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}
source('se_preparation.R')
# ===========================================================
# USER INPUTS — edit these before running
# ===========================================================
se_name       <- "breast_full_330k"          # e.g. "Breast 330k"
se_path       <- "/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/breast_330017peak.rds"        # snATAC SE .rds
trait_name    <- "BCAC_FM"            # e.g. "Schizophrenia"
trait_path    <- "/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/BCAC_FineMapping.rds"     # full path to GWAS trait GRanges .rds
gtf_path      <- "/working/lab_jonathb/lefeiW/projects/CLEAR_data/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz"
n_perm        <- 1000L
k_lineage     <- 4L
processed_dir <- file.path(dirname(dirname(getwd())), "processed_data")

# ===========================================================
# run_clear()
# Call this AFTER setting the USER INPUTS above.
# All results are written to the global environment so you
# can inspect se, gwas_mat, locus_results, etc. afterwards.
# ===========================================================
run_clear <- function() {

  message("\n========================================")
  message("CLEAR analysis (locus): ", se_name, " x ", trait_name)
  message("  SE path:    ", se_path)
  message("  Trait path: ", trait_path)
  message("  GTF path:   ", gtf_path)
  message("  n_perm:     ", n_perm)
  message("  k_lineage:  ", k_lineage)
  message("  processed_dir: ", processed_dir)
  message("========================================\n")

  # Output directories
  results_dir    <<- "results"
  plots_dir      <<- "plots"
  data_dir       <<- "data"

  for (d in c(results_dir, plots_dir, data_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Step 1: Load trait and processed/raw SE ---
  processed_se_path <- se_path

  message("\n--- Step 1: Load trait and processed/raw SE ---")
  snps <<- readRDS(trait_path)
  message("  SNPs: ", length(snps))

  if (file.exists(processed_se_path)) {
    message("  Found processed SE cache: ", processed_se_path)
    se <<- readRDS(processed_se_path)
  } else {
    message("  No processed SE cache found; loading and preparing raw SE")
    se <<- readRDS(se_path)
    se <<- prepare_clear_se(se)
    se <<- compute_metrics_se_no_cov(se)
    saveRDS(se, processed_se_path)
    message("  Saved processed SE cache: ", processed_se_path)
  }

  required_assays <- c("raw", "raw_l2", "raw_l2_rank", "cor_raw", "cor_raw_l2")
  if (!all(required_assays %in% assayNames(se))) {
    message("  Cached SE missing expected assays; reloading raw SE and refreshing cache")
    se <<- readRDS(se_path)
    se <<- prepare_clear_se(se)
    se <<- compute_metrics_se_no_cov(se)
    saveRDS(se, processed_se_path)
  }

  message("  SE:   ", nrow(se), " peaks x ", ncol(se), " cell types")

  message("  Building lineage map from genome-wide raw_l2 and plotting dendrogram")
  lineage_obj <<- build_lineage_map(assay(se, "raw_l2"), k_lineage = k_lineage)
  plot_lineage_dendrogram(lineage_obj, se_name, plots_dir)

  # --- Step 2: Extract L2 matrix ---
  message("\n--- Step 2: Extract L2 matrix ---")
  mat_l2 <<- assay(se, "raw_l2")

  # --- Step 3: GWAS overlap ---
  message("\n--- Step 3: GWAS overlap ---")
  overlap <<- compute_gwas_overlap(se, snps)
  se_gwas <<- overlap$se_gwas
  gwas_mat <<- overlap$gwas_mat
  message("  GWAS peaks: ", nrow(gwas_mat))
  if (nrow(gwas_mat) > 0) {
    n_sig_overlap <- length(unique(na.omit(unlist(strsplit(as.character(rowData(se_gwas)$signal), ",")))))
    gwas_se_name <- paste0(se_name, "_", trait_name, "_", nrow(gwas_mat), "peaks_", n_sig_overlap, "signals.rds")
    saveRDS(se_gwas, file.path(results_dir, gwas_se_name))
  }

  # --- Step 4: Correlation heatmaps and densities ---
  message("\n--- Step 4: Correlation heatmaps and densities ---")
  plot_assay_cor_heatmaps(se, se_name, plots_dir)
  plot_assay_densities(se, se_name, plots_dir, lineage_palette = lineage_obj$palette)
  if (nrow(gwas_mat) > 1) {
    plot_assay_cor_heatmaps(se_gwas, paste0(se_name, "_", trait_name, "_GWAS"), plots_dir)
    plot_assay_densities(se_gwas, paste0(se_name, "_", trait_name, "_GWAS"), plots_dir, lineage_palette = lineage_obj$palette)
  }

  # --- Step 5: Input summary ---
  message("\n--- Step 5: Input summary ---")
  summarize_inputs(se_name, trait_name, se_path, trait_path, se, snps, overlap, results_dir)

  # --- Step 6: L2 specificity summary ---
  message("\n--- Step 6: L2 specificity summary ---")
  high_thresh <- 1 / sqrt(3)
  mid_thresh <- 1 / sqrt(ncol(mat_l2))
  specificity_summary_l2 <<- addSpecificity(gwas_mat, snps, high_thresh, mid_thresh, results_dir, cell_lineage = lineage_obj$cell_lineage)
  plot_specificity_bar(specificity_summary_l2, se_name, trait_name, plots_dir)
  plot_celltype_stackbars(
    se_gwas,
    value_metric = "raw_l2",
    cumvar_threshold = 0.5,
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_cum50.pdf")),
    lineage_palette = lineage_obj$palette
  )
  plot_celltype_stackbars(
    se_gwas,
    value_metric = "raw_l2",
    cumvar_threshold = 0.8,
    outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_cum80.pdf")),
    lineage_palette = lineage_obj$palette
  )

  # --- Step 7: Dominance plots ---
  message("\n--- Step 7: Dominance plots ---")
  plot_dominance(mat_l2, gwas_mat, se_name, trait_name, plots_dir, lineage_palette = lineage_obj$palette)
  plot_domiance_lineage(mat_l2, gwas_mat, lineage_obj$cell_lineage, se_name, trait_name, plots_dir)

  # --- Step 8: Locus-level permutation concordance ---
  message("\n--- Step 8: Locus-level permutation concordance ---")
  locus_results <<- addLocusCoherence(se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir)
  if (nrow(locus_results)) {
    plot_coherence(locus_results, se_name, trait_name, plots_dir)
  }

  # --- Step 8b: Locus-level consensus concordance ---
  message("\n--- Step 8b: Locus-level consensus concordance ---")
  locus_results_consensus <<- addLocusCoherence_consensus(
    se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir,
    specificity_summary_l2 = specificity_summary_l2
  )
  if (nrow(locus_results_consensus)) {
    plot_coherence_consensus(locus_results_consensus, se_name, trait_name, plots_dir)
  }

  # --- Step 9: Global rank-based enrichment (L2) ---
  message("\n--- Step 9: Global rank-based enrichment (L2) ---")
  addGlobalEnrichment(se, snps, mat_l2, results_dir, plots_dir, se_name, trait_name)

  message("\n========================================")
  message("CLEAR analysis complete: ", se_name, " x ", trait_name)
  message("Results: ", normalizePath(results_dir))
  message("Plots:   ", normalizePath(plots_dir))
  message("========================================\n")

  invisible(NULL)
}

message("Ready. Set your USER INPUTS above, then call run_clear()")
