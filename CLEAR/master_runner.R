#!/usr/bin/env Rscript
# ============================================================
# CLEAR/master_runner.R
# Single snATAC SE x GWAS trait pair  -> runs the full CLEAR
# pipeline using the package API (post-refactor names).
#
# Usage:
#   CLI:         Rscript master_runner.R <se_name> <se_path> <trait_name> <trait_path> [gtf_path] [n_perm] [k_lineage] [processed_dir]
#   Interactive: edit USER INPUTS below, then call runCLEAR()
# Optional environment variables: CLEAR_PER_CELL=no|auto|yes,
# CLEAR_DEPTH=none|colsum|quantile, CLEAR_ACCESS_FILTER=none|q0.1|<absolute>,
# CLEAR_NCELLS_TSV=<file>, CLEAR_WRITE_HTML=1|0.
# The n_perm CLI position is retained for compatibility; no permutations run.
# ============================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
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
  gtf_path      <- if (length(args) >= 5 && nzchar(args[5])) args[5] else NULL
  n_perm        <- if (length(args) >= 6) as.integer(args[6]) else 1000L
  # k_lineage: number of lineage groups to cut the dendrogram into (user choice)
  k_lineage     <- if (length(args) >= 7) suppressWarnings(as.numeric(args[7])) else 4L
  processed_dir <- if (length(args) >= 8) args[8] else file.path(getwd(), "processed_data")
} else {
  # --- Interactive / RStudio defaults: edit these ---
  se_name       <- "tissue"
  se_path       <- "snATAC.rds"
  trait_name    <- "trait"
  trait_path    <- "GWAS.rds"
  gtf_path      <- NULL
  n_perm        <- 100L
  # k_lineage: number of lineage groups to cut the dendrogram into
  k_lineage     <- 4L
  processed_dir <- file.path(getwd(), "processed_data")
}

.accessOptions <- function(x) {
  if (identical(x, "none")) return(list(min_max = NULL, min_max_quantile = NULL))
  quantile <- grepl("^q", x)
  value <- suppressWarnings(as.numeric(sub("^q", "", x)))
  if (length(value) != 1L || !is.finite(value) || value < 0 ||
      (quantile && value >= 1)) stop("access_filter must be none, q<0..1), or a non-negative cutoff.")
  if (quantile) list(min_max = NULL, min_max_quantile = value) else
    list(min_max = value, min_max_quantile = NULL)
}

# ===========================================================
# runCLEAR()
# Writes results to the global environment so se, gwas_mat,
# locus_results, etc. are inspectable afterwards (matches the
# legacy runner's interactive workflow).
# ===========================================================
runCLEAR <- function(per_cell = Sys.getenv("CLEAR_PER_CELL", "no"),
                     depth = Sys.getenv("CLEAR_DEPTH", "none"),
                     n_cells = Sys.getenv("CLEAR_NCELLS_TSV", ""),
                     access_filter = Sys.getenv("CLEAR_ACCESS_FILTER", "none"),
                     write_html = Sys.getenv("CLEAR_WRITE_HTML", "1"),
                     write_plots = TRUE) {
  per_cell <- match.arg(per_cell, c("no", "auto", "yes"))
  depth <- match.arg(depth, c("none", "colsum", "quantile"))
  access <- .accessOptions(access_filter)
  if (identical(n_cells, "")) n_cells <- NULL
  if (length(write_html) != 1L || is.na(write_html) ||
      !write_html %in% c(TRUE, FALSE, "1", "0")) stop("write_html must be TRUE/FALSE or 1/0.")
  if (length(write_plots) != 1L || !is.logical(write_plots) || is.na(write_plots)) {
    stop("write_plots must be TRUE or FALSE.")
  }
  write_html <- write_html %in% c(TRUE, "1")

  message("\n========================================")
  message("CLEAR analysis: ", se_name, " x ", trait_name)
  message("  SE path:       ", se_path)
  message("  Trait path:    ", trait_path)
  message("  GTF path:      ", gtf_path)
  message("  per_cell/depth: ", per_cell, " / ", depth)
  message("  access_filter: ", access_filter)
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

  if (!file.exists(se_path)) stop("SE input does not exist: ", se_path)
  cells_key <- n_cells
  if (is.character(n_cells) && length(n_cells) == 1L) {
    if (!file.exists(n_cells)) stop("Cell-count file does not exist: ", n_cells)
    cells_key <- tools::md5sum(n_cells)
  }
  info <- file.info(se_path)
  cache_key <- list(schema = 1L, source = normalizePath(se_path),
                     size = info$size, modified = info$mtime,
                     per_cell = per_cell, depth = depth, n_cells = cells_key)
  processed_se_path <- file.path(processed_dir,
    paste0(se_name, "_CLEAR_pc", per_cell, "_d", depth, "_processed.rds"))
  if (normalizePath(processed_se_path, mustWork = FALSE) == normalizePath(se_path)) {
    stop("Processed cache must not overwrite the SE input.")
  }
  cached <- if (file.exists(processed_se_path)) readRDS(processed_se_path) else NULL
  cache_ok <- inherits(cached, "SummarizedExperiment") &&
    all(c("raw", "norm", "raw_l2", "raw_l2_rank") %in% assayNames(cached)) &&
    identical(metadata(cached)$CLEAR$cache_key, cache_key)
  if (cache_ok) {
    se <<- cached
    message("  Using processed cache: ", processed_se_path)
  } else {
    source_se <- readRDS(se_path)
    source_assay <- if ("counts" %in% assayNames(source_se)) "counts" else
      if ("raw" %in% assayNames(source_se)) "raw" else assayNames(source_se)[1]
    se <<- addSEMetrics(prepareSE(source_se), assay_name = source_assay,
                        per_cell = per_cell, depth = depth, n_cells = n_cells)
    metadata(se)$CLEAR$cache_key <<- cache_key
    saveRDS(se, processed_se_path)
    message("  Saved processed cache: ", processed_se_path)
  }
  se_all <<- do.call(flagAccessiblePeaks, c(list(se = se), access))
  se <<- filterAccessiblePeaks(se_all)

  message("  SE: ", nrow(se), " peaks x ", ncol(se), " cell types")

  message("  Building lineage map and plotting dendrogram")
  lineage_obj <<- addLineageMap(assay(se, "raw_l2"), k_lineage = k_lineage)
  if (write_plots) plotLineageDendrogram(lineage_obj, se_name, plots_dir)

  # ------------------------------------------------------------
  # Step 2 -- L2 matrix
  # ------------------------------------------------------------
  message("\n--- Step 2: Extract L2 matrix ---")
  mat_l2 <<- assay(se, "raw_l2")

  # ------------------------------------------------------------
  # Step 3 -- GWAS overlap
  # ------------------------------------------------------------
  message("\n--- Step 3: GWAS overlap ---")
  overlap  <<- addGWASOverlap(se, snps, se_all = se_all)
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
  if (write_plots) {
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
  }

  # ------------------------------------------------------------
  # Step 5 -- Input summary
  # ------------------------------------------------------------
  message("\n--- Step 5: Input summary ---")
  input_summary <<- getInputSummary(se_name, trait_name, se_path, trait_path,
                  se, snps, overlap, results_dir,
                  cell_lineage = lineage_obj$cell_lineage)

  # ------------------------------------------------------------
  # Step 6 -- Combined cell/lineage peak table
  # ------------------------------------------------------------
  message("\n--- Step 6: Combined peak summary ---")
  peak_summary <<- addSpecPeakTable(
    gwas_mat, snps, se_gwas, results_dir,
    cell_lineage = lineage_obj$cell_lineage, gtf_path = gtf_path
  )

  if (write_plots && nrow(gwas_mat)) {
  plotSpecificityDensity(mat_l2, gwas_mat, se_name, trait_name, plots_dir)

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

  # ------------------------------------------------------------
  # Step 7 -- Dominance plots
  # ------------------------------------------------------------
  message("\n--- Step 7: Dominance plots ---")
  plotDominance(mat_l2, gwas_mat, se_name, trait_name, plots_dir,
                lineage_palette = lineage_obj$palette)
  plotDominanceLineage(mat_l2, gwas_mat, lineage_obj$cell_lineage,
                       se_name, trait_name, plots_dir,
                       lineage_colours = lineage_obj$lineage_colours)
  }

  # ------------------------------------------------------------
  # Step 8 -- Global rank enrichment (cell + lineage)
  # ------------------------------------------------------------
  message("\n--- Step 8: Global rank enrichment (cell) ---")
  global_enrichment_cell <<- addGlobalEnrichment(
                      se, snps, mat_l2, results_dir,
                      if (write_plots) plots_dir else NULL, se_name, trait_name,
                      level = "cell",
                      lineage_palette = lineage_obj$palette)

  message("\n--- Step 8b: Global rank enrichment (lineage) ---")
  global_enrichment_lineage <<- addGlobalEnrichment(
                      se, snps, mat_l2, results_dir,
                      if (write_plots) plots_dir else NULL, se_name, trait_name,
                      level = "lineage",
                      cell_lineage = lineage_obj$cell_lineage,
                      lineage_colours = lineage_obj$lineage_colours)

  if (write_html) {
    message("\n--- Step 9: Offline HTML report ---")
    report_file <<- writeSpecHtmlReport(
      se_name, trait_name, mat_l2, gwas_mat, se_gwas, snps,
      peak_cell = peak_summary, lineage_obj = lineage_obj,
      global_enrichment_cell = global_enrichment_cell,
      global_enrichment_lineage = global_enrichment_lineage,
      input_summary = input_summary,
      output_file = file.path("report", "index.html"),
      mapping_file = file.path(results_dir, "signal_variant_peak_map.txt"),
      signal_file = file.path(results_dir, "signal_summary.txt"))
  }
  message("\n========================================")
  message("CLEAR analysis complete: ", se_name, " x ", trait_name)
  message("Results: ", normalizePath(results_dir))
  message("Plots:   ", normalizePath(plots_dir))
  message("========================================\n")
  invisible(list(peaks = peak_summary, input = input_summary,
                  report = if (write_html) report_file else NULL))
}

# Auto-run when called from CLI; print instructions interactively
if (length(commandArgs(TRUE)) >= 4) {
  runCLEAR()
} else {
  message("Ready. Edit the USER INPUTS block, then call runCLEAR()")
}
