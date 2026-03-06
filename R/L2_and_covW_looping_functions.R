# Looping to run compute_metrics_se (with quantile-norm, L2, covariance/correlation weighting) across multiple SummarizedExperiment objects.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Matrix)
  library(preprocessCore)
  library(corpcor)
  library(GenomicRanges)
  library(ggplot2)
  library(tidyr)
})

  # Defaults (override at call-site as needed)
  DEFAULT_RATIO_THRESH        <- 2
  DEFAULT_SPEC_THRESHOLD_SQ   <- 0.5
  DEFAULT_L2_HIGH_THRESHOLD   <- NULL
  DEFAULT_L2_MID_THRESHOLD    <- NULL
  DEFAULT_COV_HIGH_QUANTILE   <- 0.9
  DEFAULT_COV_MID_QUANTILE    <- 0.5

# -------- Input registries (name -> path to RDS) --------
SE_INPUTS <- c(
  union_peak                    = "data/PeakCalls_Xiong/union_GRSE.rds",
  breast_final_union_331832     = "/working/lab_jonathb/lefeiW/projects/CLEAR_breast_demo/results/union_peak331832_CLEAR_se.rds",
  breast_final_union_531279     = "/working/lab_jonathb/lefeiW/projects/CLEAR_breast_demo/results/union_peak531279_CLEAR_se.rds",
  merged                        = "data/merged_peakMatrix_perBP_normalised.rds",
  pbmc10k                       = "/working/joint_projects/bc_risk_locus_multiomics/published_scATAC/processed/pbmc10k_CLEAR_ready.se.rds",
  hten_dcis                     = "/working/joint_projects/bc_risk_locus_multiomics/published_scATAC/processed/hten_dcis.rds",
  zhang2021_mammary             = "output/Zhang2021_mammary.hg38.rds",
  zhang2021_mammary_ArchR       = "/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/output/snATAC_SE_library/zhang2021_mammary/zhang2021_mammary_ArchR_se.rds"
)

TRAIT_INPUTS <- c(
  bcac_overall_opentarget      = "data/opentarget_trait/breast_carcinoma_2017_michailidou_cs210__GCST004988/breast_carcinoma_2017_michailidou_cs210__GCST004988_CLEAR_credible_set_variants.gr.rds",
  bcac_overall_FM_Fachal_2020  = "data/BCAC_FM_GR.rds",
  t1d_santiago                 = "data/opentarget_trait/working/t1d_santiago.rds"
)

safe_l2 <- function(x) {
  s <- sum(x^2)
  if (s == 0) x else x / sqrt(s)
}

add_rank_desc <- function(M) apply(M, 2, function(x) rank(-x, ties.method = "average"))

# Ensure peak IDs exist in rownames / rowData
ensure_peakids <- function(se) {
  rr <- rowRanges(se)
  if (is.null(rownames(se)) || anyNA(rownames(se))) {
    rownames(se) <- paste0(as.character(seqnames(rr)), ":", start(rr), "-", end(rr))
  }
  if (!"peakid" %in% names(rowData(se))) rowData(se)$peakid <- rownames(se)
  se
}

#### Full metrics (L2, cov/cor weighted, QN, ranks) for one SE  (added QN and keep thinking about the composite) ####
compute_metrics_se <- function(se, assay_name = NULL, keep_top_prop = NULL) { # 0.1
  if (is.null(assay_name)) assay_name <- assayNames(se)[1] # should be the raw PeakMatrix assay

  mat <- assay(se, assay_name)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  rn <- rownames(se); if (is.null(rn)) rn <- paste0("peak_", seq_len(nrow(mat)))
  cn <- colnames(se); if (is.null(cn)) cn <- paste0("ct_",   seq_len(ncol(mat)))
  dimnames(mat) <- list(rn, cn)

  row_sums <- rowSums(mat)
  if (!is.null(keep_top_prop)) {
    thr  <- stats::quantile(row_sums, probs = 1 - keep_top_prop, na.rm = TRUE) 
    keep <- row_sums > thr
    se  <- se[keep, , drop = FALSE]
    mat <- mat[keep, , drop = FALSE]
    rn  <- rn[keep]
    dimnames(mat) <- list(rn, cn)
  }

  qn <- preprocessCore::normalize.quantiles(mat)
  dimnames(qn) <- list(rn, cn)

  raw_l2 <- t(apply(mat, 1, safe_l2)); dimnames(raw_l2) <- list(rn, cn)
  qn_l2  <- t(apply(qn,  1, safe_l2)); dimnames(qn_l2)  <- list(rn, cn)

  cov_weight <- function(m) {
    sigma <- corpcor::cov.shrink(m)
    w <- t(solve(sigma) %*% t(m))
    dimnames(w) <- dimnames(m)
    w
  }
  cor_weight <- function(m) {
    R <- corpcor::cor.shrink(m)
    w <- t(solve(R) %*% t(m))
    dimnames(w) <- dimnames(m)
    w
  }

  cov_raw    <- cov_weight(mat)
  cov_qn     <- cov_weight(qn)
  cov_raw_l2 <- cov_weight(raw_l2)
  cov_qn_l2  <- cov_weight(qn_l2)

  cor_raw    <- cor_weight(mat)
  cor_qn     <- cor_weight(qn)
  cor_raw_l2 <- cor_weight(raw_l2)
  cor_qn_l2  <- cor_weight(qn_l2)

  assays(se)[["raw"]]         <- mat
  assays(se)[["qn"]]          <- qn
  assays(se)[["raw_l2"]]      <- raw_l2
  assays(se)[["qn_l2"]]       <- qn_l2
  assays(se)[["raw_l2_rank"]] <- add_rank_desc(raw_l2)
  assays(se)[["qn_l2_rank"]]  <- add_rank_desc(qn_l2)

  assays(se)[["cov_raw"]]     <- cov_raw
  assays(se)[["cov_qn"]]      <- cov_qn
  assays(se)[["cov_raw_l2"]]  <- cov_raw_l2
  assays(se)[["cov_qn_l2"]]   <- cov_qn_l2

  assays(se)[["cor_raw"]]     <- cor_raw
  assays(se)[["cor_qn"]]      <- cor_qn
  assays(se)[["cor_raw_l2"]]  <- cor_raw_l2
  assays(se)[["cor_qn_l2"]]   <- cor_qn_l2

  metadata(se)$CLEAR <- list(
    source_assay   = assay_name,
    keep_top_prop  = keep_top_prop,
    n_peaks        = nrow(se),
    n_celltypes    = ncol(se)
  )

  se
}

process_se_list_metrics <- function(se_list, assay_name = NULL, keep_top_prop = NULL,
                                    out_dir = file.path("data", "processed_se"),
                                    plot_dir = "plots") {
  if (is.null(names(se_list))) names(se_list) <- paste0("se", seq_along(se_list))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  res <- lapply(names(se_list), function(nm) {
    message("\n===== Processing: ", nm, " =====")
    se_proc <- compute_metrics_se(se_list[[nm]], assay_name = assay_name, keep_top_prop = keep_top_prop)
    fname <- paste0(nm, "_", ncol(se_proc), "ct_", nrow(se_proc), "peaks_processed.rds")
    saveRDS(se_proc, file.path(out_dir, fname))
    message("Saved: ", file.path(out_dir, fname))
    # First-line validation plots
    plot_assay_cor_heatmaps(se_proc, nm, plot_dir = plot_dir)
    plot_assay_densities(se_proc, nm, plot_dir = plot_dir)
    se_proc
  })
  names(res) <- names(se_list)
  res
}

overlap_se_with_snps <- function(se, snps, assay_name, tissue, trait,
                                  out_dir = file.path("data", "GWAS_overlapped_se")) {
  if (!assay_name %in% assayNames(se)) stop("Assay ", assay_name, " not found in SE")
  se <- ensure_peakids(se)
  ov <- findOverlaps(se, snps, ignore.strand = TRUE)
  if (!length(ov)) {
    message("No overlaps for ", tissue, " x ", trait)
    return(list(se_overlap = se[0, ], mat = matrix(0, 0, ncol(se))))
  }
  hit_rows <- sort(unique(queryHits(ov)))
  se_sub <- se[hit_rows, , drop = FALSE]

  # Annotate each overlapped peak with its SNP IDs and GWAS signals
  snp_list  <- tapply(mcols(snps)$names[subjectHits(ov)],  queryHits(ov), paste, collapse = ",")
  locus_vec <- tapply(mcols(snps)$signal[subjectHits(ov)], queryHits(ov),
                      function(x) paste(unique(x), collapse = ","))
  rowData(se_sub)$overlap_snps    <- snp_list[as.character(hit_rows)]
  rowData(se_sub)$overlap_signals <- locus_vec[as.character(hit_rows)]

  mat <- SummarizedExperiment::assay(se_sub, assay_name)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fname <- paste0(tissue, "_", ncol(se_sub), "ct_", trait, "_", nrow(se_sub), "peaks_overlapped.rds") # naming: tissue, trait, n peaks, n cell types
  saveRDS(se_sub, file.path(out_dir, fname))
  message("Saved: ", file.path(out_dir, fname))

  list(se_overlap = se_sub, mat = mat)
}

# ---------- Assays to plot (skip ranks) ----------
PLOT_ASSAYS <- c("raw", "qn", "raw_l2", "qn_l2",
                 "cov_raw", "cov_qn", "cov_raw_l2", "cov_qn_l2",
                 "cor_raw", "cor_qn", "cor_raw_l2", "cor_qn_l2")

# Correlation heatmaps — one PDF per assay, saved under plots/{tissue}/
plot_assay_cor_heatmaps <- function(se, tissue, plot_dir = "plots") {
  out <- file.path(plot_dir, tissue)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  avail <- intersect(PLOT_ASSAYS, assayNames(se))
  for (a in avail) {
    mat <- assay(se, a)
    if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
    fname <- file.path(out, paste0(tissue, "_", a, "_cor_heatmap.pdf"))
    pdf(fname, width = 8, height = 7)
    heatmap(cor(mat), main = paste0(tissue, " \u2014 ", a, " correlation"),
            margins = c(10, 10))
    dev.off()
    message("Saved: ", fname)
  }
  invisible(out)
}

# Density plots — one PDF with faceted panels per assay
plot_assay_densities <- function(se, tissue, plot_dir = "plots") {
  out <- file.path(plot_dir, tissue)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  avail <- intersect(PLOT_ASSAYS, assayNames(se))
  long_list <- lapply(avail, function(a) {
    mat <- assay(se, a)
    if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
    df <- as.data.frame(mat)
    df <- tidyr::pivot_longer(df, cols = everything(),
                              names_to = "cell_type", values_to = "value")
    df$assay <- a
    df
  })
  long <- do.call(rbind, long_list)
  long$assay <- factor(long$assay, levels = avail)
  p <- ggplot(long, aes(x = value, colour = cell_type)) +
    geom_density() +
    facet_wrap(~ assay, scales = "free", ncol = 4) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 9)) +
    labs(title = paste0(tissue, " \u2014 per-cell-type distributions"),
         x = "Value", y = "Density", colour = "Cell type")
  fname <- file.path(out, paste0(tissue, "_assay_densities.pdf"))
  ggsave(fname, p, width = 16, height = 12)
  message("Saved: ", fname)
  invisible(fname)
}


# ---------- Assays for dominance / specificity ----------
L2_ASSAYS  <- c("raw_l2", "qn_l2")
COV_ASSAYS <- c("cov_raw", "cov_qn", "cov_raw_l2", "cov_qn_l2")
COR_ASSAYS <- c("cor_raw", "cor_qn", "cor_raw_l2", "cor_qn_l2")

#### dominance of top cell type as 2 times the second ####
compute_dominance <- function(mat, ratio_thresh = DEFAULT_RATIO_THRESH) {
  stopifnot(is.matrix(mat))
  cell_names <- colnames(mat)
  res <- t(apply(mat, 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    max_val <- x[ord[1]]
    second_val <- x[ord[2]]
    ratio <- max_val / second_val
    c(maximum_ct = cell_names[ord[1]], ratio = ratio, is_dominant = ratio > ratio_thresh)
  }))
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  res$ratio <- as.numeric(res$ratio)
  res$is_dominant <- res$is_dominant == "TRUE"
  res$peak <- rownames(mat)
  res
}


compute_specificity_summary <- function(mat_l2, cell_cor = NULL, ratio_thresh = DEFAULT_RATIO_THRESH, threshold_sq = DEFAULT_SPEC_THRESHOLD_SQ) {
  if (is.null(cell_cor)) cell_cor <- cor(mat_l2)
  results <- apply(mat_l2, 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    x_sorted <- x[ord]
    cells_sorted <- colnames(mat_l2)[ord]
    cum_sq <- cumsum(x_sorted^2)
    k <- which(cum_sq >= threshold_sq)[1]
    top_cells <- cells_sorted[seq_len(k)]
    maximum_cell  <- cells_sorted[1]
    maximum_score <- x_sorted[1]
    dominant_ratio <- if (length(x_sorted) > 1) x_sorted[1] / x_sorted[2] else NA
    coherence <- if (length(top_cells) > 1) {
      idx <- match(top_cells, colnames(mat_l2))
      sub_cor <- cell_cor[idx, idx, drop = FALSE]
      mean(sub_cor[upper.tri(sub_cor)])
    } else {
      NA
    }
    c(
      maximum_cell  = maximum_cell,
      maximum_score = maximum_score,
      dominant_ratio = dominant_ratio,
      is_dominant = dominant_ratio > ratio_thresh,
      k_cells        = k,
      top_k_cells    = paste(top_cells, collapse = ","),
      coherence      = coherence
    )
  })
  out <- as.data.frame(t(results), stringsAsFactors = FALSE)
  out$peak <- rownames(mat_l2)
  rownames(out) <- NULL
  out
}

specificity_counts_l2 <- function(mat_l2, l2_high_thresh = DEFAULT_L2_HIGH_THRESHOLD, l2_mid_thresh = DEFAULT_L2_MID_THRESHOLD) {
  high_thresh <- l2_high_thresh
  mid_thresh  <- l2_mid_thresh
  if (is.null(high_thresh)) high_thresh <- 1 / sqrt(3)
  if (is.null(mid_thresh))  mid_thresh  <- 1 / sqrt(ncol(mat_l2))
  data.frame(
    peak   = rownames(mat_l2),
    n_high = rowSums(mat_l2 > high_thresh),
    n_mid  = rowSums(mat_l2 > mid_thresh & mat_l2 <= high_thresh),
    n_low  = rowSums(mat_l2 <= mid_thresh)
  )
}

specificity_counts_quantile <- function(mat, cov_high_quantile = DEFAULT_COV_HIGH_QUANTILE, cov_mid_quantile = DEFAULT_COV_MID_QUANTILE) {
  high_thresh <- as.numeric(stats::quantile(mat, probs = cov_high_quantile, na.rm = TRUE))
  mid_thresh  <- as.numeric(stats::quantile(mat, probs = cov_mid_quantile,  na.rm = TRUE))
  data.frame(
    peak   = rownames(mat),
    n_high = rowSums(mat > high_thresh),
    n_mid  = rowSums(mat > mid_thresh & mat <= high_thresh),
    n_low  = rowSums(mat <= mid_thresh),
    high_thresh = high_thresh,
    mid_thresh  = mid_thresh
  )
}


# ---- Run dominance / specificity on one SE across its assays, save to out_dir ----
compute_se_summaries <- function(se, out_dir,
                                 l2_assays  = L2_ASSAYS,
                                 cov_assays = c(COV_ASSAYS, COR_ASSAYS),
                                 ratio_thresh = DEFAULT_RATIO_THRESH,
                                 threshold_sq = DEFAULT_SPEC_THRESHOLD_SQ,
                                 l2_high_thresh = DEFAULT_L2_HIGH_THRESHOLD,
                                 l2_mid_thresh  = DEFAULT_L2_MID_THRESHOLD,
                                 cov_high_quantile = DEFAULT_COV_HIGH_QUANTILE,
                                 cov_mid_quantile  = DEFAULT_COV_MID_QUANTILE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list()

  # L2 assays: dominance + specificity_summary + specificity_counts_l2
  for (a in intersect(l2_assays, assayNames(se))) {
    mat <- assay(se, a)
    if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
    dom  <- compute_dominance(mat, ratio_thresh)
    spec <- compute_specificity_summary(mat, ratio_thresh = ratio_thresh, threshold_sq = threshold_sq)
    cnt  <- specificity_counts_l2(mat, l2_high_thresh = l2_high_thresh, l2_mid_thresh = l2_mid_thresh)
    saveRDS(dom,  file.path(out_dir, paste0("dominance_", a, ".rds")))
    saveRDS(spec, file.path(out_dir, paste0("specificity_summary_", a, ".rds")))
    saveRDS(cnt,  file.path(out_dir, paste0("specificity_counts_", a, ".rds")))
    results[[a]] <- list(dominance = dom, specificity_summary = spec, specificity_counts = cnt)
    message("  Saved dominance + specificity for ", a)
  }

  # Cov/Cor assays: dominance + specificity_counts_quantile
  for (a in intersect(cov_assays, assayNames(se))) {
    mat <- assay(se, a)
    if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
    dom  <- compute_dominance(mat, ratio_thresh)
    cnt  <- specificity_counts_quantile(mat, cov_high_quantile = cov_high_quantile, cov_mid_quantile = cov_mid_quantile)
    saveRDS(dom, file.path(out_dir, paste0("dominance_", a, ".rds")))
    saveRDS(cnt, file.path(out_dir, paste0("specificity_counts_", a, ".rds")))
    results[[a]] <- list(dominance = dom, specificity_counts = cnt)
    message("  Saved dominance + specificity for ", a)
  }

  invisible(results)
}

# Step 6 — Dominance / specificity on full processed SEs
#   results/{tissue}/full/dominance_{assay}.rds, specificity_*_{assay}.rds
run_full_summaries <- function(se_proc_list, results_dir = "results", ...) {
  message("\n====== Step 6: Full-matrix summaries ======")
  if (is.null(names(se_proc_list))) names(se_proc_list) <- paste0("se", seq_along(se_proc_list))

  out <- lapply(names(se_proc_list), function(tissue) {
    message("\n--- ", tissue, " (full) ---")
    od <- file.path(results_dir, tissue, "full")
    compute_se_summaries(se_proc_list[[tissue]], out_dir = od, ...)
  })
  names(out) <- names(se_proc_list)
  invisible(out)
}

# Step 7 — Dominance / specificity on GWAS-overlapped SEs
#   results/{tissue}/{trait}/dominance_{assay}.rds, specificity_*_{assay}.rds
run_gwas_summaries <- function(overlap_results, results_dir = "results", ...) {
  message("\n====== Step 7: GWAS-overlap summaries ======")

  out <- lapply(names(overlap_results), function(tissue) {
    inner <- overlap_results[[tissue]]
    trait_out <- lapply(names(inner), function(trait) {
      ov <- inner[[trait]]
      if (nrow(ov$se_overlap) == 0) {
        message("\n  Skipping ", tissue, " x ", trait, " (0 overlaps)")
        return(NULL)
      }
      message("\n--- ", tissue, " x ", trait, " ---")
      od <- file.path(results_dir, tissue, trait)
      compute_se_summaries(ov$se_overlap, out_dir = od, ...)
    })
    names(trait_out) <- names(inner)
    trait_out
  })
  names(out) <- names(overlap_results)
  invisible(out)
}

# Dominance plot comparing Full vs GWAS-overlap matrices
plot_dominance_full_vs_gwas <- function(full_mat, gwas_mat, ratio_thresh = DEFAULT_RATIO_THRESH, matrix_label = "L2mat") {
  dom_full <- compute_dominance(full_mat, ratio_thresh)
  dom_gwas <- compute_dominance(gwas_mat, ratio_thresh)
  counts_full <- dom_full %>% dplyr::filter(is_dominant) %>% dplyr::count(maximum_ct)
  counts_gwas <- dom_gwas %>% dplyr::filter(is_dominant) %>% dplyr::count(maximum_ct)
  counts <- dplyr::left_join(counts_full, counts_gwas, by = "maximum_ct")
  colnames(counts) <- c("cell","Full","GWAS")
  counts <- tidyr::gather(counts, set, n, -cell) %>%
    dplyr::group_by(set) %>%
    dplyr::mutate(prop = n / sum(n, na.rm = TRUE))
  ggplot(counts, aes(x = reorder(cell, ifelse(set == "GWAS", prop, NA), na.rm = TRUE), y = prop, fill = set)) +
    geom_bar(stat = "identity", position = 'dodge') +
    theme_bw() + scale_fill_manual(values = c("GWAS" = "red", "Full" = "gray")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cell type", y = "Proportion dominant peaks (ratio > 2)",
         title = paste0("Cell-Type Dominant Peaks (", matrix_label, ")")) +
    coord_flip()
}

# ============================================================
# Pipeline steps — run one at a time
# ============================================================
#
# Step 1: load_se_inputs()          — read SE RDS files into a named list
# Step 2: process_se_list_metrics() — compute all metrics, save & plot
# Step 3: load_trait_inputs()       — read trait GR RDS files into a named list
# Step 4: overlap_all_se_traits()   — overlap every processed SE x every trait, save & plot
# Step 5: summarize_inputs()        — build summary tables of what we have
# Step 6: run_full_summaries()      — dominance/specificity on full processed SEs
# Step 7: run_gwas_summaries()      — dominance/specificity on GWAS-overlapped SEs
#
# Directory structure:
#   data/processed_se/          — processed full SEs (from Step 2)
#   data/GWAS_overlapped_se/    — GWAS-overlapped SEs (from Step 4)
#   plots/{tissue}/             — heatmaps + density plots (from Steps 2 & 4)
#   results/{tissue}/full/      — dominance + specificity results (from Step 6)
#   results/{tissue}/{trait}/   — dominance + specificity results (from Step 7)
#
# Typical usage:
#   se_list      <- load_se_inputs()
#   se_proc_list <- process_se_list_metrics(se_list)
#   trait_list   <- load_trait_inputs()
#   overlaps     <- overlap_all_se_traits(se_proc_list, trait_list)
#   summary      <- summarize_inputs(se_proc_list, trait_list, overlaps)
#   full_res     <- run_full_summaries(se_proc_list)
#   gwas_res     <- run_gwas_summaries(overlaps)

# Step 1 — Load SEs from registry
load_se_inputs <- function(se_inputs = SE_INPUTS) {
  message("\n====== Step 1: Loading SEs ======")
  se_list <- lapply(names(se_inputs), function(nm) {
    message("  Loading: ", nm, " <- ", se_inputs[[nm]])
    readRDS(se_inputs[[nm]])
  })
  names(se_list) <- names(se_inputs)
  message("  Loaded ", length(se_list), " SEs: ", paste(names(se_list), collapse = ", "))
  se_list
}

# Step 2 is process_se_list_metrics() defined above
# — computes metrics, saves to data/processed_se/, plots heatmaps + densities

# Step 3 — Load trait GRs from registry
load_trait_inputs <- function(trait_inputs = TRAIT_INPUTS) {
  message("\n====== Step 3: Loading traits ======")
  trait_list <- lapply(names(trait_inputs), function(nm) {
    message("  Loading: ", nm, " <- ", trait_inputs[[nm]])
    readRDS(trait_inputs[[nm]])
  })
  names(trait_list) <- names(trait_inputs)
  message("  Loaded ", length(trait_list), " traits: ", paste(names(trait_list), collapse = ", "))
  trait_list
}

# Step 4 — Overlap every processed SE x every trait (saves + plots)
overlap_all_se_traits <- function(se_proc_list, trait_list,
                                  overlap_assay = "qn_l2",
                                  plot_dir = "plots") {
  message("\n====== Step 4: Overlapping SE x trait ======")
  if (is.null(names(se_proc_list))) names(se_proc_list) <- paste0("se", seq_along(se_proc_list))
  if (is.null(names(trait_list)))   names(trait_list)   <- paste0("trait", seq_along(trait_list))

  overlap_results <- lapply(names(se_proc_list), function(tissue) {
    se <- se_proc_list[[tissue]]
    inner <- lapply(names(trait_list), function(trait) {
      message("\n  --- ", tissue, " x ", trait, " ---")
      ov <- overlap_se_with_snps(se, trait_list[[trait]],
                                 assay_name = overlap_assay,
                                 tissue = tissue, trait = trait)
      if (nrow(ov$se_overlap) > 0) {
        ov_label <- paste0(tissue, "_", trait)
        plot_assay_cor_heatmaps(ov$se_overlap, ov_label, plot_dir = plot_dir)
        plot_assay_densities(ov$se_overlap, ov_label, plot_dir = plot_dir)
      }
      ov
    })
    names(inner) <- names(trait_list)
    inner
  })
  names(overlap_results) <- names(se_proc_list)
  overlap_results
}

# Step 5 — Summarise what we have
summarize_inputs <- function(se_proc_list, trait_list = NULL, overlap_results = NULL) {
  message("\n====== Step 5: Input summary ======")

  se_summary <- data.frame(
    tissue      = names(se_proc_list),
    n_celltypes = sapply(se_proc_list, ncol),
    n_peaks     = sapply(se_proc_list, nrow),
    stringsAsFactors = FALSE
  )

  if (!is.null(trait_list)) {
    trait_summary <- data.frame(
      trait  = names(trait_list),
      n_snps = sapply(trait_list, length),
      stringsAsFactors = FALSE
    )
  } else {
    trait_summary <- NULL
  }

  if (!is.null(overlap_results)) {
    ov_rows <- do.call(rbind, lapply(names(overlap_results), function(tissue) {
      inner <- overlap_results[[tissue]]
      do.call(rbind, lapply(names(inner), function(trait) {
        data.frame(tissue = tissue, trait = trait,
                   n_overlapped = nrow(inner[[trait]]$se_overlap),
                   stringsAsFactors = FALSE)
      }))
    }))
  } else {
    ov_rows <- NULL
  }

  message("\n-- SE datasets --")
  print(se_summary)
  if (!is.null(trait_summary)) { message("\n-- Traits --"); print(trait_summary) }
  if (!is.null(ov_rows))       { message("\n-- Overlaps --"); print(ov_rows) }

  invisible(list(se_summary = se_summary, trait_summary = trait_summary,
                 overlap_summary = ov_rows))
}





