#!/usr/bin/env Rscript
# ============================================================

# 26_02_2026 Fixed stackbar plot, save processed SE under data
# ============================================================
# run_clear_pair.R    ONE snATAC SE × GWAS trait pair. HPC accessible.

# incorporate covW ideas, specificity summary (dominance ratio idea) ; not including the per signal coherence in this version 
# 
# Output lands in the current working directory under:
#   results/                     tsv
#   plots/                       PDFs
#   data/                      processed SE / GWAS SE
# ============================================================

# Load
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
  library(preprocessCore)
  library(ggsci)
})

# arguments 
args <- commandArgs(TRUE)

se_name    <- args[1] # List of tissue + path
se_path    <- args[2]
trait_name <- args[3] # List of trait + path    
trait_path <- args[4]

message("\n========================================")
message("CLEAR analysis: ", se_name, " x ", trait_name)
message("  SE path:    ", se_path)
message("  Trait path: ", trait_path)
message("========================================\n")

# Output directories bash cd ${WD}
results_dir    <- "results"
plots_dir      <- "plots"
data_dir       <- "data"           # store processed SE per pair (but would be the same FULL SE each tissue jsut different GWAS-overlap subset SE)
# create qn directories (could be used as benchmark or tests done during development, supplementary decision making)
qn_results_dir <- file.path(results_dir, "qn")
qn_plots_dir   <- file.path(plots_dir, "qn")
dir.create(results_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(qn_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qn_plots_dir,   recursive = TRUE, showWarnings = FALSE)

# ===========================================================
# Helper functions (from L2_and_covW.R)
# ===========================================================

# GRanges from rowData columns (seqnames, start, end) to make RSE (some aren't)
# Also ensures rownames and peakid are set
prepare_clear_se <- function(se, genome = "hg38") {
  rd <- as.data.frame(rowData(se))
  
  # Build GRanges from rowData (ArchR groupSE always has seqnames/start/end)
  if (!all(c("seqnames", "start", "end") %in% names(rd))) {
    stop("rowData must contain seqnames, start, end columns (ArchR PeakMatrix format expected)") # Check what Signac does
  }
  
  gr <- GRanges(
    seqnames = rd$seqnames,
    ranges   = IRanges::IRanges(start = as.integer(rd$start), end = as.integer(rd$end))
  )
  
  # Keep rowData columns 
  extra_cols <- setdiff(names(rd), c("seqnames", "start", "end", "width", "strand"))
  if (length(extra_cols)) mcols(gr) <- rd[, extra_cols, drop = FALSE]
  
  # Rebuild as RSE
  se <- SummarizedExperiment(
    assays    = as.list(assays(se)),
    rowRanges = gr,
    colData   = colData(se)
  )
  
  # Rownames as genomic coords; peakid in rowData
  rownames(se) <- paste0(as.character(seqnames(gr)), ":", start(gr), "-", end(gr))
  rowData(se)$peakid <- rownames(se)
  
  try(seqlevelsStyle(rowRanges(se)) <- "UCSC", silent = TRUE)
  try(GenomeInfoDb::genome(rowRanges(se)) <- genome, silent = TRUE)
  
  se
}

safe_l2 <- function(x) {
  s <- sum(x^2)
  if (s == 0) x else x / sqrt(s)
}

add_rank_desc <- function(M) apply(M, 2, function(x) rank(-x, ties.method = "average"))

L2_mat <- function(mat) {
  row_norms <- sqrt(rowSums(mat^2))
  row_norms[row_norms == 0] <- 1
  mat / row_norms
}
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
  
  # Worth adding a l2_squared_mat (stacked bar and dominance)
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

PLOT_ASSAYS <- c("raw", "cov_raw_l2", "raw_l2", "cov_raw",
                 "qn", "qn_l2", "cov_qn", "cov_qn_l2",
                 "cor_raw", "cor_qn", "cor_raw_l2", "cor_qn_l2")

# Correlation heatmaps, heatmap show dendrograms
plot_assay_cor_heatmaps <- function(se, tissue, plot_dir = "plots") {
  out <- file.path(plot_dir, tissue)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  avail <- intersect(PLOT_ASSAYS, assayNames(se))
  for (a in avail) {
    mat <- assay(se, a)
    if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
    if (ncol(mat) < 2) next
    cor_mat <- cor(mat, use = "pairwise.complete.obs")
    fname <- file.path(out, paste0(tissue, "_", a, "_cor_heatmap.pdf"))
    pdf(fname, width = 8, height = 7)
    ht <- ComplexHeatmap::Heatmap(
      cor_mat,
      name                  = "Correlation",
      column_title          = paste0(tissue, " \u2014 ", a, " correlation"),
      clustering_method_rows    = "ward.D2",
      clustering_method_columns = "ward.D2",
      clustering_distance_rows    = "euclidean",
      clustering_distance_columns = "euclidean",
      show_row_dend    = TRUE,
      show_column_dend = TRUE,
      row_names_gp     = grid::gpar(fontsize = 7),
      column_names_gp  = grid::gpar(fontsize = 7),
      col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    )
    ComplexHeatmap::draw(ht)
    dev.off()
    message("  Saved: ", fname)
  }
  invisible(out)
}

# Density plots  one PDF with faceted panels per assay
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
  # Choose ggsci colour palette based on number of cell types
  n_ct_dens <- length(unique(long$cell_type))
  if (n_ct_dens <= 10) {
    colour_scale <- ggsci::scale_color_npg()
  } else if (n_ct_dens <= 20) {
    colour_scale <- ggsci::scale_color_d3(palette = "category20")
  } else {
    colour_scale <- ggsci::scale_color_igv()
  }
  p <- ggplot(long, aes(x = value, colour = cell_type)) +
    geom_density() +
    facet_wrap(~ assay, scales = "free", ncol = 4) +
    colour_scale +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 9)) +
    labs(title = paste0(tissue, " \u2014 per-cell-type distributions"),
         x = "Value", y = "Density", colour = "Cell type")
  fname <- file.path(out, paste0(tissue, "_assay_densities.pdf"))
  ggsave(fname, p, width = 16, height = max(12, ceiling(length(avail) / 4) * 3))
  message("  Saved: ", fname)
  invisible(fname)
}

# Per peak dominance
# dominated as 2 times of the next in the specificity matrix (row-wise)
compute_dominance <- function(mat, ratio_thresh = 2) { 
  message("  identifying dominant cell types...")
  
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  if (n_rows == 0L) {
    return(data.frame(maximum_ct = character(), ratio = numeric(), is_dominant = logical(),
                      peak = character(), stringsAsFactors = FALSE))
  }
  
  if (is.null(colnames(mat))) colnames(mat) <- paste0("ct", seq_len(n_cols))
  
  if (n_cols < 2L) {
    warning("compute_dominance: fewer than 2 cell types; ratios set to NA")
    return(data.frame(
      maximum_ct  = rep(colnames(mat)[1], n_rows),
      ratio       = rep(NA_real_, n_rows),
      is_dominant = rep(FALSE, n_rows),
      peak        = rownames(mat),
      stringsAsFactors = FALSE
    ))
  }
  
  cell_names <- colnames(mat)
  
  # Vectorised: find max and second-max per row (avoids slow t(apply()) and
  # the name-mangling bug where c(ratio = named_value) → "ratio.CellName")
  max_idx  <- max.col(mat, ties.method = "first")
  max_val  <- mat[cbind(seq_len(n_rows), max_idx)]
  
  mat2 <- mat
  mat2[cbind(seq_len(n_rows), max_idx)] <- -Inf
  second_val <- mat2[cbind(seq_len(n_rows), max.col(mat2, ties.method = "first"))]
  
  ratio <- ifelse(second_val == 0, Inf, max_val / second_val)
  
  data.frame(
    maximum_ct  = cell_names[max_idx],
    ratio       = ratio,
    is_dominant  = ratio > ratio_thresh,
    peak        = rownames(mat),
    stringsAsFactors = FALSE
  )
}

# Adapted from clear_test
# Selects peaks where the top k (≤ max_k) cell types cumulatively
# explain ≥ cumvar_threshold of total variance (L2² sums to 1),
# and the smallest of those k cell types has L2 ≥ 1/sqrt(N/2)
plot_celltype_stackbars <- function(
    se_overlap,
    value_metric       = "raw_l2",
    cumvar_threshold   = 0.8, # cumulative cutoff, e.g. 3 cell types can have dominant if 40/20/20 or 60/10/10 (how does this show in tracks?) or 30/30/30 could also be interesting to check lineages (composition)
    max_k              = 3, # Consider using half of the number of cell types, N_lineages could be better?
    outfile,
    max_peaks = 200
) {
  stopifnot(!missing(outfile))
  
  if (!value_metric %in% assayNames(se_overlap)) {
    warning("Assay '", value_metric, "' not found; skipping stack bar plots.")
    return(invisible(FALSE))
  }
  
  val_mat <- SummarizedExperiment::assay(se_overlap, value_metric)
  if (inherits(val_mat, "Matrix")) val_mat <- as.matrix(val_mat)
  storage.mode(val_mat) <- "double"
  
  # My thoughts for specificity cutoff: 
  n_celltypes <- ncol(val_mat)
  min_l2 <- 1 / sqrt(n_celltypes / 2)
  val_sq <- val_mat^2
  
  peak_info <- lapply(seq_len(nrow(val_mat)), function(i) {
    row_sq  <- val_sq[i, ]
    row_l2  <- val_mat[i, ]
    ord     <- order(row_sq, decreasing = TRUE)
    cum_sq  <- cumsum(row_sq[ord])
    
    k <- which(cum_sq >= cumvar_threshold)[1]
    if (is.na(k) || k > max_k) return(NULL)
    
    top_k_idx <- ord[seq_len(k)]
    if (row_l2[top_k_idx[k]] < min_l2) return(NULL)
    
    list(
      idx       = i,
      k         = k,
      cum_var   = cum_sq[k],
      top_cells = colnames(val_mat)[top_k_idx],
      min_l2    = row_l2[top_k_idx[k]]
    )
  })
  peak_info <- Filter(Negate(is.null), peak_info)
  
  if (!length(peak_info)) {
    warning("No peaks meet criteria (top \u2264", max_k, " cells cumulative L2\u00B2 \u2265 ",
            cumvar_threshold, ", min L2 \u2265 ", round(min_l2, 3), ")")
    grDevices::pdf(outfile, width = 10, height = 6)
    plot.new(); text(0.5, 0.5, "No peaks meeting criteria")
    dev.off()
    return(invisible(FALSE))
  }
  
  sel <- sapply(peak_info, `[[`, "idx")
  sel_k <- sapply(peak_info, `[[`, "k")
  
  if (length(sel) > max_peaks) {
    ord <- order(sel_k, -apply(val_sq[sel, , drop = FALSE], 1, max))
    keep <- ord[seq_len(max_peaks)]
    sel <- sel[keep]
    sel_k <- sel_k[keep]
    peak_info <- peak_info[keep]
  }
  
  # Labels
  rd <- as.data.frame(rowData(se_overlap))
  peakid     <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se_overlap)
  ccv_label  <- if ("CCVs"  %in% names(rd)) rd$CCVs  else ""
  sig_label  <- if ("signal" %in% names(rd)) rd$signal else ""
  base_label <- ifelse(nzchar(ccv_label), ccv_label, peakid)
  labels     <- ifelse(nzchar(sig_label), paste0(base_label, " (", sig_label, ")"), base_label)
  
  submat <- val_sq[sel, , drop = FALSE]
  
  # Identify the dominant cell type for each peak (for ordering)
  dom_ct <- apply(val_mat[sel, , drop = FALSE], 1, function(x) names(which.max(x)))
  
  # ---- Build long-format with per-peak ordering ----
  # For each peak, rank cell types by L2² (descending) and assign a
  # stacking order so the biggest contributor sits at the bottom (y = 0)
  long_list <- lapply(seq_along(sel), function(j) {
    i <- sel[j]
    row_sq <- val_sq[i, ]
    ct_ord <- order(row_sq, decreasing = TRUE)   # biggest first
    data.frame(
      peak_label = labels[i],
      celltype   = colnames(val_sq)[ct_ord],
      value      = row_sq[ct_ord],
      stack_rank = seq_along(ct_ord),             # 1 = biggest → plotted first (bottom)
      k          = sel_k[j],
      dom_ct     = dom_ct[j],
      max_sq     = max(row_sq),
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, long_list)
  
  # Peak y-axis ordering: by k, then dominant cell type, then max score
  peak_order <- long %>%
    dplyr::distinct(peak_label, k, dom_ct, max_sq) %>%
    dplyr::arrange(k, dom_ct, -max_sq) %>%
    dplyr::pull(peak_label)
  long$peak_label <- factor(long$peak_label, levels = rev(peak_order))
  
  # Cell type fill ordering: use stack_rank so within each bar the biggest
  
  # value is placed first (bottom).  We achieve this by converting celltype
  # to a factor whose levels are ordered by *overall* total (for the legend)
  # and then using geom_bar's group aesthetic driven by stack_rank.
  ct_totals <- tapply(long$value, long$celltype, sum, na.rm = TRUE)
  ct_levels <- names(sort(ct_totals, decreasing = TRUE))
  long$celltype <- factor(long$celltype, levels = ct_levels)
  
  # To force per-peak stacking order (biggest at bottom), we use
  # fct_reorder inside each peak via geom_col + position_stack(reverse=FALSE)
  # with the data pre-sorted so biggest comes first per peak.
  long <- long %>% dplyr::arrange(peak_label, -value)
  
  threshold_pct <- cumvar_threshold * 100
  n_sel <- length(sel)
  n_ct <- length(ct_levels)
  
  # Choose ggsci palette  d3 has 20 colours, igv has 51
  if (n_ct <= 10) {
    fill_scale <- ggsci::scale_fill_npg()
  } else if (n_ct <= 20) {
    fill_scale <- ggsci::scale_fill_d3(palette = "category20")
  } else {
    fill_scale <- ggsci::scale_fill_igv()
  }
  
  p <- ggplot(long, aes(x = peak_label, y = value, fill = celltype,
                        group = factor(stack_rank))) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.85) +
    geom_hline(yintercept = cumvar_threshold, linetype = "dashed",
               color = "black", linewidth = 0.3) +
    coord_flip() +
    fill_scale +
    labs(
      title = paste0("GWAS peaks: top \u2264", max_k, " cell types explain \u2265",
                     threshold_pct, "% variance",
                     " (min L2 \u2265 ", round(min_l2, 3), ")"),
      subtitle = paste0(n_sel, " peaks selected"),
      x = "Peak / CCV",
      y = "Variance explained (L2\u00B2, cumulative = 1.0)",
      fill = "Cell type"
    ) +
    theme_bw() +
    theme(
      axis.text.y  = element_text(size = 4),
      axis.text.x  = element_text(size = 6),
      legend.text  = element_text(size = 6),
      legend.title = element_text(size = 8),
      plot.title   = element_text(size = 10)
    )
  
  pdf_height <- max(6, n_sel * 0.12 + 3)
  grDevices::pdf(outfile, width = 12, height = min(pdf_height, 50))
  print(p)
  grDevices::dev.off()
  
  message("  \u2713 Saved stackbar plot (", n_sel, " peaks): ", outfile)
  invisible(TRUE)
}

# ===========================================================
# Constants for specificity
# ===========================================================
DEFAULT_RATIO_THRESH      <- 2
DEFAULT_SPEC_THRESHOLD_SQ <- 0.8 # or 0.5?

# ===========================================================
# compute_specificity_summary
# ===========================================================
compute_specificity_summary <- function(mat_l2, cell_cor = NULL,
                                        ratio_thresh  = DEFAULT_RATIO_THRESH,
                                        threshold_sq  = DEFAULT_SPEC_THRESHOLD_SQ) {
  if (is.null(cell_cor)) cell_cor <- cor(mat_l2, use = "pairwise.complete.obs") # specify the cor_raw or cor_raw_l2 ???
  results <- apply(mat_l2, 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    x_sorted     <- x[ord] 
    cells_sorted <- colnames(mat_l2)[ord] # cell names sorted by specificity score for this peak
    cum_sq <- cumsum(x_sorted^2) # L2^2
    k <- which(cum_sq >= threshold_sq)[1] # Number of top cell types needed to reach threshold
    if (is.na(k)) k <- length(x_sorted)
    top_cells <- cells_sorted[seq_len(k)] # 
    maximum_cell   <- cells_sorted[1]
    maximum_score  <- x_sorted[1]
    dominant_ratio <- if (length(x_sorted) > 1) x_sorted[1] / x_sorted[2] else NA
    coherence <- if (length(top_cells) > 1) {
      idx <- match(top_cells, colnames(mat_l2))
      sub_cor <- cell_cor[idx, idx, drop = FALSE]
      mean(sub_cor[upper.tri(sub_cor)])
    } else {
      NA
    }
    c(
      maximum_cell   = maximum_cell,
      maximum_score  = as.numeric(maximum_score),
      dominant_ratio = as.numeric(dominant_ratio),
      is_dominant    = as.character(dominant_ratio > ratio_thresh),
      k_cells        = as.numeric(k),
      top_k_cells    = paste(top_cells, collapse = ","),
      coherence      = as.numeric(coherence)
    )
  })
  out <- as.data.frame(t(results), stringsAsFactors = FALSE)
  out$peak <- rownames(mat_l2)
  # Coerce numeric columns (apply returns character matrix)
  out$maximum_score  <- as.numeric(out$maximum_score)
  out$dominant_ratio <- as.numeric(out$dominant_ratio)
  out$is_dominant    <- as.logical(out$is_dominant)
  out$k_cells        <- as.integer(out$k_cells)
  out$coherence      <- as.numeric(out$coherence)
  rownames(out) <- NULL
  out
}

# ===========================================================
#   Compare genome-wide vs GWAS dominant cell-type proportions (dodged plot)
# ===========================================================
dominance_dodge_plot <- function(genome_mat, gwas_mat,
                                 ratio_thresh = DEFAULT_RATIO_THRESH,
                                 title_label  = "") {
  dom_genome <- compute_dominance(genome_mat, ratio_thresh = ratio_thresh)
  dom_gwas   <- compute_dominance(gwas_mat,   ratio_thresh = ratio_thresh)

  n_dom_genome <- sum(dom_genome$is_dominant)
  n_dom_gwas   <- sum(dom_gwas$is_dominant)

  # Proportion within dominant peaks per set
  df <- dplyr::bind_rows(
    dom_genome %>% dplyr::filter(is_dominant) %>% dplyr::count(maximum_ct) %>%
      dplyr::mutate(proportion = n / n_dom_genome, set = "Genome-wide"),
    dom_gwas   %>% dplyr::filter(is_dominant) %>% dplyr::count(maximum_ct) %>%
      dplyr::mutate(proportion = n / n_dom_gwas,   set = "GWAS")
  )

  # Fill missing cell type × set combos with 0
  df <- tidyr::complete(df, maximum_ct, set,
                        fill = list(n = 0L, proportion = 0))

  df$set       <- factor(df$set, levels = c("Genome-wide", "GWAS"))
  df$bar_label <- ifelse(df$n > 0,
                         paste0(df$n, " (", round(df$proportion * 100, 1), "%)"), "")

  # Order cell types by GWAS proportion
  ct_order <- df %>% dplyr::filter(set == "GWAS") %>%
    dplyr::arrange(proportion) %>% dplyr::pull(maximum_ct)
  df$maximum_ct <- factor(df$maximum_ct, levels = ct_order)

  ggplot(df, aes(x = maximum_ct, y = proportion, fill = set)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = bar_label),
              position = position_dodge(width = 0.8), hjust = -0.05, size = 2.8) +
    coord_flip() +
    scale_fill_manual(values = c("Genome-wide" = "grey60", "GWAS" = "firebrick")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    theme_bw() +
    labs(
      title    = paste0("Dominant cell types (ratio > ", ratio_thresh, ") ", title_label),
      subtitle = paste0("Genome-wide: ", nrow(genome_mat), " peaks (", n_dom_genome,
                        " dominant);  GWAS: ", nrow(gwas_mat), " peaks (", n_dom_gwas, " dominant)"),
      x    = "Cell type",
      y    = "Proportion within dominant peaks",
      fill = "Peak set"
    )
}

# ===========================================================
# CLEAR full run through 
message("\n--- Running CLEAR combination stated from bash script  ---")
# ===========================================================

# ===========================================================
# Step 1: Load inputs
# ===========================================================
message("\n--- Loading inputs ---")
se   <- readRDS(se_path)
snps <- readRDS(trait_path)

# Ensure SE has proper rowRanges 
se <- prepare_clear_se(se)

message("  SE:    ", nrow(se), " peaks x ", ncol(se), " cell types")
message("  SNPs:  ", length(snps))


# ===========================================================
# Step 2: Compute all metrics (L2, QN, cov/cor weighted)
# ===========================================================
message("\n--- Computing all metrics (L2, QN, cov/cor weighted) ---")
se <- compute_metrics_se(se)
mat_l2    <- assay(se, "raw_l2")
mat_qn_l2 <- assay(se, "qn_l2")

# Save processed full SE but the same for each tissue across different traits
se_save_name <- paste0(se_name, "_", ncol(se), "celltypes_", nrow(se), "peaks.rds")
saveRDS(se, file.path(data_dir, se_save_name))
message("  \u2713 Saved processed SE: ", file.path(data_dir, se_save_name))

# ===========================================================
# Step 3: Intersect GWAS SNPs with L2 matrix
# ===========================================================
message("\n--- GWAS SNP intersection ---")
hits_se  <- findOverlaps(snps, rowRanges(se))
gwas_idx <- unique(subjectHits(hits_se))

if (!length(gwas_idx)) {
  message("  No GWAS overlaps; downstream GWAS-specific outputs will be empty.")
  se_gwas    <- se[0, , drop = FALSE]
  gwas_mat   <- mat_l2[0, , drop = FALSE]
  gwas_qn_l2 <- mat_qn_l2[0, , drop = FALSE]
} else {
  se_gwas <- se[gwas_idx, , drop = FALSE]
  
  ann <- data.frame( # annotated GWAS overlaps for locus-level summaries
    peak_idx = subjectHits(hits_se),
    CCVs     = mcols(snps)$names[queryHits(hits_se)],
    signal   = mcols(snps)$signal[queryHits(hits_se)],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(peak_idx) %>%
    dplyr::summarise(
      CCVs   = paste(unique(na.omit(CCVs)), collapse = ","),
      signal = paste(unique(na.omit(signal)), collapse = ","),
      .groups = "drop"
    )
  
  rd <- as.data.frame(rowData(se_gwas))
  if (!"CCVs" %in% names(rd))   rd$CCVs   <- ""
  if (!"signal" %in% names(rd)) rd$signal <- ""
  pos <- match(ann$peak_idx, gwas_idx)
  rd$CCVs[pos]   <- ann$CCVs
  rd$signal[pos] <- ann$signal
  rowData(se_gwas) <- S4Vectors::DataFrame(rd, row.names = rownames(se_gwas))
  
  gwas_mat   <- assay(se_gwas, "raw_l2")
  gwas_qn_l2 <- assay(se_gwas, "qn_l2")
}
message("  GWAS peaks: ", nrow(gwas_mat))

# Save GWAS-overlapped SE
if (nrow(gwas_mat) > 0) {
  n_sig_overlap <- length(unique(na.omit(unlist(strsplit(
    as.character(rowData(se_gwas)$signal), ",")))))
  gwas_se_name <- paste0(se_name, "_", trait_name, "_",
                         nrow(gwas_mat), "peaks_", n_sig_overlap, "signals.rds")
  saveRDS(se_gwas, file.path(results_dir, gwas_se_name))
  message("  \u2713 Saved GWAS SE: ", file.path(results_dir, gwas_se_name))
}

# ===========================================================
# Step 4: Correlation heatmaps (panel across all assays)
# ===========================================================
message("\n--- Correlation heatmaps (panel) ---")

plot_assay_cor_heatmaps(se, se_name, plots_dir)
plot_assay_densities(se, se_name, plots_dir)
if (nrow(gwas_mat) > 1) {
  plot_assay_cor_heatmaps(se_gwas, paste0(se_name, "_", trait_name, "_GWAS"), plots_dir)
  plot_assay_densities(se_gwas, paste0(se_name, "_", trait_name, "_GWAS"), plots_dir)
}

# ===========================================================
# Step 5: Input summary
# ===========================================================
n_unique_signals <- length(unique(na.omit(mcols(snps)$signal)))

# Calculate number of unique signals overlapping peaks
if (length(gwas_idx) > 0) {
  overlapping_signals <- mcols(snps)$signal[queryHits(hits_se)]
  n_unique_signals_in_peaks <- length(unique(na.omit(overlapping_signals)))
} else {
  n_unique_signals_in_peaks <- 0
}

# Raw PeakMatrix column sums: check comparability across cell types
raw_mat <- assay(se, "raw")
raw_col_sums <- colSums(raw_mat)
raw_colsum_summary <- data.frame(
  celltype         = names(raw_col_sums),
  total_raw_signal = as.numeric(raw_col_sums),
  stringsAsFactors = FALSE
)
raw_colsum_summary <- raw_colsum_summary[order(-raw_colsum_summary$total_raw_signal), ]
write.table(raw_colsum_summary, file.path(results_dir, "raw_signal_per_celltype.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("\n--- Raw PeakMatrix column sums ---")
print(raw_colsum_summary)

input_summary <- data.frame(
  se_name                    = se_name,
  trait_name                 = trait_name,
  GWAS_trait                 = basename(trait_path),
  snATAC_dataset             = basename(se_path),
  n_peaks                    = nrow(mat_l2),
  n_celltypes                = ncol(mat_l2),
  n_snps                     = length(snps),
  n_unique_signals           = n_unique_signals,
  n_gwas_peaks               = nrow(gwas_mat),
  n_unique_signals_in_peaks  = n_unique_signals_in_peaks,
  mean_total_raw_signal      = mean(raw_col_sums),
  var_total_raw_signal       = var(as.numeric(raw_col_sums)),
  max_signal_celltype        = names(which.max(raw_col_sums)),
  max_signal_value           = max(raw_col_sums),
  min_signal_celltype        = names(which.min(raw_col_sums)),
  min_signal_value           = min(raw_col_sums),
  stringsAsFactors = FALSE
)
write.table(input_summary, file.path(results_dir, "input_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("\n--- Input summary ---")
print(input_summary)

# ===========================================================
# Step 6: L2 specificity (GWAS peaks) - COMBINED
# ===========================================================
message("\n--- L2 specificity summary (combined) ---")
high_thresh <- 1 / sqrt(3)  # 3 or less cell types matter. Threshold breaks when correlated CTs more than 3 and split signal 
mid_thresh  <- 1 / sqrt(ncol(mat_l2)) # Uniform null
if (nrow(gwas_mat) > 0) {
  
  # Compute detailed metrics first
  specificity_summary_l2_detailed <- compute_specificity_summary(gwas_mat)
  
  # Compute basic metrics
  specificity_basic <- data.frame(
    peak   = rownames(gwas_mat),
    n_high = rowSums(gwas_mat > high_thresh),
    n_mid  = rowSums(gwas_mat > mid_thresh & gwas_mat <= high_thresh),
    n_low  = rowSums(gwas_mat <= mid_thresh)
  )
  high_cells <- apply(gwas_mat, 1, function(x) names(x)[x > high_thresh])
  specificity_basic$high_cells    <- sapply(high_cells, paste, collapse = ", ")
  specificity_basic$multiple_high <- lengths(high_cells) > 1
  
  # Merge: detailed + basic + GWAS annotation
  specificity_summary_l2 <- specificity_summary_l2_detailed %>%
    left_join(specificity_basic, by = "peak")
  
  # Add GWAS SNP annotation
  peak_gr <- GRanges(specificity_summary_l2$peak)
  hits <- findOverlaps(peak_gr, snps)
  snp_list  <- tapply(mcols(snps)$names[subjectHits(hits)],  queryHits(hits), paste, collapse = ",")
  locus_vec <- tapply(mcols(snps)$signal[subjectHits(hits)], queryHits(hits), function(x) paste(unique(x), collapse = ", "))
  specificity_summary_l2$gwas_snps  <- NA_character_
  specificity_summary_l2$gwas_locus <- NA_character_
  specificity_summary_l2$gwas_snps[as.numeric(names(snp_list))]   <- snp_list
  specificity_summary_l2$gwas_locus[as.numeric(names(locus_vec))] <- locus_vec
  
  # Reorder columns for clarity
  specificity_summary_l2 <- specificity_summary_l2 %>%
    select(peak, gwas_snps, gwas_locus, 
           maximum_cell, maximum_score, dominant_ratio, is_dominant,
           n_high, n_mid, n_low, high_cells, multiple_high,
           k_cells, top_k_cells, coherence)
  
  write.table(specificity_summary_l2, file.path(results_dir, "specificity_summary_l2.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Specificity bar: order 1, 2, ... first, 0 at end; add count labels
  count_df <- as.data.frame(table(specificity_summary_l2$n_high), stringsAsFactors = FALSE)
  colnames(count_df) <- c("n_high", "count")
  count_df$proportion <- count_df$count / sum(count_df$count)
  # Reorder: 1, 2, 3, ... then 0 at the end
  lvls <- as.character(sort(as.integer(count_df$n_high)))
  if ("0" %in% lvls) lvls <- c(setdiff(lvls, "0"), "0")
  count_df$n_high <- factor(count_df$n_high, levels = lvls)
  count_df$bar_label <- paste0(count_df$count, " (", round(count_df$proportion * 100, 1), "%)")
  p_spec_l2 <- ggplot(count_df, aes(x = "GWAS peaks", y = proportion, fill = n_high)) +
    geom_bar(stat = "identity", width = 0.6, position = "stack") +
    geom_text(aes(label = bar_label), position = position_stack(vjust = 0.5), size = 3) +
    theme_bw() +
    labs(title = paste0(se_name, " x ", trait_name, " - High-specific cell counts (L2)"),
         x = "", y = "Proportion", fill = "# Highly specific\ncell types") + coord_flip()
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_l2_specificity_bar.pdf")),
         p_spec_l2, width = 8, height = 4)
  
  # Single stackbar plot: peaks with ≤3 cell types at L2² ≥ 0.1, min L2 ≥ 1/sqrt(N/2)
  plot_celltype_stackbars(se_gwas, "raw_l2", cumvar_threshold = 0.5,
                          outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_cum50.pdf")))
  plot_celltype_stackbars(se_gwas, "raw_l2", cumvar_threshold = 0.8,
                          outfile = file.path(plots_dir, paste0(se_name, "_", trait_name, "_stackbars_l2_cum80.pdf")))
  
} else {
  specificity_summary_l2 <- data.frame()
  message("  Skipping L2 specificity summaries (no GWAS-overlapping peaks)")
}

# ===========================================================
# Step 6b: QN L2 specificity (GWAS peaks) - COMBINED
# ===========================================================
message("\n--- QN L2 specificity summary (combined) ---")
if (nrow(gwas_qn_l2) > 0) {
  
  specificity_summary_qn_l2_detailed <- compute_specificity_summary(gwas_qn_l2)
  
  specificity_basic_qn <- data.frame(
    peak   = rownames(gwas_qn_l2),
    n_high = rowSums(gwas_qn_l2 > high_thresh),
    n_mid  = rowSums(gwas_qn_l2 > mid_thresh & gwas_qn_l2 <= high_thresh),
    n_low  = rowSums(gwas_qn_l2 <= mid_thresh)
  )
  high_cells_qn <- apply(gwas_qn_l2, 1, function(x) names(x)[x > high_thresh])
  specificity_basic_qn$high_cells    <- sapply(high_cells_qn, paste, collapse = ", ")
  specificity_basic_qn$multiple_high <- lengths(high_cells_qn) > 1
  
  specificity_summary_qn_l2 <- specificity_summary_qn_l2_detailed %>%
    left_join(specificity_basic_qn, by = "peak")
  
  peak_gr_qn <- GRanges(specificity_summary_qn_l2$peak)
  hits_qn <- findOverlaps(peak_gr_qn, snps)
  snp_list_qn  <- tapply(mcols(snps)$names[subjectHits(hits_qn)],  queryHits(hits_qn), paste, collapse = ",")
  locus_vec_qn <- tapply(mcols(snps)$signal[subjectHits(hits_qn)], queryHits(hits_qn), function(x) paste(unique(x), collapse = ", "))
  specificity_summary_qn_l2$gwas_snps  <- NA_character_
  specificity_summary_qn_l2$gwas_locus <- NA_character_
  specificity_summary_qn_l2$gwas_snps[as.numeric(names(snp_list_qn))]   <- snp_list_qn
  specificity_summary_qn_l2$gwas_locus[as.numeric(names(locus_vec_qn))] <- locus_vec_qn
  
  specificity_summary_qn_l2 <- specificity_summary_qn_l2 %>%
    select(peak, gwas_snps, gwas_locus,
           maximum_cell, maximum_score, dominant_ratio, is_dominant,
           n_high, n_mid, n_low, high_cells, multiple_high,
           k_cells, top_k_cells, coherence)
  
  write.table(specificity_summary_qn_l2, file.path(qn_results_dir, "specificity_summary_qn_l2.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Specificity bar: order 1, 2, ... first, 0 at end; add count labels
  count_df_qn <- as.data.frame(table(specificity_summary_qn_l2$n_high), stringsAsFactors = FALSE)
  colnames(count_df_qn) <- c("n_high", "count")
  count_df_qn$proportion <- count_df_qn$count / sum(count_df_qn$count)
  lvls_qn <- as.character(sort(as.integer(count_df_qn$n_high)))
  if ("0" %in% lvls_qn) lvls_qn <- c(setdiff(lvls_qn, "0"), "0")
  count_df_qn$n_high <- factor(count_df_qn$n_high, levels = lvls_qn)
  count_df_qn$bar_label <- paste0(count_df_qn$count, " (", round(count_df_qn$proportion * 100, 1), "%)")
  p_spec_qn <- ggplot(count_df_qn, aes(x = "GWAS peaks", y = proportion, fill = n_high)) +
    geom_bar(stat = "identity", width = 0.6, position = "stack") +
    geom_text(aes(label = bar_label), position = position_stack(vjust = 0.5), size = 3) +
    theme_bw() +
    labs(title = paste0(se_name, " x ", trait_name, " - High-specific cell counts (QN-L2)"),
         x = "", y = "Proportion", fill = "# Highly specific\ncell types") + coord_flip()
  ggsave(file.path(qn_plots_dir, paste0(se_name, "_", trait_name, "_qn_l2_specificity_bar.pdf")),
         p_spec_qn, width = 8, height = 4)
  
  # Single stackbar plot for QN
  plot_celltype_stackbars(se_gwas, "qn_l2", cumvar_threshold = 0.5,
                          outfile = file.path(qn_plots_dir, paste0(se_name, "_", trait_name, "_stackbars_qn_l2_cum50.pdf")))
  plot_celltype_stackbars(se_gwas, "qn_l2", cumvar_threshold = 0.8,
                          outfile = file.path(qn_plots_dir, paste0(se_name, "_", trait_name, "_stackbars_qn_l2_cum80.pdf")))
  
} else {
  specificity_summary_qn_l2 <- data.frame()
  message("  Skipping QN L2 specificity summaries (no GWAS-overlapping peaks)")
}

# ===========================================================
# Step 7: L2 dominance  genome-wide vs GWAS
# ===========================================================
message("\n--- L2 dominance: genome-wide vs GWAS ---")
p_dom_l2 <- dominance_dodge_plot(mat_l2, gwas_mat, title_label = "L2mat")
ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_dominance_l2.pdf")),
       p_dom_l2, width = 8, height = 6)

# --- QN L2 dominance ---
message("\n--- QN L2 dominance: genome-wide vs GWAS ---")
p_dom_qn_l2 <- dominance_dodge_plot(mat_qn_l2, gwas_qn_l2, title_label = "QN-L2mat")
ggsave(file.path(qn_plots_dir, paste0(se_name, "_", trait_name, "_dominance_qn_l2.pdf")),
       p_dom_qn_l2, width = 8, height = 6)

# ===========================================================
# Step 8: Covariance-weighted matrices (already computed by compute_metrics_se)
# ===========================================================
message("\n--- Extracting cov-weighted matrices ---")
wmat      <- assay(se, "cov_raw_l2")
wmat_qn   <- assay(se, "cov_qn_l2")
gwas_wmat    <- if (nrow(se_gwas)) assay(se_gwas, "cov_raw_l2") else wmat[0, , drop = FALSE]
gwas_wmat_qn <- if (nrow(se_gwas)) assay(se_gwas, "cov_qn_l2")  else wmat_qn[0, , drop = FALSE]
message("  GWAS peaks (cov-weighted): ", nrow(gwas_wmat))

# ===========================================================
# Step 9: Cov-weighted specificity (GWAS peaks)
# ===========================================================
message("\n--- Cov-weighted specificity summary ---")
if (nrow(gwas_wmat) > 0) {
  high_thresh_cw <- as.numeric(quantile(gwas_wmat, 0.9))
  mid_thresh_cw  <- as.numeric(quantile(gwas_wmat, 0.5))
  
  specificity_summary_cw <- data.frame(
    peak   = rownames(gwas_wmat),
    n_high = rowSums(gwas_wmat > high_thresh_cw),
    n_mid  = rowSums(gwas_wmat > mid_thresh_cw & gwas_wmat <= high_thresh_cw),
    n_low  = rowSums(gwas_wmat <= mid_thresh_cw)
  )
  specificity_summary_cw$dominant_ratio <- apply(gwas_wmat, 1, function(x) {
    xs <- sort(x, decreasing = TRUE)
    if (length(xs) < 2) NA_real_ else xs[1] / xs[2]
  })
  specificity_summary_cw$maximum_cell   <- apply(gwas_wmat, 1, function(x) paste(names(x)[x == max(x)], collapse = ", "))
  high_cells_cw <- apply(gwas_wmat, 1, function(x) names(x)[x > high_thresh_cw])
  specificity_summary_cw$high_cells       <- sapply(high_cells_cw, paste, collapse = ", ")
  specificity_summary_cw$is_max_dominant  <- specificity_summary_cw$dominant_ratio > 2
  specificity_summary_cw$multiple_high    <- lengths(high_cells_cw) > 1
  
  peak_gr_cw <- GRanges(rownames(specificity_summary_cw))
  hits_cw    <- findOverlaps(peak_gr_cw, snps)
  snp_list_cw  <- tapply(mcols(snps)$names[subjectHits(hits_cw)],  queryHits(hits_cw), paste, collapse = ",")
  locus_vec_cw <- tapply(mcols(snps)$signal[subjectHits(hits_cw)], queryHits(hits_cw), function(x) paste(unique(x), collapse = ", "))
  specificity_summary_cw$gwas_snps[as.numeric(names(snp_list_cw))]   <- snp_list_cw
  specificity_summary_cw$gwas_locus[as.numeric(names(locus_vec_cw))] <- locus_vec_cw
  
  write.table(specificity_summary_cw, file.path(results_dir, "specificity_summary_covweighted.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  specificity_summary_cw <- data.frame()
  message("  Skipping cov-weighted specificity summaries (no GWAS-overlapping peaks)")
}

# ===========================================================
# Step 10: Cov-weighted dominance  genome-wide vs GWAS
# ===========================================================
message("\n--- Cov-weighted dominance: genome-wide vs GWAS ---")
p_dom_cw <- dominance_dodge_plot(wmat, gwas_wmat, title_label = "CW-L2mat")
ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_dominance_covweighted.pdf")),
       p_dom_cw, width = 8, height = 6)

# Cov-weighted heatmap
if (nrow(gwas_wmat) > 1) {
  pdf(file.path(plots_dir, paste0(se_name, "_", trait_name, "_gwas_covweighted_cor.pdf")), width = 8, height = 7)
  heatmap(cor(gwas_wmat), main = paste0(se_name, " x ", trait_name, " - GWAS cov-weighted correlation"), margins = c(10, 10))
  dev.off()
} else {
  message("  Skipping GWAS cov-weighted correlation heatmap (", nrow(gwas_wmat), " GWAS peaks)")
}

# ===========================================================
# Step 11: Multi-panel summary plot (cov-weighted)
# ===========================================================
message("\n--- Multi-panel summary plot ---")
if (nrow(specificity_summary_cw)) {
  prop_df_cw <- specificity_summary_cw %>%
    summarise(total_high = sum(n_high), total_mid = sum(n_mid), total_low = sum(n_low)) %>%
    pivot_longer(everything(), names_to = "category", values_to = "count") %>%
    mutate(prop = count / sum(count))
  prop_df_cw$category <- factor(prop_df_cw$category,
                                levels = c("total_high", "total_mid", "total_low"))
  
  pA <- ggplot(prop_df_cw, aes(x = "GWAS peaks", y = prop, fill = category)) +
    geom_bar(stat = "identity", width = 0.5) + theme_bw() +
    labs(x = "", y = "Proportion", fill = "Specificity", title = "A. Overall Specificity Distribution")
  
  pB <- ggplot(specificity_summary_cw, aes(x = n_high)) +
    geom_histogram(binwidth = 1, boundary = -1, fill = "steelblue") + theme_bw() +
    labs(title = "B. # High-Specific Cell Types per Peak", x = "# High-specific cells", y = "# Peaks")
  
  pC <- ggplot(specificity_summary_cw, aes(x = factor(n_high), y = dominant_ratio)) +
    geom_boxplot(fill = "darkgreen") + theme_bw() +
    labs(title = "C. Dominant Ratio by # High Cells", x = "N high-specific cells", y = "Dominant ratio")
  
  all_cells <- unique(unlist(strsplit(specificity_summary_cw$high_cells, ", ")))
  high_mat <- sapply(all_cells, function(cell) as.integer(grepl(cell, specificity_summary_cw$high_cells)))
  rownames(high_mat) <- specificity_summary_cw$peak
  pD <- Heatmap(high_mat, name = "High", show_row_names = FALSE, show_column_names = TRUE,
                cluster_rows = TRUE, cluster_columns = TRUE, col = c("white", "red"),
                column_title = "D. High-specific Cell Types per Peak")
  
  pE <- Heatmap(gwas_wmat, name = "CW-L2", show_row_names = FALSE, show_column_names = TRUE,
                cluster_rows = TRUE, cluster_columns = TRUE, col = c("white", "blue"),
                column_title = "E. Cov-Weighted GWAS Peak Matrix")
  
  pdf(file.path(plots_dir, paste0(se_name, "_", trait_name, "_summary_panels_ABC.pdf")), width = 16, height = 5)
  grid.arrange(pA, pB, pC, ncol = 3)
  dev.off()
  
  pdf(file.path(plots_dir, paste0(se_name, "_", trait_name, "_summary_panel_D.pdf")), width = 10, height = 8)
  draw(pD)
  dev.off()
  
  pdf(file.path(plots_dir, paste0(se_name, "_", trait_name, "_summary_panel_E.pdf")), width = 10, height = 8)
  draw(pE)
  dev.off()
} else {
  message("  Skipping cov-weighted summary panels (no GWAS-overlapping peaks)")
}

# ===========================================================
# Step 12: Global rank-based enrichment (on L2 matrix, NOT cov-weighted)
# ===========================================================
message("\n--- Global rank-based enrichment ---")
nmat_r <- apply(mat_l2, 2, function(x) rank(x, ties.method = "min"))
nse_rank <- SummarizedExperiment(nmat_r, rowRanges = rowRanges(se))
gwas_rank_idx <- unique(subjectHits(findOverlaps(snps, rowRanges(nse_rank))))
gwas_rank_mat <- if (length(gwas_rank_idx)) assay(nse_rank[gwas_rank_idx, , drop = FALSE], 1) else nmat_r[0, , drop = FALSE]
if (nrow(gwas_rank_mat) > 0) {
  N <- nrow(nse_rank)
  n <- nrow(gwas_rank_mat)
  exp_mean <- (N + 1) / 2
  se_z     <- sqrt((N^2 - 1) / (12 * n))
  
  global_enrichment <- data.frame(
    observed_mean_rank = colMeans(gwas_rank_mat),
    stringsAsFactors = FALSE
  )
  global_enrichment$n_peaks  <- n
  global_enrichment$exp_rank <- exp_mean
  global_enrichment$obs.exp  <- global_enrichment$observed_mean_rank / exp_mean
  global_enrichment$z        <- (global_enrichment$observed_mean_rank - exp_mean) / se_z
  global_enrichment$p        <- pnorm(global_enrichment$z, lower.tail = FALSE)
  global_enrichment$fdr      <- p.adjust(global_enrichment$p, method = "BH")
  global_enrichment$positive <- global_enrichment$z > 0 & global_enrichment$fdr < 0.05
  global_enrichment$celltype <- rownames(global_enrichment)
  
  write.table(global_enrichment, file.path(results_dir, "global_enrichment.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  p_enrich <- ggplot(global_enrichment, aes(x = reorder(celltype, z), y = z, fill = positive)) +
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    labs(x = "Cell type", y = "Z-score (rank enrichment)",
         title = paste0(se_name, " x ", trait_name, " - Global GWAS Enrichment"),
         fill = "FDR < 0.05")
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_global_enrichment.pdf")),
         p_enrich, width = 8, height = 6)
  
  # FDR bar plot
  p_fdr <- ggplot(global_enrichment, aes(x = reorder(celltype, -log10(fdr)),
                                         y = -log10(fdr), fill = positive)) +
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    labs(x = "Cell type", y = "-log10(FDR)",
         title = paste0(se_name, " x ", trait_name, " - Global Enrichment (FDR)"),
         fill = "FDR < 0.05")
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_global_enrichment_fdr.pdf")),
         p_fdr, width = 8, height = 6)
  
} else {
  global_enrichment <- data.frame()
  message("  Skipping global enrichment (no GWAS-overlapping peaks)")
}

# ===========================================================
# Step 12b: Global rank-based enrichment (QN L2)
# ===========================================================
message("\n--- Global rank-based enrichment (QN L2) ---")
nmat_qn_r <- apply(mat_qn_l2, 2, function(x) rank(x, ties.method = "min"))
nse_qn_rank <- SummarizedExperiment(nmat_qn_r, rowRanges = rowRanges(se))
gwas_qn_rank_idx <- unique(subjectHits(findOverlaps(snps, rowRanges(nse_qn_rank))))
gwas_qn_rank_mat <- if (length(gwas_qn_rank_idx)) assay(nse_qn_rank[gwas_qn_rank_idx, , drop = FALSE], 1) else nmat_qn_r[0, , drop = FALSE]
if (nrow(gwas_qn_rank_mat) > 0) {
  N_qn <- nrow(nse_qn_rank)
  n_qn <- nrow(gwas_qn_rank_mat)
  exp_mean_qn <- (N_qn + 1) / 2
  se_z_qn     <- sqrt((N_qn^2 - 1) / (12 * n_qn))
  
  global_enrichment_qn <- data.frame(
    observed_mean_rank = colMeans(gwas_qn_rank_mat),
    stringsAsFactors = FALSE
  )
  global_enrichment_qn$n_peaks  <- n_qn
  global_enrichment_qn$exp_rank <- exp_mean_qn
  global_enrichment_qn$obs.exp  <- global_enrichment_qn$observed_mean_rank / exp_mean_qn
  global_enrichment_qn$z        <- (global_enrichment_qn$observed_mean_rank - exp_mean_qn) / se_z_qn
  global_enrichment_qn$p        <- pnorm(global_enrichment_qn$z, lower.tail = FALSE)
  global_enrichment_qn$fdr      <- p.adjust(global_enrichment_qn$p, method = "BH")
  global_enrichment_qn$positive <- global_enrichment_qn$z > 0 & global_enrichment_qn$fdr < 0.05
  global_enrichment_qn$celltype <- rownames(global_enrichment_qn)
  
  write.table(global_enrichment_qn, file.path(qn_results_dir, "global_enrichment_qn_l2.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Z-score
  p_enrich_qn <- ggplot(global_enrichment_qn, aes(x = reorder(celltype, z), y = z, fill = positive)) +
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    labs(x = "Cell type", y = "Z-score (rank enrichment)",
         title = paste0(se_name, " x ", trait_name, " - Global Enrichment QN-L2 (Z)"),
         fill = "FDR < 0.05")
  ggsave(file.path(qn_plots_dir, paste0(se_name, "_", trait_name, "_global_enrichment_qn_l2.pdf")),
         p_enrich_qn, width = 8, height = 6)
  
  # FDR
  p_fdr_qn <- ggplot(global_enrichment_qn, aes(x = reorder(celltype, -log10(fdr)),
                                               y = -log10(fdr), fill = positive)) +
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    labs(x = "Cell type", y = "-log10(FDR)",
         title = paste0(se_name, " x ", trait_name, " - Global Enrichment QN-L2 (FDR)"),
         fill = "FDR < 0.05")
  ggsave(file.path(qn_plots_dir, paste0(se_name, "_", trait_name, "_global_enrichment_qn_l2_fdr.pdf")),
         p_fdr_qn, width = 8, height = 6)
  
}