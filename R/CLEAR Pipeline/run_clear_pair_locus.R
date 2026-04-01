#!/usr/bin/env Rscript
# ============================================================
# run_clear_pair_locus.R
# One snATAC SE x GWAS trait pair:
# - L2 specificity + dominance on cell type and lineage (add lineage before dominance)
# - Global rank-based enrichment (Cell type specific)
# - Locus-level rank concordance with matched-bin permutation null (test whether locus-level peaks converge on a shared program)
#   (no QN / cov-weighted matrices)
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
# Arguments
# ---------------------------
args <- commandArgs(TRUE)
if (length(args) < 4) {
  stop(
    "Usage: Rscript run_clear_pair_locus.R <se_name> <se_path> <trait_name> <trait_path> [gtf_path] [n_perm] [k_lineage] [processed_dir]"
  )
}

se_name    <- args[1]
se_path    <- args[2]
trait_name <- args[3]
trait_path <- args[4]
gtf_path   <- if (length(args) >= 5) args[5] else "/working/lab_jonathb/lefeiW/projects/CLEAR_data/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz"
n_perm     <- if (length(args) >= 6) as.integer(args[6]) else 100L
k_lineage  <- if (length(args) >= 7) as.integer(args[7]) else 4L # 
processed_dir <- if (length(args) >= 8) args[8] else file.path(dirname(dirname(getwd())), "processed_data")

message("\n========================================")
message("CLEAR analysis (locus): ", se_name, " x ", trait_name)
message("  SE path:    ", se_path)
message("  Trait path: ", trait_path)
message("  GTF path:   ", gtf_path)
message("  n_perm:     ", n_perm)
message("  k_lineage:  ", k_lineage)
message("  processed_dir: ", processed_dir)
message("========================================\n")

# ---------------------------
# Output directories
# ---------------------------
results_dir    <- "results"
plots_dir      <- "plots"
data_dir       <- "data"

for (d in c(results_dir, plots_dir, data_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Helpers
# ---------------------------
prepare_clear_se <- function(se, genome = "hg38") {
  rd <- as.data.frame(rowData(se))
  
  has_rowdata_coords <- all(c("seqnames", "start", "end") %in% names(rd))
  
  if (has_rowdata_coords) {
    # ArchR-style object: coordinates stored in rowData columns.
    gr <- GRanges(
      seqnames = rd$seqnames,
      ranges   = IRanges::IRanges(start = as.integer(rd$start), end = as.integer(rd$end))
    )
    
    extra_cols <- setdiff(names(rd), c("seqnames", "start", "end", "width", "strand"))
    if (length(extra_cols)) mcols(gr) <- rd[, extra_cols, drop = FALSE]
  } else if (inherits(rowRanges(se), "GenomicRanges")) {
    # Some inputs keep genomic coordinates only in rowRanges.
    gr <- rowRanges(se)
    if (ncol(rd) > 0) {
      for (nm in names(rd)) {
        if (!nm %in% names(mcols(gr))) {
          mcols(gr)[[nm]] <- rd[[nm]]
        }
      }
    }
  } else {
    stop(
      "Input SE has no usable genomic coordinates. Provide either rowData seqnames/start/end or a valid rowRanges(se)."
    )
  }
  
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

compute_metrics_se_no_cov <- function(se, assay_name = NULL, keep_top_prop = NULL) {
  if (is.null(assay_name)) assay_name <- assayNames(se)[1]
  
  mat <- assay(se, assay_name)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  
  rn <- rownames(se); if (is.null(rn)) rn <- paste0("peak_", seq_len(nrow(mat)))
  cn <- colnames(se); if (is.null(cn)) cn <- paste0("ct_", seq_len(ncol(mat)))
  dimnames(mat) <- list(rn, cn)
  
  row_sums <- rowSums(mat)
  if (!is.null(keep_top_prop)) {
    thr <- stats::quantile(row_sums, probs = 1 - keep_top_prop, na.rm = TRUE)
    keep <- row_sums > thr
    se <- se[keep, , drop = FALSE]
    mat <- mat[keep, , drop = FALSE]
    rn <- rn[keep]
    dimnames(mat) <- list(rn, cn)
  }
  
  raw_l2 <- t(apply(mat, 1, safe_l2)); dimnames(raw_l2) <- list(rn, cn)
  
  cor_weight <- function(m) {
    R <- corpcor::cor.shrink(m)
    w <- t(solve(R) %*% t(m))
    dimnames(w) <- dimnames(m)
    w
  }
  
  cor_raw <- cor_weight(mat)
  cor_raw_l2 <- cor_weight(raw_l2)
  
  assays(se)[["raw"]] <- mat
  assays(se)[["raw_l2"]] <- raw_l2
  assays(se)[["raw_l2_rank"]] <- add_rank_desc(raw_l2)
  assays(se)[["cor_raw"]] <- cor_raw
  assays(se)[["cor_raw_l2"]] <- cor_raw_l2
  
  metadata(se)$CLEAR <- list(
    source_assay = assay_name,
    keep_top_prop = keep_top_prop,
    n_peaks = nrow(se),
    n_celltypes = ncol(se),
    use_l2_only = TRUE
  )
  
  se
}

PLOT_ASSAYS <- c("raw", "raw_l2", "cor_raw", "cor_raw_l2")

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
      name = "Correlation",
      column_title = paste0(tissue, " - ", a, " correlation"),
      clustering_method_rows = "ward.D2",
      clustering_method_columns = "ward.D2",
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      row_names_gp = grid::gpar(fontsize = 7),
      column_names_gp = grid::gpar(fontsize = 7),
      col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    )
    ComplexHeatmap::draw(ht)
    dev.off()
  }
  invisible(out)
}

plot_assay_densities <- function(se, tissue, plot_dir = "plots", lineage_palette = NULL) {
  out <- file.path(plot_dir, tissue)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  avail <- intersect(PLOT_ASSAYS, assayNames(se))
  long_list <- lapply(avail, function(a) {
    mat <- assay(se, a)
    if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
    df <- as.data.frame(mat)
    df <- tidyr::pivot_longer(df, cols = everything(), names_to = "cell_type", values_to = "value")
    df$assay <- a
    df
  })
  long <- do.call(rbind, long_list)
  long$assay <- factor(long$assay, levels = avail)
  
  ct_names <- unique(long$cell_type)
  if (!is.null(lineage_palette) && all(ct_names %in% names(lineage_palette))) {
    colour_scale <- scale_color_manual(values = lineage_palette)
  } else {
    n_ct_dens <- length(ct_names)
    if (n_ct_dens <= 10) {
      colour_scale <- ggsci::scale_color_npg()
    } else if (n_ct_dens <= 20) {
      colour_scale <- ggsci::scale_color_d3(palette = "category20")
    } else {
      colour_scale <- ggsci::scale_color_igv()
    }
  }
  
  p <- ggplot(long, aes(x = value, colour = cell_type)) +
    geom_density() +
    facet_wrap(~assay, scales = "free", ncol = 4) +
    colour_scale +
    theme_bw() +
    theme(legend.position = "bottom", strip.text = element_text(size = 9)) +
    labs(
      title = paste0(tissue, " - per-cell-type distributions"),
      x = "Value", y = "Density", colour = "Cell type"
    )
  
  fname <- file.path(out, paste0(tissue, "_assay_densities.pdf"))
  ggsave(fname, p, width = 16, height = max(12, ceiling(length(avail) / 4) * 3))
  invisible(fname)
}

compute_dominance <- function(mat, ratio_thresh = 2) {
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  if (n_rows == 0L) {
    return(data.frame(maximum_ct = character(), ratio = numeric(), is_dominant = logical(), peak = character()))
  }
  if (is.null(colnames(mat))) colnames(mat) <- paste0("ct", seq_len(n_cols))
  
  if (n_cols < 2L) {
    return(data.frame(
      maximum_ct = rep(colnames(mat)[1], n_rows),
      ratio = rep(NA_real_, n_rows),
      is_dominant = rep(FALSE, n_rows),
      peak = rownames(mat)
    ))
  }
  
  cell_names <- colnames(mat)
  max_idx <- max.col(mat, ties.method = "first")
  max_val <- mat[cbind(seq_len(n_rows), max_idx)]
  
  mat2 <- mat
  mat2[cbind(seq_len(n_rows), max_idx)] <- -Inf
  second_val <- mat2[cbind(seq_len(n_rows), max.col(mat2, ties.method = "first"))]
  
  ratio <- ifelse(second_val == 0, Inf, max_val / second_val)
  
  data.frame(
    maximum_ct = cell_names[max_idx],
    ratio = ratio,
    is_dominant = ratio > ratio_thresh,
    peak = rownames(mat),
    stringsAsFactors = FALSE
  )
}

DEFAULT_RATIO_THRESH <- 2

build_lineage_map <- function(mat_l2, k_lineage = 4L) {
  if (inherits(mat_l2, "Matrix")) mat_l2 <- as.matrix(mat_l2)
  storage.mode(mat_l2) <- "double"
  if (is.null(colnames(mat_l2))) {
    colnames(mat_l2) <- paste0("ct", seq_len(ncol(mat_l2)))
  }

  cell_cor <- cor(mat_l2, method = "pearson")
  hc <- hclust(as.dist(1 - cell_cor), method = "average")
  lineage_ids <- cutree(hc, k = k_lineage)
  lineage_labels <- tapply(
    names(lineage_ids),
    lineage_ids,
    function(x) paste(sort(x), collapse = ":")
  )
  cell_lineage <- setNames(as.character(lineage_labels[as.character(lineage_ids)]), names(lineage_ids))

  palette <- build_lineage_palette(lineage_ids)

  list(
    cell_cor = cell_cor,
    hc = hc,
    lineage_ids = lineage_ids,
    lineage_labels = lineage_labels,
    cell_lineage = cell_lineage,
    k_lineage = k_lineage,
    palette = palette
  )
}

# Hierarchical colour palette: one master hue per lineage, shades per cell type
# Supports up to 8 lineages x 8 cell types each
build_lineage_palette <- function(lineage_ids) {
  # Master hues (up to 8 lineages)
  master_hues <- c(
    "#2166AC",  # blue
    "#B2182B",  # red
    "#1B7837",  # green
    "#E08214",  # orange
    "#6A3D9A",  # purple
    "#A6761D",  # brown
    "#E7298A",  # pink
    "#66C2A5"   # teal
  )

  groups <- split(names(lineage_ids), lineage_ids)
  n_groups <- length(groups)
  if (n_groups > length(master_hues)) {
    master_hues <- c(master_hues, grDevices::hcl.colors(n_groups - length(master_hues), "Set2"))
  }

  pal <- character()
  for (i in seq_along(groups)) {
    cts <- sort(groups[[i]])
    n_ct <- length(cts)
    if (n_ct == 1L) {
      shades <- master_hues[i]
    } else {
      # Generate shades from dark to light within each master hue
      base_rgb <- grDevices::col2rgb(master_hues[i]) / 255
      alphas <- seq(1.0, 0.35, length.out = n_ct)
      shades <- vapply(alphas, function(a) {
        blended <- base_rgb * a + 1 * (1 - a)  # blend toward white
        grDevices::rgb(blended[1], blended[2], blended[3])
      }, character(1))
    }
    names(shades) <- cts
    pal <- c(pal, shades)
  }

  pal
}

plot_lineage_dendrogram <- function(lineage_obj, se_name, plot_dir = "plots") {
  fname <- file.path(plot_dir, paste0(se_name, "_celltype_lineage_dendrogram.pdf"))
  hc <- lineage_obj$hc
  k_lineage <- lineage_obj$k_lineage

  grDevices::pdf(fname, width = 10, height = 6)
  op <- par(mar = c(10, 4, 4, 2))
  plot(
    hc,
    main = paste0(se_name, " - Cell-type lineage dendrogram (k=", k_lineage, ")"),
    xlab = "",
    sub = "Distance = 1 - Pearson correlation on genome-wide raw_l2",
    cex = 0.7
  )
  rect.hclust(hc, k = k_lineage, border = "firebrick")
  par(op)
  grDevices::dev.off()

  invisible(fname)
}

plot_celltype_stackbars <- function(
    se_overlap,
    value_metric = "raw_l2",
    cumvar_threshold = 0.8,
    max_k = 3,
    outfile,
    max_peaks = 200,
    lineage_palette = NULL) {
  stopifnot(!missing(outfile))
  
  if (!value_metric %in% assayNames(se_overlap)) {
    warning("Assay '", value_metric, "' not found; skipping stack bar plots.")
    return(invisible(FALSE))
  }
  
  val_mat <- SummarizedExperiment::assay(se_overlap, value_metric)
  if (inherits(val_mat, "Matrix")) val_mat <- as.matrix(val_mat)
  storage.mode(val_mat) <- "double"
  
  n_celltypes <- ncol(val_mat)
  min_l2 <- 1 / sqrt(n_celltypes / 2)
  val_sq <- val_mat^2
  
  peak_info <- lapply(seq_len(nrow(val_mat)), function(i) {
    row_sq <- val_sq[i, ]
    row_l2 <- val_mat[i, ]
    ord <- order(row_sq, decreasing = TRUE)
    cum_sq <- cumsum(row_sq[ord])
    
    k <- which(cum_sq >= cumvar_threshold)[1]
    if (is.na(k) || k > max_k) return(NULL)
    
    top_k_idx <- ord[seq_len(k)]
    if (row_l2[top_k_idx[k]] < min_l2) return(NULL)
    
    list(
      idx = i,
      k = k,
      cum_var = cum_sq[k],
      top_cells = colnames(val_mat)[top_k_idx],
      min_l2 = row_l2[top_k_idx[k]]
    )
  })
  peak_info <- Filter(Negate(is.null), peak_info)
  
  if (!length(peak_info)) {
    warning(
      "No peaks meet criteria (top <=", max_k,
      " cells cumulative L2^2 >= ", cumvar_threshold,
      ", min L2 >= ", round(min_l2, 3), ")"
    )
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
  
  rd <- as.data.frame(rowData(se_overlap))
  peakid <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se_overlap)
  ccv_label <- if ("CCVs" %in% names(rd)) rd$CCVs else ""
  sig_label <- if ("signal" %in% names(rd)) rd$signal else ""
  base_label <- ifelse(nzchar(ccv_label), ccv_label, peakid)
  labels <- ifelse(nzchar(sig_label), paste0(base_label, " (", sig_label, ")"), base_label)
  
  dom_ct <- apply(val_mat[sel, , drop = FALSE], 1, function(x) names(which.max(x)))
  
  long_list <- lapply(seq_along(sel), function(j) {
    i <- sel[j]
    row_sq <- val_sq[i, ]
    ct_ord <- order(row_sq, decreasing = TRUE)
    data.frame(
      peak_label = labels[i],
      celltype = colnames(val_sq)[ct_ord],
      value = row_sq[ct_ord],
      stack_rank = seq_along(ct_ord),
      k = sel_k[j],
      dom_ct = dom_ct[j],
      max_sq = max(row_sq),
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, long_list)
  
  peak_order <- long %>%
    dplyr::distinct(peak_label, k, dom_ct, max_sq) %>%
    dplyr::arrange(k, dom_ct, -max_sq) %>%
    dplyr::pull(peak_label)
  long$peak_label <- factor(long$peak_label, levels = rev(peak_order))
  
  ct_totals <- tapply(long$value, long$celltype, sum, na.rm = TRUE)
  ct_levels <- names(sort(ct_totals, decreasing = TRUE))
  long$celltype <- factor(long$celltype, levels = ct_levels)
  long <- long %>% dplyr::arrange(peak_label, -value)
  
  threshold_pct <- cumvar_threshold * 100
  n_sel <- length(sel)
  n_ct <- length(ct_levels)
  
  if (!is.null(lineage_palette) && all(ct_levels %in% names(lineage_palette))) {
    fill_scale <- scale_fill_manual(values = lineage_palette)
  } else if (n_ct <= 10) {
    fill_scale <- ggsci::scale_fill_npg()
  } else if (n_ct <= 20) {
    fill_scale <- ggsci::scale_fill_d3(palette = "category20")
  } else {
    fill_scale <- ggsci::scale_fill_igv()
  }
  
  p <- ggplot(long, aes(x = peak_label, y = value, fill = celltype, group = factor(stack_rank))) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.85) +
    geom_hline(yintercept = cumvar_threshold, linetype = "dashed", color = "black", linewidth = 0.3) +
    coord_flip() +
    fill_scale +
    labs(
      title = paste0(
        "GWAS peaks: top <=", max_k, " cell types explain >=", threshold_pct,
        "% variance (min L2 >= ", round(min_l2, 3), ")"
      ),
      subtitle = paste0(n_sel, " peaks selected"),
      x = "Peak / CCV",
      y = "Variance explained (L2^2, cumulative = 1.0)",
      fill = "Cell type"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 4),
      axis.text.x = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 8),
      plot.title = element_text(size = 10)
    )
  
  pdf_height <- max(6, n_sel * 0.12 + 3)
  grDevices::pdf(outfile, width = 12, height = min(pdf_height, 50))
  print(p)
  grDevices::dev.off()
  
  message("  Saved stackbar plot (", n_sel, " peaks): ", outfile)
  invisible(TRUE)
}

compute_specificity_summary <- function(mat_l2, cell_cor = NULL, ratio_thresh = 2, threshold_sq = 0.8) {
  if (is.null(cell_cor)) cell_cor <- cor(mat_l2, use = "pairwise.complete.obs")
  results <- apply(mat_l2, 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    x_sorted <- x[ord]
    cells_sorted <- colnames(mat_l2)[ord]
    cum_sq <- cumsum(x_sorted^2)
    k <- which(cum_sq >= threshold_sq)[1]
    if (is.na(k)) k <- length(x_sorted)
    top_cells <- cells_sorted[seq_len(k)]
    
    maximum_cell <- cells_sorted[1]
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
      maximum_cell = maximum_cell,
      maximum_score = as.numeric(maximum_score),
      dominant_ratio = as.numeric(dominant_ratio),
      is_dominant = as.character(dominant_ratio > ratio_thresh),
      k_cells = as.numeric(k),
      top_k_cells = paste(top_cells, collapse = ","),
      coherence = as.numeric(coherence)
    )
  })
  
  out <- as.data.frame(t(results), stringsAsFactors = FALSE)
  out$peak <- rownames(mat_l2)
  out$maximum_score <- as.numeric(out$maximum_score)
  out$dominant_ratio <- as.numeric(out$dominant_ratio)
  out$is_dominant <- as.logical(out$is_dominant)
  out$k_cells <- as.integer(out$k_cells)
  out$coherence <- as.numeric(out$coherence)
  rownames(out) <- NULL
  out
}

dominance_dodge_plot <- function(genome_mat, gwas_mat,
                                 ratio_thresh = DEFAULT_RATIO_THRESH,
                                 title_label  = "",
                                 lineage_palette = NULL) {
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
  
  p <- ggplot(df, aes(x = maximum_ct, y = proportion, fill = set)) +
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

  # Add lineage-coloured annotation strip on the y-axis
  if (!is.null(lineage_palette)) {
    ct_cols <- lineage_palette[levels(df$maximum_ct)]
    ct_cols <- ct_cols[!is.na(ct_cols)]
    if (length(ct_cols)) {
      p <- p + theme(axis.text.y = element_text(colour = ct_cols))
    }
  }
  p
}

compute_dominance_lineage <- function(mat, cell_lineage, ratio_thresh = 2) {
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  if (n_rows == 0L) {
    return(data.frame(maximum_lineage = character(), ratio = numeric(), is_dominant = logical(), peak = character()))
  }
  if (is.null(colnames(mat))) colnames(mat) <- paste0("ct", seq_len(n_cols))

  lineage_vec <- cell_lineage[colnames(mat)]
  missing_lineage <- is.na(lineage_vec) | !nzchar(lineage_vec)
  if (any(missing_lineage)) {
    lineage_vec[missing_lineage] <- colnames(mat)[missing_lineage]
  }

  lineage_levels <- unique(as.character(lineage_vec))
  mat_sq <- mat^2
  lineage_mat <- sapply(lineage_levels, function(lg) {
    rowSums(mat_sq[, lineage_vec == lg, drop = FALSE])
  })
  lineage_mat <- as.matrix(lineage_mat)
  if (is.null(dim(lineage_mat))) {
    lineage_mat <- matrix(lineage_mat, ncol = 1)
    colnames(lineage_mat) <- lineage_levels[1]
  }
  rownames(lineage_mat) <- rownames(mat)

  if (ncol(lineage_mat) < 2L) {
    return(data.frame(
      maximum_lineage = rep(colnames(lineage_mat)[1], nrow(lineage_mat)),
      ratio = rep(NA_real_, nrow(lineage_mat)),
      is_dominant = rep(FALSE, nrow(lineage_mat)),
      peak = rownames(lineage_mat),
      stringsAsFactors = FALSE
    ))
  }

  max_idx <- max.col(lineage_mat, ties.method = "first")
  max_val <- lineage_mat[cbind(seq_len(nrow(lineage_mat)), max_idx)]

  lineage_mat2 <- lineage_mat
  lineage_mat2[cbind(seq_len(nrow(lineage_mat2)), max_idx)] <- -Inf
  second_val <- lineage_mat2[cbind(seq_len(nrow(lineage_mat2)), max.col(lineage_mat2, ties.method = "first"))]

  ratio <- ifelse(second_val == 0, Inf, max_val / second_val)

  data.frame(
    maximum_lineage = colnames(lineage_mat)[max_idx],
    ratio = ratio,
    is_dominant = ratio > ratio_thresh,
    peak = rownames(lineage_mat),
    stringsAsFactors = FALSE
  )
}

dominance_lineage_dodge_plot <- function(genome_mat, gwas_mat, cell_lineage,
                                         ratio_thresh = DEFAULT_RATIO_THRESH,
                                         title_label  = "") {
  dom_genome <- compute_dominance_lineage(genome_mat, cell_lineage, ratio_thresh = ratio_thresh)
  dom_gwas <- compute_dominance_lineage(gwas_mat, cell_lineage, ratio_thresh = ratio_thresh)

  n_dom_genome <- sum(dom_genome$is_dominant)
  n_dom_gwas <- sum(dom_gwas$is_dominant)

  df <- dplyr::bind_rows(
    dom_genome %>% dplyr::filter(is_dominant) %>% dplyr::count(maximum_lineage) %>%
      dplyr::mutate(proportion = if (n_dom_genome > 0) n / n_dom_genome else 0, set = "Genome-wide"),
    dom_gwas %>% dplyr::filter(is_dominant) %>% dplyr::count(maximum_lineage) %>%
      dplyr::mutate(proportion = if (n_dom_gwas > 0) n / n_dom_gwas else 0, set = "GWAS")
  )

  if (!nrow(df)) {
    df <- data.frame(
      maximum_lineage = character(),
      n = integer(),
      proportion = numeric(),
      set = character(),
      stringsAsFactors = FALSE
    )
  }

  df <- tidyr::complete(df, maximum_lineage, set,
                        fill = list(n = 0L, proportion = 0))

  df$set <- factor(df$set, levels = c("Genome-wide", "GWAS"))
  df$bar_label <- ifelse(df$n > 0,
                         paste0(df$n, " (", round(df$proportion * 100, 1), "%)"), "")

  ct_order <- df %>% dplyr::filter(set == "GWAS") %>%
    dplyr::arrange(proportion) %>% dplyr::pull(maximum_lineage)
  df$maximum_lineage <- factor(df$maximum_lineage, levels = ct_order)

  ggplot(df, aes(x = maximum_lineage, y = proportion, fill = set)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = bar_label),
              position = position_dodge(width = 0.8), hjust = -0.05, size = 2.8) +
    coord_flip() +
    scale_fill_manual(values = c("Genome-wide" = "grey60", "GWAS" = "firebrick")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    theme_bw() +
    labs(
      title = paste0("Dominant lineages (ratio > ", ratio_thresh, ") ", title_label),
      subtitle = paste0("Genome-wide: ", nrow(genome_mat), " peaks (", n_dom_genome,
                        " dominant);  GWAS: ", nrow(gwas_mat), " peaks (", n_dom_gwas, " dominant)"),
      x = "Lineage",
      y = "Proportion within dominant peaks",
      fill = "Peak set"
    )
}

get_tss <- function(gtf) {
  genes <- gtf[gtf$type == "gene"]
  genes <- genes[genes$gene_type %in% c("lncRNA", "protein_coding")]
  tss_gr <- GenomicRanges::promoters(genes, upstream = 0, downstream = 1)
  S4Vectors::mcols(tss_gr) <- S4Vectors::mcols(tss_gr)[, c("gene_type", "gene_name")]
  tss_gr
}

compute_peak_density <- function(peaks_gr, window = 50000) {
  counts <- countOverlaps(peaks_gr, resize(peaks_gr, width = window, fix = "center"))
  counts - 1
}

make_density_bins <- function(density_vec, nbins = 4) {
  cut(
    density_vec,
    breaks = unique(quantile(density_vec, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)),
    include.lowest = TRUE
  )
}

prepare_peak_annotations <- function(se, gtf, gwas_mat, specificity_summary_GWAS = NULL,
                                     density_window = 50000,
                                     max_locus_distance = 100000) {
  peaks_gr <- rowRanges(se)
  peak_names <- names(peaks_gr)

  tss_gr <- get_tss(gtf)
  nearest_hits <- distanceToNearest(peaks_gr, tss_gr)

  dist_to_tss <- mcols(nearest_hits)$distance
  names(dist_to_tss) <- peak_names[queryHits(nearest_hits)]
  mcols(peaks_gr)$dist_to_tss <- dist_to_tss[peak_names]

  tss_bins <- cut(
    dist_to_tss,
    breaks = c(0, 1000, 5000, 50000, 200000, Inf),
    labels = c("core_prom", "prox_prom", "near_reg", "distal", "long_range"),
    include.lowest = TRUE
  )
  names(tss_bins) <- peak_names
  tss_bins_gwas <- tss_bins[rownames(gwas_mat)]

  density_vec <- compute_peak_density(peaks_gr, window = density_window)
  names(density_vec) <- peak_names

  density_bins <- cut(
    density_vec,
    breaks = unique(quantile(density_vec, probs = seq(0, 1, length.out = 10))),
    include.lowest = TRUE
  )
  names(density_bins) <- peak_names

  joint_bins <- interaction(tss_bins, density_bins, drop = TRUE)
  joint_bins <- as.character(joint_bins)
  names(joint_bins) <- peak_names

  peaks_sorted <- sort(peaks_gr)
  same_chr <- as.character(seqnames(peaks_sorted)[-1]) ==
    as.character(seqnames(peaks_sorted)[-length(peaks_sorted)])
  close_peaks <- diff(start(peaks_sorted)) <= max_locus_distance
  cluster_id <- cumsum(c(1, !(close_peaks & same_chr)))
  locus_ids <- paste0("locus_", cluster_id)
  names(locus_ids) <- names(peaks_sorted)
  locus_ids_full <- locus_ids[peak_names]
  peak_to_locus_full <- setNames(locus_ids_full, peak_names)

  if (!is.null(specificity_summary_GWAS) && "gwas_locus" %in% names(specificity_summary_GWAS)) {
    peak_to_locus_gwas <- specificity_summary_GWAS$gwas_locus
    names(peak_to_locus_gwas) <- specificity_summary_GWAS$peak
  } else {
    peak_to_locus_gwas <- setNames(as.character(rowData(se[rownames(gwas_mat),])$signal), rownames(gwas_mat))
  }

  peak_annot <- data.frame(
    peak = peak_names,
    dist_to_tss = dist_to_tss[peak_names],
    density = density_vec[peak_names],
    tss_bin = as.character(tss_bins[peak_names]),
    density_bin = as.character(density_bins[peak_names]),
    joint_bin = joint_bins[peak_names],
    locus = locus_ids_full,
    stringsAsFactors = FALSE
  )

  list(
    peak_annot = peak_annot,
    peak_to_locus_full = peak_to_locus_full,
    peak_to_locus_gwas = peak_to_locus_gwas,
    tss_bins_full = tss_bins,
    tss_bins_gwas = tss_bins_gwas,
    density_vec = density_vec,
    density_bins = density_bins,
    joint_bins = joint_bins
  )
}

dominance_consensus <- function(specificity_matrix, peak_weights, group_labels) {
  cell_types <- colnames(specificity_matrix)
  groups     <- group_labels[cell_types]

  group_sums <- lapply(unique(groups), function(g) {
    cols <- cell_types[groups == g]
    rowSums(specificity_matrix[, cols, drop = FALSE])
  }) %>%
    setNames(unique(groups)) %>% as.data.frame()

  top_group_per_peak <- apply(group_sums, 1, function(x) names(which.max(x)))
  peak_spec_weight <- apply(group_sums, 1, max)
  peak_weighted_spec <- peak_spec_weight * peak_weights

  group_totals  <- tapply(peak_weighted_spec, top_group_per_peak, sum)
  top_group     <- names(which.max(group_totals))
  consensus     <- group_totals[top_group] / sum(group_totals)

  list(
    consensus    = unname(consensus),
    top_group    = top_group,
    group_totals = sort(group_totals / sum(group_totals), decreasing = TRUE)
  )
}

compute_locus_consensus_concordance_vec <- function(
    gwas_mat,
    full_mat,
    l2_weights,
    peak_to_locus,
    peak_to_locus_full,
    joint_bin_full,
    cell_lineage,
    n_perm = 100
) {
  locus_to_peaks_full <- split(names(peak_to_locus_full), peak_to_locus_full)
  all_bins <- sort(unique(joint_bin_full))

  message("precomputing bin counts per locus (consensus)...")
  locus_bin_mat <- t(sapply(locus_to_peaks_full, function(peaks) {
    tab <- table(factor(joint_bin_full[peaks], levels = all_bins))
    as.numeric(tab)
  }))
  rownames(locus_bin_mat) <- names(locus_to_peaks_full)
  locus_bin_prop <- locus_bin_mat / rowSums(locus_bin_mat)

  gwas_loci_full <- unique(peak_to_locus_full[rownames(gwas_mat)])
  non_gwas_loci <- setdiff(names(locus_to_peaks_full), gwas_loci_full)

  locus_split <- split(rownames(gwas_mat), peak_to_locus[rownames(gwas_mat)])

  group_labels <- setNames(colnames(gwas_mat), colnames(gwas_mat))

  results <- lapply(names(locus_split), function(locus) {
    message(paste0("calculating consensus coherence for ", locus, "..."))
    peaks <- locus_split[[locus]]
    k <- length(peaks)

    if (k < 2) {
      return(data.frame(
        locus = locus, n_peaks = k,
        consensus_fine = NA_real_, consensus_lineage = NA_real_,
        perm_mean = NA_real_, perm_sd = NA_real_,
        z_score = NA_real_, p_value = NA_real_,
        dominant_cell_consensus = NA_character_,
        dominant_lineage_consensus = NA_character_,
        stringsAsFactors = FALSE))
    }

    fine_obs <- dominance_consensus(gwas_mat[peaks, , drop = FALSE], l2_weights[peaks], group_labels)
    fine_obs_conc <- fine_obs$consensus
    fine_top_cell <- fine_obs$top_group

    lineage_obs <- dominance_consensus(gwas_mat[peaks, , drop = FALSE], l2_weights[peaks], cell_lineage)
    lineage_obs_conc <- lineage_obs$consensus
    lineage_top_grp  <- lineage_obs$top_group

    locus_joint_bins <- joint_bin_full[peaks]

    perm_vals <- replicate(n_perm, {
      valid_loci <- non_gwas_loci[sapply(locus_to_peaks_full[non_gwas_loci], length) >= k]
      if (length(valid_loci) == 0) return(NA)

      lb <- table(factor(locus_joint_bins, levels = all_bins))
      lb <- as.numeric(lb) / sum(lb)
      candidate_mat <- locus_bin_prop[valid_loci, , drop = FALSE]
      dists <- rowSums(abs(candidate_mat - matrix(lb, nrow = nrow(candidate_mat), ncol = length(lb), byrow = TRUE)))

      top_n <- min(20, length(valid_loci))
      best_loci <- valid_loci[order(dists)][seq_len(top_n)]
      sampled_locus <- sample(best_loci, 1)
      sampled_peaks <- sample(locus_to_peaks_full[[sampled_locus]], k)

      m <- full_mat[sampled_peaks, , drop = FALSE]
      dominance_consensus(m, l2_weights[sampled_peaks], group_labels)$consensus
    })

    perm_mean <- mean(perm_vals, na.rm = TRUE)
    perm_sd   <- sd(perm_vals, na.rm = TRUE)
    n_perm_eff <- sum(!is.na(perm_vals))
    z_score <- if (is.na(perm_sd) || perm_sd == 0) NA_real_ else (fine_obs_conc - perm_mean) / perm_sd
    p_value <- if (n_perm_eff > 0) sum(perm_vals >= fine_obs_conc, na.rm = TRUE) / n_perm_eff else NA_real_

    data.frame(
      locus = locus, n_peaks = k,
      consensus_fine = fine_obs_conc,
      consensus_lineage = lineage_obs_conc,
      perm_mean = perm_mean, perm_sd = perm_sd,
      z_score = z_score, p_value = p_value,
      dominant_cell_consensus = fine_top_cell,
      dominant_lineage_consensus = lineage_top_grp,
      stringsAsFactors = FALSE)
  })

  do.call(rbind, results)
}

compute_locus_rank_concordance_vec <- function( # Whether peaks in the same GWAS locus show more silimar cell type ranking patterns 
  gwas_mat,
  full_mat,
  peak_to_locus, # signal vector 
  joint_bin_full, # density bin and TSS distance
  cell_lineage,
  n_perm = 100) {
  
  rank_gwas <- t(apply(gwas_mat, 1, rank))
  rank_full <- t(apply(full_mat, 1, rank))
  
  valid_locus <- peak_to_locus[rownames(gwas_mat)]
  keep <- !is.na(valid_locus) & nzchar(valid_locus)
  if (!any(keep)) return(data.frame())
  
  gwas_peaks <- rownames(gwas_mat)[keep]
  locus_split <- split(gwas_peaks, valid_locus[keep])
  full_bins <- split(names(joint_bin_full), joint_bin_full)
  
  results <- lapply(names(locus_split), function(locus) {
    peaks <- locus_split[[locus]]
    k <- length(peaks)
    
    if (k < 2) {
      return(data.frame(
        locus = locus,
        n_peaks = k,
        concordance = NA_real_,
        perm_mean = NA_real_,
        perm_sd = NA_real_,
        z_score = NA_real_,
        p_value = NA_real_,
        dominant_cell_locus = NA_character_,
        dominant_score_locus = NA_real_,
        dominant_lineage = NA_character_,
        lineage_fraction = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    obs_mat <- rank_gwas[peaks, , drop = FALSE]
    obs_cor <- cor(t(obs_mat), method = "spearman")
    obs_conc <- mean(obs_cor[upper.tri(obs_cor)])
    
    mean_profile <- colMeans(gwas_mat[peaks, , drop = FALSE])
    dominant_cell_locus <- names(which.max(mean_profile))
    dominant_score_locus <- max(mean_profile)
    
    lineage_profile <- tapply(mean_profile, cell_lineage[names(mean_profile)], sum)
    dominant_lineage <- names(which.max(lineage_profile))
    lineage_fraction <- max(lineage_profile) / sum(lineage_profile)
    
    locus_joint_bins <- joint_bin_full[peaks]
    locus_bins <- table(locus_joint_bins)
    
    if (any(!names(locus_bins) %in% names(full_bins))) {
      return(data.frame(
        locus = locus,
        n_peaks = k,
        concordance = obs_conc,
        perm_mean = NA_real_,
        perm_sd = NA_real_,
        z_score = NA_real_,
        p_value = NA_real_,
        dominant_cell_locus = dominant_cell_locus,
        dominant_score_locus = dominant_score_locus,
        dominant_lineage = dominant_lineage,
        lineage_fraction = lineage_fraction,
        stringsAsFactors = FALSE
      ))
    }
    
    sampled_matrix <- matrix(NA_character_, nrow = k, ncol = n_perm)
    row_index <- 1
    
    for (bin in names(locus_bins)) {
      n_bin <- as.integer(locus_bins[[bin]])
      pool <- full_bins[[bin]]
      if (is.null(pool) || length(pool) == 0) next
      
      draws <- replicate(
        n_perm,
        sample(pool, n_bin, replace = length(pool) < n_bin),
        simplify = "matrix"
      )
      
      sampled_matrix[row_index:(row_index + n_bin - 1), ] <- draws
      row_index <- row_index + n_bin
    }
    
    perm_vals <- apply(sampled_matrix, 2, function(sampled_peaks) {
      perm_mat <- rank_full[sampled_peaks, , drop = FALSE]
      perm_cor <- cor(t(perm_mat), method = "spearman")
      mean(perm_cor[upper.tri(perm_cor)])
    })
    
    perm_mean <- mean(perm_vals, na.rm = TRUE)
    perm_sd <- sd(perm_vals, na.rm = TRUE)
    n_perm_eff <- sum(!is.na(perm_vals))
    z_score <- ifelse(is.na(perm_sd) || perm_sd == 0, NA_real_, (obs_conc - perm_mean) / perm_sd)
    p_value <- if (n_perm_eff > 0) {
      # Empirical p-value with +1 correction prevents exact zeros.
      (sum(perm_vals >= obs_conc, na.rm = TRUE) + 1) / (n_perm_eff + 1)
    } else {
      NA_real_
    }
    
    data.frame(
      locus = locus,
      n_peaks = k,
      concordance = obs_conc,
      perm_mean = perm_mean,
      perm_sd = perm_sd,
      z_score = z_score,
      p_value = p_value,
      dominant_cell_locus = dominant_cell_locus,
      dominant_score_locus = dominant_score_locus,
      dominant_lineage = dominant_lineage,
      lineage_fraction = lineage_fraction,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, results)
}

addLocusCoherence <- function(se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir) {
  if (!nrow(gwas_mat)) return(data.frame())
  if (!file.exists(gtf_path)) {
    message("  Skipping locus analysis: GTF file not found at ", gtf_path)
    return(data.frame())
  }
  
  peak_to_locus <- setNames(as.character(rowData(se_gwas)$signal), rownames(se_gwas))
  
  gtf <- rtracklayer::import(gtf_path)
  tss_gr <- get_tss(gtf)
  
  peaks_gr <- rowRanges(se)
  nearest_hits <- distanceToNearest(peaks_gr, tss_gr)
  dist2tss <- mcols(nearest_hits)$distance
  
  tss_bin <- cut(
    dist2tss,
    breaks = c(0, 1000, 5000, 50000, 200000, Inf),
    labels = c("core_prom", "prox_prom", "near_reg", "distal", "long_range"),
    include.lowest = TRUE
  )
  names(tss_bin) <- rownames(mat_l2)
  
  density <- compute_peak_density(rowRanges(se))
  density_bin <- make_density_bins(density)
  names(density_bin) <- rownames(mat_l2)
  
  joint_bin_full <- interaction(tss_bin, density_bin, drop = TRUE)
  joint_bin_full <- as.character(joint_bin_full)
  names(joint_bin_full) <- names(tss_bin)

  cell_lineage <- lineage_obj$cell_lineage
  
  locus_results <- compute_locus_rank_concordance_vec(
    gwas_mat = gwas_mat,
    full_mat = mat_l2,
    peak_to_locus = peak_to_locus,
    joint_bin_full = joint_bin_full,
    cell_lineage = cell_lineage,
    n_perm = n_perm
  )
  
  if (!nrow(locus_results)) return(data.frame())
  
  locus_results$significant <- locus_results$p_value < 0.05
  locus_results$category <- dplyr::case_when(
    locus_results$n_peaks == 1 ~ "Single_peak",
    locus_results$z_score > 2 & locus_results$p_value < 0.05 & locus_results$dominant_score_locus > 0.5 ~ "Highly_coherent_specific",
    locus_results$z_score > 2 & locus_results$p_value < 0.05 ~ "Coherent_moderate",
    locus_results$z_score < 1 ~ "Incoherent",
    TRUE ~ "Intermediate"
  )
  
  write.table(locus_results, file.path(results_dir, "locus_concordance_permutation.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  locus_results
}

addLocusCoherence_consensus <- function(se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir, specificity_summary_l2 = NULL) {
  if (!nrow(gwas_mat)) return(data.frame())
  if (!file.exists(gtf_path)) {
    message("  Skipping consensus locus analysis: GTF file not found at ", gtf_path)
    return(data.frame())
  }

  gtf <- rtracklayer::import(gtf_path)

  peak_annotation <- prepare_peak_annotations(
    se = se, gtf = gtf, gwas_mat = gwas_mat,
    specificity_summary_GWAS = specificity_summary_l2
  )

  # Use signal-based locus mapping when specificity_summary not available
  peak_to_locus_gwas <- peak_annotation$peak_to_locus_gwas
  if (is.null(peak_to_locus_gwas) || all(is.na(peak_to_locus_gwas))) {
    peak_to_locus_gwas <- setNames(as.character(rowData(se_gwas)$signal), rownames(se_gwas))
  }

  raw_mat <- assay(se, "raw")
  if (inherits(raw_mat, "Matrix")) raw_mat <- as.matrix(raw_mat)
  l2_weights <- rowMeans(raw_mat)

  cell_lineage <- lineage_obj$cell_lineage

  locus_results_consensus <- compute_locus_consensus_concordance_vec(
    gwas_mat = gwas_mat,
    full_mat = mat_l2,
    l2_weights = l2_weights,
    peak_to_locus = peak_to_locus_gwas,
    peak_to_locus_full = peak_annotation$peak_to_locus_full,
    joint_bin_full = peak_annotation$joint_bins,
    cell_lineage = cell_lineage,
    n_perm = n_perm
  )

  if (!nrow(locus_results_consensus)) return(data.frame())

  locus_results_consensus$significant <- locus_results_consensus$p_value < 0.05
  locus_results_consensus$category <- dplyr::case_when(
    locus_results_consensus$n_peaks == 1 ~ "Single_peak",
    locus_results_consensus$z_score > 2 & locus_results_consensus$p_value < 0.05 &
      locus_results_consensus$consensus_lineage > 0.5 ~ "Highly_coherent_specific",
    locus_results_consensus$z_score > 2 & locus_results_consensus$p_value < 0.05 ~ "Coherent_moderate",
    locus_results_consensus$z_score < 1 ~ "Incoherent",
    TRUE ~ "Intermediate"
  )

  write.table(locus_results_consensus, file.path(results_dir, "locus_concordance_consensus_permutation.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  locus_results_consensus
}

add_gwas_peak_annotation <- function(summary_df, snps) {
  if (!nrow(summary_df)) return(summary_df)
  
  peak_gr <- GRanges(summary_df$peak)
  hits <- findOverlaps(peak_gr, snps)
  snps_names <- if ("names" %in% colnames(mcols(snps))) mcols(snps)$names else names(snps)
  signal_col <- if ("signal" %in% colnames(mcols(snps))) mcols(snps)$signal else rep(NA_character_, length(snps))
  
  snp_list <- tapply(snps_names[subjectHits(hits)], queryHits(hits), paste, collapse = ",")
  locus_vec <- tapply(signal_col[subjectHits(hits)], queryHits(hits), function(x) paste(unique(na.omit(x)), collapse = ", "))
  
  summary_df$gwas_snps <- NA_character_
  summary_df$gwas_locus <- NA_character_
  summary_df$gwas_snps[as.numeric(names(snp_list))] <- snp_list
  summary_df$gwas_locus[as.numeric(names(locus_vec))] <- locus_vec
  summary_df
}

# ---------------------------
# Package-style high-level functions
# ---------------------------
load_clear_inputs <- function(se_path, trait_path) {
  se <- readRDS(se_path)
  snps <- readRDS(trait_path)
  se <- prepare_clear_se(se)
  list(se = se, snps = snps)
}

compute_gwas_overlap <- function(se, snps) {
  hits_se <- findOverlaps(snps, rowRanges(se))
  gwas_idx <- unique(subjectHits(hits_se))
  
  if (!length(gwas_idx)) {
    se_gwas <- se[0, , drop = FALSE]
  } else {
    se_gwas <- se[gwas_idx, , drop = FALSE]
    
    snp_name_col <- if ("names" %in% colnames(mcols(snps))) mcols(snps)$names else names(snps)
    signal_col <- if ("signal" %in% colnames(mcols(snps))) mcols(snps)$signal else rep(NA_character_, length(snps))
    
    ann <- data.frame(
      peak_idx = subjectHits(hits_se),
      CCVs = snp_name_col[queryHits(hits_se)],
      signal = signal_col[queryHits(hits_se)],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::group_by(peak_idx) %>%
      dplyr::summarise(
        CCVs = paste(unique(na.omit(CCVs)), collapse = ","),
        signal = paste(unique(na.omit(signal)), collapse = ","),
        .groups = "drop"
      )
    
    rd <- as.data.frame(rowData(se_gwas))
    if (!"CCVs" %in% names(rd)) rd$CCVs <- ""
    if (!"signal" %in% names(rd)) rd$signal <- ""
    
    pos <- match(ann$peak_idx, gwas_idx)
    rd$CCVs[pos] <- ann$CCVs
    rd$signal[pos] <- ann$signal
    rowData(se_gwas) <- S4Vectors::DataFrame(rd, row.names = rownames(se_gwas))
  }
  
  list(
    hits_se = hits_se,
    gwas_idx = gwas_idx,
    se_gwas = se_gwas,
    gwas_mat = assay(se_gwas, "raw_l2")
  )
}

summarize_inputs <- function(se_name, trait_name, se_path, trait_path, se, snps, overlap, results_dir) {
  n_unique_signals <- if ("signal" %in% colnames(mcols(snps))) length(unique(na.omit(mcols(snps)$signal))) else NA_integer_
  if (length(overlap$gwas_idx) > 0 && "signal" %in% colnames(mcols(snps))) {
    overlapping_signals <- mcols(snps)$signal[queryHits(overlap$hits_se)]
    n_unique_signals_in_peaks <- length(unique(na.omit(overlapping_signals)))
  } else {
    n_unique_signals_in_peaks <- 0
  }
  
  raw_mat <- assay(se, "raw")
  raw_col_sums <- colSums(raw_mat)
  
  raw_colsum_summary <- data.frame(
    celltype = names(raw_col_sums),
    total_raw_signal = as.numeric(raw_col_sums),
    stringsAsFactors = FALSE
  )
  raw_colsum_summary <- raw_colsum_summary[order(-raw_colsum_summary$total_raw_signal), ]
  write.table(raw_colsum_summary, file.path(results_dir, "raw_signal_per_celltype.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  input_summary <- data.frame(
    se_name = se_name,
    trait_name = trait_name,
    GWAS_trait = basename(trait_path),
    snATAC_dataset = basename(se_path),
    n_peaks = nrow(se),
    n_celltypes = ncol(se),
    n_snps = length(snps),
    n_unique_signals = n_unique_signals,
    n_gwas_peaks = nrow(overlap$gwas_mat),
    n_unique_signals_in_peaks = n_unique_signals_in_peaks,
    mean_total_raw_signal = mean(raw_col_sums),
    var_total_raw_signal = var(as.numeric(raw_col_sums)),
    max_signal_celltype = names(which.max(raw_col_sums)),
    max_signal_value = max(raw_col_sums),
    min_signal_celltype = names(which.min(raw_col_sums)),
    min_signal_value = min(raw_col_sums),
    stringsAsFactors = FALSE
  )
  write.table(input_summary, file.path(results_dir, "input_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  input_summary
}

plot_specificity_bar <- function(specificity_summary_l2, se_name, trait_name, plots_dir) {
  if (!nrow(specificity_summary_l2)) return(invisible(NULL))
  
  count_df <- as.data.frame(table(specificity_summary_l2$n_high), stringsAsFactors = FALSE)
  colnames(count_df) <- c("n_high", "count")
  count_df$proportion <- count_df$count / sum(count_df$count)
  lvls <- as.character(sort(as.integer(count_df$n_high)))
  if ("0" %in% lvls) lvls <- c(setdiff(lvls, "0"), "0")
  count_df$n_high <- factor(count_df$n_high, levels = lvls)
  count_df$bar_label <- paste0(count_df$count, " (", round(count_df$proportion * 100, 1), "%)")
  
  p_spec_l2 <- ggplot(count_df, aes(x = "GWAS peaks", y = proportion, fill = n_high)) +
    geom_bar(stat = "identity", width = 0.6, position = "stack") +
    geom_text(aes(label = bar_label), position = position_stack(vjust = 0.5), size = 3) +
    theme_bw() +
    labs(
      title = paste0(se_name, " x ", trait_name, " - High-specific cell counts (L2)"),
      x = "", y = "Proportion", fill = "# Highly specific\ncell types"
    ) +
    coord_flip()
  
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_l2_specificity_bar.pdf")), p_spec_l2, width = 8, height = 4)
  invisible(p_spec_l2)
}

addSpecificity <- function(gwas_mat, snps, high_thresh, mid_thresh, results_dir,
                           cell_lineage = NULL, se_gwas = NULL) {
  if (!nrow(gwas_mat)) return(data.frame())
  
  specificity_summary_l2_detailed <- compute_specificity_summary(gwas_mat)

  # --- Cell-type level: high / mid / low counts and names ---
  specificity_basic <- data.frame(
    peak = rownames(gwas_mat),
    n_high = rowSums(gwas_mat > high_thresh),
    n_mid = rowSums(gwas_mat > mid_thresh & gwas_mat <= high_thresh),
    n_low = rowSums(gwas_mat <= mid_thresh)
  )
  
  high_cells <- apply(gwas_mat, 1, function(x) names(x)[x > high_thresh])
  specificity_basic$high_cells <- sapply(high_cells, paste, collapse = ", ")
  specificity_basic$multiple_high <- lengths(high_cells) > 1

  mid_cells <- apply(gwas_mat, 1, function(x) {
    names(x)[x > mid_thresh & x <= high_thresh]
  })
  specificity_basic$mid_cells <- sapply(mid_cells, paste, collapse = ", ")
  
  specificity_summary_l2 <- specificity_summary_l2_detailed %>%
    left_join(specificity_basic, by = "peak")

  # --- Lineage-level dominance (ratio-based) ---
  if (!is.null(cell_lineage)) {
    dom_lineage <- compute_dominance_lineage(gwas_mat, cell_lineage, ratio_thresh = DEFAULT_RATIO_THRESH) %>%
      dplyr::select(
        peak,
        dominant_lineage = maximum_lineage,
        dominant_lineage_ratio = ratio,
        dominant_lineage_is_dominant = is_dominant
      )
    specificity_summary_l2 <- specificity_summary_l2 %>%
      left_join(dom_lineage, by = "peak")

    # --- Lineage-level high / mid counts and names ---
    lineage_vec <- cell_lineage[colnames(gwas_mat)]
    lineage_levels <- unique(as.character(lineage_vec))
    mat_sq <- gwas_mat^2
    lineage_mat <- sapply(lineage_levels, function(lg) {
      rowSums(mat_sq[, lineage_vec == lg, drop = FALSE])
    })
    lineage_mat <- as.matrix(lineage_mat)
    rownames(lineage_mat) <- rownames(gwas_mat)

    lin_high_thresh <- high_thresh^2
    lin_mid_thresh  <- mid_thresh^2

    specificity_summary_l2$n_high_lineage <- rowSums(lineage_mat > lin_high_thresh)
    specificity_summary_l2$n_mid_lineage  <- rowSums(lineage_mat > lin_mid_thresh & lineage_mat <= lin_high_thresh)
    specificity_summary_l2$n_low_lineage  <- rowSums(lineage_mat <= lin_mid_thresh)

    high_lin <- apply(lineage_mat, 1, function(x) names(x)[x > lin_high_thresh])
    mid_lin  <- apply(lineage_mat, 1, function(x) names(x)[x > lin_mid_thresh & x <= lin_high_thresh])
    specificity_summary_l2$high_lineages <- sapply(high_lin, paste, collapse = ", ")
    specificity_summary_l2$mid_lineages  <- sapply(mid_lin,  paste, collapse = ", ")

    specificity_summary_l2$top_lineage <- colnames(lineage_mat)[max.col(lineage_mat, ties.method = "first")]
  } else {
    specificity_summary_l2$dominant_lineage <- NA_character_
    specificity_summary_l2$dominant_lineage_ratio <- NA_real_
    specificity_summary_l2$dominant_lineage_is_dominant <- NA
    specificity_summary_l2$n_high_lineage <- NA_integer_
    specificity_summary_l2$n_mid_lineage  <- NA_integer_
    specificity_summary_l2$n_low_lineage  <- NA_integer_
    specificity_summary_l2$high_lineages  <- NA_character_
    specificity_summary_l2$mid_lineages   <- NA_character_
    specificity_summary_l2$top_lineage    <- NA_character_
  }

  # --- Top cell type (argmax L2, always populated) ---
  specificity_summary_l2$top_cell <- colnames(gwas_mat)[max.col(gwas_mat, ties.method = "first")]

  # --- Weight: mean raw accessibility across cell types ---
  if (!is.null(se_gwas) && "raw" %in% assayNames(se_gwas)) {
    raw_gwas <- assay(se_gwas, "raw")
    if (inherits(raw_gwas, "Matrix")) raw_gwas <- as.matrix(raw_gwas)
    specificity_summary_l2$weight <- rowMeans(raw_gwas[rownames(gwas_mat), , drop = FALSE])
  } else {
    specificity_summary_l2$weight <- NA_real_
  }

  specificity_summary_l2 <- add_gwas_peak_annotation(specificity_summary_l2, snps)
  specificity_summary_l2$locus <- ifelse(
    is.na(specificity_summary_l2$gwas_locus),
    "",
    as.character(specificity_summary_l2$gwas_locus)
  )
  
  specificity_summary_l2 <- specificity_summary_l2 %>%
    select(
      locus, peak, gwas_snps, gwas_locus,
      top_cell, maximum_cell, maximum_score, dominant_ratio, is_dominant,
      n_high, n_mid, n_low, high_cells, mid_cells, multiple_high,
      top_lineage, dominant_lineage, dominant_lineage_ratio, dominant_lineage_is_dominant,
      n_high_lineage, n_mid_lineage, n_low_lineage, high_lineages, mid_lineages,
      weight, k_cells, top_k_cells, coherence
    )
  
  write.table(specificity_summary_l2, file.path(results_dir, "specificity_summary_l2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  specificity_summary_l2
}

plot_dominance <- function(mat_l2, gwas_mat, se_name, trait_name, plots_dir, lineage_palette = NULL) {
  p_dom_l2 <- dominance_dodge_plot(mat_l2, gwas_mat, title_label = "L2mat", lineage_palette = lineage_palette)
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_dominance_l2.pdf")), p_dom_l2, width = 8, height = 6)
  invisible(p_dom_l2)
}

plot_dominance_lineage <- function(mat_l2, gwas_mat, cell_lineage, se_name, trait_name, plots_dir) {
  p_dom_lineage <- dominance_lineage_dodge_plot(mat_l2, gwas_mat, cell_lineage = cell_lineage, title_label = "lineage-summed L2^2")
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_dominance_lineage_l2.pdf")), p_dom_lineage, width = 8, height = 6)
  invisible(p_dom_lineage)
}

plot_domiance <- plot_dominance
plot_domiance_lineage <- plot_dominance_lineage

plot_coherence <- function(locus_results, se_name, trait_name, plots_dir) {
  if (!nrow(locus_results)) return(invisible(NULL))
  
  plot_df <- locus_results %>%
    dplyr::filter(n_peaks > 1) %>%
    dplyr::mutate(
      neg_log10_p = -log10(p_value),
      dominant_lineage = as.character(dominant_lineage),
      category = factor(
        category,
        levels = c("Highly_coherent_specific", "Coherent_moderate", "Intermediate", "Incoherent", "Single_peak")
      )
    ) %>%
    dplyr::filter(is.finite(neg_log10_p))
  
  if (!nrow(plot_df)) return(invisible(NULL))
  
  y_cap <- as.numeric(stats::quantile(plot_df$neg_log10_p, probs = 0.99, na.rm = TRUE))
  y_cap <- max(y_cap, -log10(0.05) + 0.2)
  
  p_locus_cat <- ggplot(plot_df, aes(x = concordance, y = neg_log10_p)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.55, color = "gray35") +
    geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.5, color = "gray45") +
    geom_point(aes(size = n_peaks, color = category, shape = dominant_lineage), alpha = 0.88) +
    scale_color_manual(
      values = c(
        Highly_coherent_specific = "#C43C39",
        Coherent_moderate = "#E7965A",
        Intermediate = "#5A5A5A",
        Incoherent = "#2F78B7",
        Single_peak = "#bdbdbd"
      ),
      drop = FALSE
    ) +
    scale_shape_discrete(name = "Dominant lineage") +
    scale_size_continuous(range = c(2.3, 8.2), breaks = c(2, 5, 10, 20, 30)) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.25), expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(ylim = c(0, y_cap)) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "gray30", size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "gray88", linewidth = 0.3),
      panel.border = element_rect(fill = NA, color = "gray40", linewidth = 0.45)
    ) +
    labs(
      x = "Mean pairwise rank concordance",
      y = "-log10(p-value)",
      title = paste0(se_name, " x ", trait_name, " - Locus concordance categories"),
      subtitle = paste0("Loci with >1 peak: n = ", nrow(plot_df), " (y-axis capped at 99th percentile)"),
      color = "Category",
      shape = "Dominant lineage",
      size = "N peaks"
    )
  
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_locus_concordance_category.pdf")), p_locus_cat, width = 8, height = 6)
  invisible(p_locus_cat)
}

plot_coherence_consensus <- function(locus_results_consensus, se_name, trait_name, plots_dir) {
  if (!nrow(locus_results_consensus)) return(invisible(NULL))

  plot_df <- locus_results_consensus %>%
    dplyr::filter(n_peaks > 1) %>%
    dplyr::mutate(
      neg_log10_p = -log10(p_value),
      dominant_lineage_consensus = as.character(dominant_lineage_consensus),
      category = factor(
        category,
        levels = c("Highly_coherent_specific", "Coherent_moderate", "Intermediate", "Incoherent", "Single_peak")
      )
    ) %>%
    dplyr::filter(is.finite(neg_log10_p))

  if (!nrow(plot_df)) return(invisible(NULL))

  y_cap <- as.numeric(stats::quantile(plot_df$neg_log10_p, probs = 0.99, na.rm = TRUE))
  y_cap <- max(y_cap, -log10(0.05) + 0.2)

  # Consensus fine vs p-value (category coloured)
  p_consensus_cat <- ggplot(plot_df, aes(x = consensus_fine, y = neg_log10_p)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.55, color = "gray35") +
    geom_point(aes(size = n_peaks, color = category, shape = dominant_lineage_consensus), alpha = 0.88) +
    scale_color_manual(
      values = c(
        Highly_coherent_specific = "#C43C39",
        Coherent_moderate = "#E7965A",
        Intermediate = "#5A5A5A",
        Incoherent = "#2F78B7",
        Single_peak = "#bdbdbd"
      ),
      drop = FALSE
    ) +
    scale_shape_discrete(name = "Dominant lineage") +
    scale_size_continuous(range = c(2.3, 8.2), breaks = c(2, 5, 10, 20, 30)) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "gray30", size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray40", linewidth = 0.45)
    ) +
    labs(
      x = "Weighted consensus (cell type)",
      y = "-log10(p-value)",
      title = paste0(se_name, " x ", trait_name, " - Consensus locus concordance"),
      subtitle = paste0("Loci with >1 peak: n = ", nrow(plot_df)),
      color = "Category", shape = "Dominant lineage", size = "N peaks"
    )
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_locus_concordance_consensus.pdf")),
         p_consensus_cat, width = 8, height = 6)

  # Lineage consensus vs p-value
  p_consensus_lineage <- ggplot(plot_df, aes(x = consensus_lineage, y = neg_log10_p)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_point(aes(size = n_peaks, color = category), alpha = 0.8) +
    scale_color_manual(
      values = c(
        Highly_coherent_specific = "#C43C39",
        Coherent_moderate = "#E7965A",
        Intermediate = "#5A5A5A",
        Incoherent = "#2F78B7",
        Single_peak = "#bdbdbd"
      ),
      drop = FALSE
    ) +
    theme_bw() + ggsci::scale_color_igv() +
    labs(
      x = "Weighted consensus (lineage)",
      y = "-log10(p-value)",
      title = paste0(se_name, " x ", trait_name, " - Lineage consensus concordance")
    )
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_locus_concordance_consensus_lineage.pdf")),
         p_consensus_lineage, width = 8, height = 6)

  # Z-score vs n_peaks
  p_consensus_z <- ggplot(plot_df, aes(y = z_score, x = n_peaks)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(aes(color = dominant_lineage_consensus, shape = category), alpha = 0.8) +
    theme_bw() + ggsci::scale_color_igv() +
    labs(
      x = "Number of GWAS peaks",
      y = "Z score (consensus)",
      title = paste0(se_name, " x ", trait_name, " - Consensus Z-scores")
    )
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_locus_zscore_consensus.pdf")),
         p_consensus_z, width = 8, height = 6)

  invisible(list(category = p_consensus_cat, lineage = p_consensus_lineage, zscore = p_consensus_z))
}


addGlobalEnrichment <- function(se, snps, mat_l2, results_dir, plots_dir, se_name, trait_name) {
  nmat_r <- apply(mat_l2, 2, function(x) rank(x, ties.method = "min"))
  nse_rank <- SummarizedExperiment(nmat_r, rowRanges = rowRanges(se))
  gwas_rank_idx <- unique(subjectHits(findOverlaps(snps, rowRanges(nse_rank))))
  gwas_rank_mat <- if (length(gwas_rank_idx)) assay(nse_rank[gwas_rank_idx, , drop = FALSE], 1) else nmat_r[0, , drop = FALSE]
  
  if (!nrow(gwas_rank_mat)) return(data.frame())
  
  N <- nrow(nse_rank)
  n <- nrow(gwas_rank_mat)
  exp_mean <- (N + 1) / 2
  se_z <- sqrt((N^2 - 1) / (12 * n))
  
  global_enrichment <- data.frame(observed_mean_rank = colMeans(gwas_rank_mat), stringsAsFactors = FALSE)
  global_enrichment$n_peaks <- n
  global_enrichment$exp_rank <- exp_mean
  global_enrichment$obs.exp <- global_enrichment$observed_mean_rank / exp_mean
  global_enrichment$z <- (global_enrichment$observed_mean_rank - exp_mean) / se_z
  global_enrichment$p <- pnorm(global_enrichment$z, lower.tail = FALSE)
  global_enrichment$fdr <- p.adjust(global_enrichment$p, method = "BH")
  global_enrichment$positive <- global_enrichment$z > 0 & global_enrichment$fdr < 0.05
  global_enrichment$celltype <- rownames(global_enrichment)
  
  write.table(global_enrichment, file.path(results_dir, "global_enrichment.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  p_enrich <- ggplot(global_enrichment, aes(x = reorder(celltype, z), y = z, fill = positive)) +
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    labs(
      x = "Cell type", y = "Z-score (rank enrichment)",
      title = paste0(se_name, " x ", trait_name, " - Global GWAS Enrichment"),
      fill = "FDR < 0.05"
    )
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_global_enrichment.pdf")), p_enrich, width = 8, height = 6)
  
  p_fdr <- ggplot(global_enrichment, aes(x = reorder(celltype, -log10(fdr)), y = -log10(fdr), fill = positive)) +
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    labs(
      x = "Cell type", y = "-log10(FDR)",
      title = paste0(se_name, " x ", trait_name, " - Global Enrichment (FDR)"),
      fill = "FDR < 0.05"
    )
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_global_enrichment_fdr.pdf")), p_fdr, width = 8, height = 6)
  
  global_enrichment
}

# ---------------------------
# Main (step-by-step)
# ---------------------------
processed_se_path <- file.path(processed_dir, paste0(se_name, "_processed_se.rds"))

message("\n--- Step 1: Load trait and processed/raw SE ---")
snps <- readRDS(trait_path)
message("  SNPs: ", length(snps))

if (file.exists(processed_se_path)) {
  message("  Found processed SE cache: ", processed_se_path)
  se <- readRDS(processed_se_path)
} else {
  message("  No processed SE cache found; loading and preparing raw SE")
  se <- readRDS(se_path)
  se <- prepare_clear_se(se)
  se <- compute_metrics_se_no_cov(se)
  saveRDS(se, processed_se_path)
  message("  Saved processed SE cache: ", processed_se_path)
}

required_assays <- c("raw", "raw_l2", "raw_l2_rank")
if (!all(required_assays %in% assayNames(se))) {
  message("  Cached SE missing expected assays; reloading raw SE and refreshing cache")
  se <- readRDS(se_path)
  se <- prepare_clear_se(se)
  se <- compute_metrics_se_no_cov(se)
  saveRDS(se, processed_se_path)
}

message("  SE:   ", nrow(se), " peaks x ", ncol(se), " cell types")

message("  Building lineage map from genome-wide raw_l2 and plotting dendrogram")
lineage_obj <- build_lineage_map(assay(se, "raw_l2"), k_lineage = k_lineage)
plot_lineage_dendrogram(lineage_obj, se_name, plots_dir)

message("\n--- Step 2: Extract L2 matrix ---")
mat_l2 <- assay(se, "raw_l2")

message("\n--- Step 3: GWAS overlap ---")
overlap <- compute_gwas_overlap(se, snps)
se_gwas <- overlap$se_gwas
gwas_mat <- overlap$gwas_mat
message("  GWAS peaks: ", nrow(gwas_mat))
if (nrow(gwas_mat) > 0) {
  n_sig_overlap <- length(unique(na.omit(unlist(strsplit(as.character(rowData(se_gwas)$signal), ",")))))
  gwas_se_name <- paste0(se_name, "_", trait_name, "_", nrow(gwas_mat), "peaks_", n_sig_overlap, "signals.rds")
  saveRDS(se_gwas, file.path(results_dir, gwas_se_name))
}

message("\n--- Step 4: Correlation heatmaps and densities ---")
plot_assay_cor_heatmaps(se, se_name, plots_dir)
plot_assay_densities(se, se_name, plots_dir, lineage_palette = lineage_obj$palette)
if (nrow(gwas_mat) > 1) {
  plot_assay_cor_heatmaps(se_gwas, paste0(se_name, "_", trait_name, "_GWAS"), plots_dir)
  plot_assay_densities(se_gwas, paste0(se_name, "_", trait_name, "_GWAS"), plots_dir, lineage_palette = lineage_obj$palette)
}

message("\n--- Step 5: Input summary ---")
summarize_inputs(se_name, trait_name, se_path, trait_path, se, snps, overlap, results_dir)

message("\n--- Step 6: L2 specificity summary ---")
high_thresh <- 1 / sqrt(3)
mid_thresh <- 1 / sqrt(ncol(mat_l2))
specificity_summary_l2 <- addSpecificity(gwas_mat, snps, high_thresh, mid_thresh, results_dir, cell_lineage = lineage_obj$cell_lineage)
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

message("\n--- Step 7: Dominance plots ---")
plot_dominance(mat_l2, gwas_mat, se_name, trait_name, plots_dir, lineage_palette = lineage_obj$palette)
plot_domiance_lineage(mat_l2, gwas_mat, lineage_obj$cell_lineage, se_name, trait_name, plots_dir)

message("\n--- Step 8: Locus-level permutation concordance ---")
locus_results <- addLocusCoherence(se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir)
if (nrow(locus_results)) {
  plot_coherence(locus_results, se_name, trait_name, plots_dir)
}

message("\n--- Step 8b: Locus-level consensus concordance ---")
locus_results_consensus <- addLocusCoherence_consensus(
  se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir,
  specificity_summary_l2 = specificity_summary_l2
)
if (nrow(locus_results_consensus)) {
  plot_coherence_consensus(locus_results_consensus, se_name, trait_name, plots_dir)
}

message("\n--- Step 9: Global rank-based enrichment (L2) ---")
addGlobalEnrichment(se, snps, mat_l2, results_dir, plots_dir, se_name, trait_name)

message("\n========================================")
message("CLEAR analysis complete: ", se_name, " x ", trait_name)
message("Results: ", normalizePath(results_dir))
message("Plots:   ", normalizePath(plots_dir))
message("========================================\n")
