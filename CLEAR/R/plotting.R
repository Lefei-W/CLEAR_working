# ============================================================
# plotting.R
# All visualisation functions for CLEAR pipeline
# ============================================================

.pvalStars <- function(p) {
  dplyr::case_when(
    is.na(p)    ~ "",
    p < 0.001   ~ "***",
    p < 0.01    ~ "**",
    p < 0.05    ~ "*",
    TRUE        ~ ""
  )
}

# cor_raw and cor_raw_l2 are kept in se_preparation as optional assays
# (corpcor-weighted versions), but not used in the current pipeline.
# Only raw and raw_l2 are plotted for QC.
PLOT_ASSAYS <- c("raw", "raw_l2")

plotAssayCorHeatmaps <- function(se, tissue, plot_dir = "plots") {
  out <- file.path(plot_dir, tissue)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  avail <- intersect(PLOT_ASSAYS, assayNames(se))
  
  for (a in avail) {
    mat <- assay(se, a)
    if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
    if (ncol(mat) < 2) next
    cor_mat <- cor(mat, use = "pairwise.complete.obs")
    rng <- range(cor_mat, na.rm = TRUE)
    # For L2-based correlations, fix the colour scale at [-1, 0, 1] so
    # the cell-type heatmap is directly comparable to the lineage
    # heatmap (and across tissues). Raw-signal correlations sit in a
    # narrow positive band where a [-1, 1] scale washes out structure,
    # so they keep a dynamic range.
    if (a == "raw_l2") {
      col_fun <- circlize::colorRamp2(c(-1, 0, 1),
                                      c("#2166AC", "white", "#B2182B"))
    } else {
      mid <- (rng[1] + rng[2]) / 2
      col_fun <- circlize::colorRamp2(c(rng[1], mid, rng[2]),
                                      c("#2166AC", "white", "#B2182B"))
    }
    fname <- file.path(out, paste0(tissue, "_", a, "_cor_heatmap.pdf"))
    pdf(fname, width = 8, height = 7)
    ht <- ComplexHeatmap::Heatmap(
      cor_mat,
      name = "Correlation",
      column_title = paste0(tissue, " - ", a, " correlation",
                            "  [range: ", round(rng[1], 2), " to ", round(rng[2], 2), "]"),
      clustering_method_rows = "ward.D2",
      clustering_method_columns = "ward.D2",
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      row_names_gp = grid::gpar(fontsize = 7),
      column_names_gp = grid::gpar(fontsize = 7),
      col = col_fun
    )
    ComplexHeatmap::draw(ht)
    dev.off()
  }
  invisible(out)
}

plotAssayDensities <- function(se, tissue, plot_dir = "plots", lineage_palette = NULL) {
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

plotLineageDendrogram <- function(lineage_obj, se_name, plot_dir = "plots",
                                  k_preview = 3:6) {
  fname <- file.path(plot_dir, paste0(se_name, "_celltype_lineage_dendrogram.pdf"))
  hc <- lineage_obj$hc
  k_lineage <- lineage_obj$k_lineage
  n_ct <- length(hc$order)

  # Restrict preview ks to a sensible range
  k_preview <- sort(unique(as.integer(k_preview)))
  k_preview <- k_preview[k_preview >= 2L & k_preview < n_ct]

  # Cut height for k clusters = midpoint between the (n - k)th and
  # (n - k + 1)th merge heights (sorted ascending). rect.hclust uses
  # the same convention via cutree(k=).
  heights_asc <- sort(hc$height)
  cut_height <- function(k) {
    idx <- length(heights_asc) - k + 1L
    if (idx < 1L || idx > length(heights_asc)) return(NA_real_)
    if (idx == length(heights_asc)) return(heights_asc[idx] * 1.05)
    mean(c(heights_asc[idx], heights_asc[idx + 1L]))
  }
  preview_h <- setNames(vapply(k_preview, cut_height, numeric(1)), k_preview)

  # Colour palette for preview lines (distinct from the highlight colour)
  pal <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a", "#e31a1c", "#b15928")
  line_cols <- setNames(pal[seq_along(k_preview)], as.character(k_preview))

  grDevices::pdf(fname, width = 11, height = 6.5)
  op <- par(mar = c(10, 4, 4, 8), xpd = NA)
  plot(
    hc,
    main = paste0(se_name, " - Cell-type lineage dendrogram (chosen k=", k_lineage, ")"),
    xlab = "",
    sub  = "Distance = 1 - Pearson correlation on genome-wide raw_l2",
    cex  = 0.7
  )

  # Preview cut heights for several k values
  for (kk in as.character(k_preview)) {
    h <- preview_h[[kk]]
    if (is.finite(h)) {
      abline(h = h, col = line_cols[[kk]], lty = 2, lwd = 1.2)
    }
  }

  # Highlight the chosen k with rectangles
  rect.hclust(hc, k = k_lineage, border = "firebrick")

  # Legend: preview ks + chosen
  legend("topright",
         inset  = c(-0.13, 0),
         legend = c(paste0("k = ", k_preview, " (h=", signif(preview_h, 3), ")"),
                    paste0("chosen k = ", k_lineage)),
         col    = c(line_cols, "firebrick"),
         lty    = c(rep(2, length(k_preview)), 1),
         lwd    = c(rep(1.2, length(k_preview)), 2),
         bty    = "n",
         cex    = 0.75)
  par(op)
  grDevices::dev.off()

  invisible(fname)
}

plotCelltypeStackbars <- function(
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
  # Each of the top-k contributors must exceed the uniform-vector reference
  # 1/sqrt(n); checking the k-th (smallest of top-k, since sorted desc) is
  # sufficient to guarantee this for all top-k cell types.
  min_l2 <- 1 / sqrt(n_celltypes)
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

# ---------------------------------------------------------------
# Plot: Cell-type stacked bar chart — ALL variant-peaks
# Same logic as plotCelltypeStackbars but does not filter peaks.
# Ordering:
#   1. Dominant (passing) peaks first, grouped by dominant cell type.
#   2. Cell-type groups ordered by number of passing peaks they dominate (desc).
#   3. Within each cell-type group, peaks ordered by decreasing max L2^2.
#   4. Non-passing peaks follow, same cell-type priority order.
# A dotted line marks the pass/fail boundary.
# ---------------------------------------------------------------
plotCelltypeStackbarsAll <- function(
    se_overlap,
    value_metric    = "raw_l2",
    cumvar_threshold = 0.8,
    max_k           = 3,
    outfile,
    lineage_palette = NULL) {
  stopifnot(!missing(outfile))

  if (!value_metric %in% assayNames(se_overlap)) {
    warning("Assay '", value_metric, "' not found; skipping all-peak stackbar.")
    return(invisible(FALSE))
  }

  val_mat <- SummarizedExperiment::assay(se_overlap, value_metric)
  if (inherits(val_mat, "Matrix")) val_mat <- as.matrix(val_mat)
  storage.mode(val_mat) <- "double"

  n_celltypes <- ncol(val_mat)
  val_sq      <- val_mat^2

  # Dominance: peak is "dominant" iff max1(L2) > ratio_thresh * max2(L2)
  # using the canonical pipeline definition (addDominance).
  ratio_thresh <- 2
  dom_df  <- addDominance(val_mat, ratio_thresh = ratio_thresh)
  pass_v  <- as.logical(dom_df$is_dominant)
  dom_ct_v <- as.character(dom_df$maximum_ct)
  max_sq_v <- vapply(seq_len(nrow(val_sq)), function(i) max(val_sq[i, ]),
                     numeric(1))

  # Peak labels
  rd         <- as.data.frame(rowData(se_overlap))
  peakid     <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se_overlap)
  ccv_label  <- if ("CCVs"   %in% names(rd)) rd$CCVs   else ""
  sig_label  <- if ("signal" %in% names(rd)) rd$signal else ""
  base_label <- ifelse(nzchar(ccv_label), ccv_label, peakid)
  labels     <- ifelse(nzchar(sig_label), paste0(base_label, " (", sig_label, ")"),
                       base_label)

  meta_df <- tibble::tibble(
    idx        = seq_len(nrow(val_mat)),
    pass       = pass_v,
    dom_ct     = dom_ct_v,
    max_sq     = max_sq_v,
    peak_label = labels
  )

  # Ordering: dominant (passing) peaks first, grouped by dominant cell
  # type. Cell-type group priority = number of passing peaks dominated by
  # that cell type (desc). Within a cell type, peaks ordered by decreasing
  # max L2^2. Non-passing peaks follow, in the SAME cell-type priority
  # order (cell types not represented in passing peaks appended last,
  # ordered by their non-passing count desc).
  pass_meta <- meta_df %>% dplyr::filter(pass)
  fail_meta <- meta_df %>% dplyr::filter(!pass)

  pass_ct_counts <- sort(table(pass_meta$dom_ct), decreasing = TRUE)
  fail_ct_counts <- sort(table(fail_meta$dom_ct), decreasing = TRUE)
  ct_priority    <- c(names(pass_ct_counts),
                      setdiff(names(fail_ct_counts), names(pass_ct_counts)))

  pass_df <- pass_meta %>%
    dplyr::mutate(ct_rank = match(dom_ct, ct_priority)) %>%
    dplyr::arrange(ct_rank, -max_sq)
  fail_df <- fail_meta %>%
    dplyr::mutate(ct_rank = match(dom_ct, ct_priority)) %>%
    dplyr::arrange(ct_rank, -max_sq)

  peak_order <- c(pass_df$peak_label, fail_df$peak_label)

  # Long table covering ALL peaks
  long_list <- lapply(seq_len(nrow(val_mat)), function(i) {
    row_sq <- val_sq[i, ]
    ct_ord <- order(row_sq, decreasing = TRUE)
    data.frame(
      peak_label = labels[i],
      celltype   = colnames(val_sq)[ct_ord],
      value      = row_sq[ct_ord],
      stack_rank = seq_along(ct_ord),
      pass       = meta_df$pass[i],
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, long_list)
  long$peak_label <- factor(long$peak_label, levels = rev(peak_order))

  ct_totals <- tapply(long$value, long$celltype, sum, na.rm = TRUE)
  ct_levels <- names(sort(ct_totals, decreasing = TRUE))
  long$celltype <- factor(long$celltype, levels = ct_levels)
  long <- long %>% dplyr::arrange(peak_label, -value)

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

  threshold_pct <- cumvar_threshold * 100
  n_total <- nrow(meta_df)
  n_pass  <- sum(meta_df$pass)

  # Boundary line between passing and non-passing peaks (in flipped coords).
  # peak_label factor levels are rev(peak_order): failing peaks first, then
  # passing peaks at the top. Boundary sits between them.
  boundary_y <- if (n_pass > 0 && n_pass < n_total) nrow(fail_df) + 0.5 else NA_real_

  p <- ggplot(long, aes(x = peak_label, y = value, fill = celltype,
                        group = factor(stack_rank))) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.85) +
    geom_hline(yintercept = cumvar_threshold,
               linetype = "dashed", color = "black", linewidth = 0.3) +
    coord_flip() +
    fill_scale +
    labs(
      title = paste0(
        "All variant-peaks: cell-type L2^2 decomposition (", n_total, " peaks)"
      ),
      subtitle = paste0(
        n_pass, " peaks dominant (max L2 > ", ratio_thresh,
        "x second; canonical addDominance rule)",
        " (top of plot, grouped by dominant cell type by # dominant peaks); ",
        "remaining ",
        n_total - n_pass, " peaks below"
      ),
      x = "Peak / CCV",
      y = "Variance explained (L2^2, cumulative = 1.0)",
      fill = "Cell type"
    ) +
    theme_bw() +
    theme(
      axis.text.y  = element_text(size = 3),
      axis.text.x  = element_text(size = 6),
      legend.text  = element_text(size = 6),
      legend.title = element_text(size = 8),
      plot.title   = element_text(size = 10)
    )

  if (!is.na(boundary_y)) {
    p <- p + geom_vline(xintercept = boundary_y,
                        linetype = "dotted", color = "grey30", linewidth = 0.4)
  }

  pdf_height <- max(6, n_total * 0.10 + 3)
  grDevices::pdf(outfile, width = 12, height = min(pdf_height, 80))
  print(p)
  grDevices::dev.off()

  message("  Saved all-peak cell-type stackbar (", n_total,
          " peaks; ", n_pass, " pass): ", outfile)
  # Return the peak ordering (top-to-bottom) so callers can reuse it
  # for the lineage-level all-peak plot, ensuring identical row order
  # for side-by-side visual comparison.
  invisible(peak_order)
}

# ---------------------------------------------------------------
# Plot: Lineage-level stacked bar chart
# Filters peaks where top <=max_k lineages explain >=cumvar_threshold
# of L2^2 variance.  Bars ordered high-to-low within each peak.
# ---------------------------------------------------------------
plotLineageStackbars <- function(
    se_overlap,
    value_metric = "raw_l2",
    cumvar_threshold = 0.8,
    max_k = 2,
    outfile,
    cell_lineage,
    lineage_colours = NULL,
    max_peaks = 200) {
  stopifnot(!missing(outfile))

  if (!value_metric %in% assayNames(se_overlap)) {
    warning("Assay '", value_metric, "' not found; skipping lineage stackbar.")
    return(invisible(FALSE))
  }

  val_mat <- SummarizedExperiment::assay(se_overlap, value_metric)
  if (inherits(val_mat, "Matrix")) val_mat <- as.matrix(val_mat)
  storage.mode(val_mat) <- "double"

  # Square the lineage scores: lineage_mat contains sqrt(rowSums(L2²)) per group,
  # so squaring gives the fraction of L2² variance per lineage — these sum to 1
  # per peak and are the correct values to stack and filter against cumvar_threshold.
  lineage_mat_sqrt <- getLineageMatrix(val_mat, cell_lineage)
  if (ncol(lineage_mat_sqrt) < 2) {
    warning("Fewer than 2 lineages; skipping lineage stackbar.")
    return(invisible(FALSE))
  }
  lineage_mat <- lineage_mat_sqrt^2  # L2² per lineage: rows sum to 1

  # Filter: top <=max_k lineages must explain >=cumvar_threshold of L2² variance
  peak_info <- lapply(seq_len(nrow(lineage_mat)), function(i) {
    row_vals <- lineage_mat[i, ]
    ord <- order(row_vals, decreasing = TRUE)
    cum_vals <- cumsum(row_vals[ord])
    k <- which(cum_vals >= cumvar_threshold)[1]
    if (is.na(k) || k > max_k) return(NULL)
    list(idx = i, k = k, cum_var = cum_vals[k])
  })
  peak_info <- Filter(Negate(is.null), peak_info)

  if (!length(peak_info)) {
    threshold_pct <- cumvar_threshold * 100
    warning("No peaks where top <=", max_k, " lineages explain >=", threshold_pct, "% L2^2")
    grDevices::pdf(outfile, width = 10, height = 6)
    plot.new(); text(0.5, 0.5, "No peaks meeting criteria")
    dev.off()
    return(invisible(FALSE))
  }

  sel   <- sapply(peak_info, `[[`, "idx")
  sel_k <- sapply(peak_info, `[[`, "k")

  if (length(sel) > max_peaks) {
    ord <- order(sel_k, -apply(lineage_mat[sel, , drop = FALSE], 1, max))
    keep <- ord[seq_len(max_peaks)]
    sel <- sel[keep]; sel_k <- sel_k[keep]; peak_info <- peak_info[keep]
  }

  # Peak labels
  rd <- as.data.frame(rowData(se_overlap))
  peakid <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se_overlap)
  ccv_label <- if ("CCVs" %in% names(rd)) rd$CCVs else ""
  sig_label <- if ("signal" %in% names(rd)) rd$signal else ""
  base_label <- ifelse(nzchar(ccv_label), ccv_label, peakid)
  labels <- ifelse(nzchar(sig_label), paste0(base_label, " (", sig_label, ")"), base_label)

  dom_lin <- apply(lineage_mat[sel, , drop = FALSE], 1, function(x) names(which.max(x)))

  # Build long data with stack_rank (high-to-low within each peak)
  long_list <- lapply(seq_along(sel), function(j) {
    i <- sel[j]
    row_vals <- lineage_mat[i, ]
    lin_ord <- order(row_vals, decreasing = TRUE)
    data.frame(
      peak_label = labels[i],
      lineage    = colnames(lineage_mat)[lin_ord],
      value      = row_vals[lin_ord],
      stack_rank = seq_along(lin_ord),
      k          = sel_k[j],
      dom_lin    = dom_lin[j],
      max_val    = max(row_vals),
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, long_list)

  peak_order <- long %>%
    dplyr::distinct(peak_label, k, dom_lin, max_val) %>%
    dplyr::arrange(k, dom_lin, -max_val) %>%
    dplyr::pull(peak_label)
  long$peak_label <- factor(long$peak_label, levels = rev(peak_order))

  lin_totals <- tapply(long$value, long$lineage, sum, na.rm = TRUE)
  lin_levels <- names(sort(lin_totals, decreasing = TRUE))
  long$lineage <- factor(long$lineage, levels = lin_levels)
  long <- long %>% dplyr::arrange(peak_label, -value)

  if (!is.null(lineage_colours) && all(lin_levels %in% names(lineage_colours))) {
    fill_scale <- scale_fill_manual(values = lineage_colours)
  } else {
    fill_scale <- ggsci::scale_fill_npg()
  }

  threshold_pct <- cumvar_threshold * 100
  n_sel <- length(sel)
  p <- ggplot(long, aes(x = peak_label, y = value, fill = lineage, group = factor(stack_rank))) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.85) +
    geom_hline(yintercept = cumvar_threshold, linetype = "dashed", color = "black", linewidth = 0.3) +
    coord_flip() +
    fill_scale +
    labs(
      title = paste0("GWAS peaks: top <=", max_k, " lineages explain >=", threshold_pct, "% L2^2"),
      subtitle = paste0(n_sel, " peaks selected"),
      x = "Peak / CCV",
      y = "Lineage L2^2 (sums to 1.0)",
      fill = "Lineage"
    ) +
    theme_bw() +
    theme(
      axis.text.y  = element_text(size = 4),
      axis.text.x  = element_text(size = 6),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.title   = element_text(size = 10)
    )

  pdf_height <- max(6, n_sel * 0.12 + 3)
  grDevices::pdf(outfile, width = 12, height = min(pdf_height, 50))
  print(p)
  grDevices::dev.off()

  message("  Saved lineage stackbar plot (", n_sel, " peaks): ", outfile)
  invisible(TRUE)
}

# ---------------------------------------------------------------
# Plot: Lineage-level stacked bar chart — ALL variant-peaks
# Counterpart to plotCelltypeStackbarsAll at lineage resolution.
# Plots EVERY variant-peak using a LINEAGE-derived ordering (independent
# of the cell-type plot's ordering):
#   1. Dominant (passing) peaks first, grouped by dominant lineage.
#   2. Lineage groups ordered by number of passing peaks they dominate (desc).
#   3. Within each lineage group, peaks ordered by decreasing max lineage L2^2.
#   4. Non-passing peaks follow, same lineage priority order.
# A dotted line marks the pass/fail boundary.
# ---------------------------------------------------------------
plotLineageStackbarsAll <- function(
    se_overlap,
    value_metric    = "raw_l2",
    cumvar_threshold = 0.8,
    max_k           = 2,
    outfile,
    cell_lineage,
    lineage_colours = NULL) {
  stopifnot(!missing(outfile))

  if (!value_metric %in% assayNames(se_overlap)) {
    warning("Assay '", value_metric, "' not found; skipping all-peak lineage stackbar.")
    return(invisible(FALSE))
  }

  val_mat <- SummarizedExperiment::assay(se_overlap, value_metric)
  if (inherits(val_mat, "Matrix")) val_mat <- as.matrix(val_mat)
  storage.mode(val_mat) <- "double"

  lineage_mat_sqrt <- getLineageMatrix(val_mat, cell_lineage)
  if (ncol(lineage_mat_sqrt) < 2) {
    warning("Fewer than 2 lineages; skipping all-peak lineage stackbar.")
    return(invisible(FALSE))
  }
  lineage_mat <- lineage_mat_sqrt^2  # rows sum to 1

  n_lineages <- ncol(lineage_mat)

  # Dominance at the lineage level: max1(lineage L2) > ratio_thresh * max2.
  # addDominanceLineage builds an internal lineage matrix using the
  # same sqrt(rowSums(L2^2)) construction as getLineageMatrix, so
  # its dominance flags align with the values we're plotting here.
  ratio_thresh <- 2
  dom_df    <- addDominanceLineage(val_mat, cell_lineage,
                                         ratio_thresh = ratio_thresh)
  pass_v    <- as.logical(dom_df$is_dominant)
  dom_lin_v <- as.character(dom_df$maximum_lineage)
  max_val_v <- vapply(seq_len(nrow(lineage_mat)),
                      function(i) max(lineage_mat[i, ]), numeric(1))

  rd         <- as.data.frame(rowData(se_overlap))
  peakid     <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se_overlap)
  ccv_label  <- if ("CCVs"   %in% names(rd)) rd$CCVs   else ""
  sig_label  <- if ("signal" %in% names(rd)) rd$signal else ""
  base_label <- ifelse(nzchar(ccv_label), ccv_label, peakid)
  labels     <- ifelse(nzchar(sig_label), paste0(base_label, " (", sig_label, ")"),
                       base_label)

  meta_df <- tibble::tibble(
    idx        = seq_len(nrow(lineage_mat)),
    pass       = pass_v,
    dom_lin    = dom_lin_v,
    max_val    = max_val_v,
    peak_label = labels
  )

  # Ordering: dominant (passing) peaks first, grouped by dominant lineage.
  # Lineage group priority = number of passing peaks dominated by that
  # lineage (desc). Within a lineage, peaks ordered by decreasing max
  # lineage L2^2. Non-passing peaks follow with the SAME lineage priority
  # (lineages absent from passing block appended last by their non-pass
  # count desc).  This is computed independently from the cell-type plot
  # so the lineage ordering is lineage-derived end-to-end.
  pass_meta <- meta_df %>% dplyr::filter(pass)
  fail_meta <- meta_df %>% dplyr::filter(!pass)

  pass_lin_counts <- sort(table(pass_meta$dom_lin), decreasing = TRUE)
  fail_lin_counts <- sort(table(fail_meta$dom_lin), decreasing = TRUE)
  lin_priority    <- c(names(pass_lin_counts),
                       setdiff(names(fail_lin_counts), names(pass_lin_counts)))

  pass_df <- pass_meta %>%
    dplyr::mutate(lin_rank = match(dom_lin, lin_priority)) %>%
    dplyr::arrange(lin_rank, -max_val)
  fail_df <- fail_meta %>%
    dplyr::mutate(lin_rank = match(dom_lin, lin_priority)) %>%
    dplyr::arrange(lin_rank, -max_val)

  peak_order <- c(pass_df$peak_label, fail_df$peak_label)
  n_fail     <- nrow(fail_df)

  long_list <- lapply(seq_len(nrow(lineage_mat)), function(i) {
    row_vals <- lineage_mat[i, ]
    lin_ord  <- order(row_vals, decreasing = TRUE)
    data.frame(
      peak_label = labels[i],
      lineage    = colnames(lineage_mat)[lin_ord],
      value      = row_vals[lin_ord],
      stack_rank = seq_along(lin_ord),
      pass       = meta_df$pass[i],
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, long_list)
  long$peak_label <- factor(long$peak_label, levels = rev(peak_order))

  lin_totals <- tapply(long$value, long$lineage, sum, na.rm = TRUE)
  lin_levels <- names(sort(lin_totals, decreasing = TRUE))
  long$lineage <- factor(long$lineage, levels = lin_levels)
  long <- long %>% dplyr::arrange(peak_label, -value)

  if (!is.null(lineage_colours) && all(lin_levels %in% names(lineage_colours))) {
    fill_scale <- scale_fill_manual(values = lineage_colours)
  } else {
    fill_scale <- ggsci::scale_fill_npg()
  }

  threshold_pct <- cumvar_threshold * 100
  n_total <- nrow(meta_df)
  n_pass  <- sum(meta_df$pass)

  boundary_y <- if (n_pass > 0 && n_pass < n_total) n_fail + 0.5 else NA_real_

  subtitle_txt <- paste0(
    n_pass, " peaks dominant (max lineage L2 > ", ratio_thresh,
    "x second; canonical addDominanceLineage rule)",
    " (top of plot, grouped by dominant lineage by # dominant peaks); ",
    "remaining ", n_total - n_pass, " peaks below."
  )

  p <- ggplot(long, aes(x = peak_label, y = value, fill = lineage,
                        group = factor(stack_rank))) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.85) +
    geom_hline(yintercept = cumvar_threshold,
               linetype = "dashed", color = "black", linewidth = 0.3) +
    coord_flip() +
    fill_scale +
    labs(
      title = paste0(
        "All variant-peaks: lineage L2^2 decomposition (", n_total, " peaks)"
      ),
      subtitle = subtitle_txt,
      x = "Peak / CCV",
      y = "Lineage L2^2 (sums to 1.0)",
      fill = "Lineage"
    ) +
    theme_bw() +
    theme(
      axis.text.y  = element_text(size = 3),
      axis.text.x  = element_text(size = 6),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.title   = element_text(size = 10)
    )

  if (!is.na(boundary_y)) {
    p <- p + geom_vline(xintercept = boundary_y,
                        linetype = "dotted", color = "grey30", linewidth = 0.4)
  }

  pdf_height <- max(6, n_total * 0.10 + 3)
  grDevices::pdf(outfile, width = 12, height = min(pdf_height, 80))
  print(p)
  grDevices::dev.off()

  message("  Saved all-peak lineage stackbar (", n_total,
          " peaks; ", n_pass, " pass): ", outfile)
  invisible(peak_order)
}

# ---------------------------------------------------------------
# Plot: Cell-type stacked bar chart — grouped by LOCUS (signal)
# Variant-peaks at the same GWAS locus (rowData$signal) are kept
# together so multi-peak loci can be visually inspected for coherence
# across their constituent peaks.
#
# Ordering:
#   1. Loci ordered by n_peaks (desc); ties broken by max peak L2^2 (desc).
#   2. Within each locus, peaks ordered by decreasing max L2^2.
#   3. Optional multi_peak_only=TRUE drops single-peak loci (default FALSE).
# A faint dotted separator is drawn between consecutive loci so groups
# are visually distinct.
# ---------------------------------------------------------------
plotCelltypeStackbarsByLocus <- function(
    se_overlap,
    value_metric    = "raw_l2",
    cumvar_threshold = 0.8,
    outfile,
    multi_peak_only  = FALSE,
    lineage_palette  = NULL) {
  stopifnot(!missing(outfile))

  if (!value_metric %in% assayNames(se_overlap)) {
    warning("Assay '", value_metric, "' not found; skipping locus stackbar.")
    return(invisible(FALSE))
  }

  val_mat <- SummarizedExperiment::assay(se_overlap, value_metric)
  if (inherits(val_mat, "Matrix")) val_mat <- as.matrix(val_mat)
  storage.mode(val_mat) <- "double"
  val_sq <- val_mat^2

  rd <- as.data.frame(rowData(se_overlap))
  if (!"signal" %in% names(rd)) {
    warning("rowData missing 'signal' (locus id); skipping locus stackbar.")
    return(invisible(FALSE))
  }
  locus_id <- as.character(rd$signal)
  locus_id[is.na(locus_id) | !nzchar(locus_id)] <- "unassigned"

  peakid     <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se_overlap)
  ccv_label  <- if ("CCVs" %in% names(rd)) rd$CCVs else ""
  base_label <- ifelse(nzchar(ccv_label), ccv_label, peakid)
  # Peak label retains the CCV id; locus is encoded separately and used
  # for grouping, ordering and the visible group strip.
  labels <- base_label

  max_sq_v <- vapply(seq_len(nrow(val_sq)), function(i) max(val_sq[i, ]),
                     numeric(1))

  meta_df <- tibble::tibble(
    idx        = seq_len(nrow(val_mat)),
    locus      = locus_id,
    max_sq     = max_sq_v,
    peak_label = labels
  )

  # Locus-level summary: n_peaks per locus, max L2^2 across its peaks
  locus_summary <- meta_df %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(n_peaks = dplyr::n(),
                     locus_max = max(max_sq),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_peaks), dplyr::desc(locus_max))

  if (multi_peak_only) {
    locus_summary <- dplyr::filter(locus_summary, n_peaks >= 2)
    if (!nrow(locus_summary)) {
      warning("No multi-peak loci; skipping locus stackbar.")
      grDevices::pdf(outfile, width = 10, height = 6)
      plot.new(); text(0.5, 0.5, "No multi-peak loci")
      grDevices::dev.off()
      return(invisible(FALSE))
    }
    meta_df <- dplyr::filter(meta_df, locus %in% locus_summary$locus)
  }

  locus_levels <- locus_summary$locus
  meta_df <- meta_df %>%
    dplyr::mutate(locus_rank = match(locus, locus_levels)) %>%
    dplyr::arrange(locus_rank, -max_sq)

  # Disambiguate any duplicated peak_labels by appending the locus id
  # (rowData labels can repeat when the same CCV string is reused).
  dup_lab <- duplicated(meta_df$peak_label) |
             duplicated(meta_df$peak_label, fromLast = TRUE)
  meta_df$peak_label[dup_lab] <- paste0(meta_df$peak_label[dup_lab],
                                        " [", meta_df$locus[dup_lab], "]")

  peak_order <- meta_df$peak_label

  long_list <- lapply(seq_len(nrow(meta_df)), function(j) {
    i <- meta_df$idx[j]
    row_sq <- val_sq[i, ]
    ct_ord <- order(row_sq, decreasing = TRUE)
    data.frame(
      peak_label = meta_df$peak_label[j],
      locus      = meta_df$locus[j],
      celltype   = colnames(val_sq)[ct_ord],
      value      = row_sq[ct_ord],
      stack_rank = seq_along(ct_ord),
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, long_list)
  long$peak_label <- factor(long$peak_label, levels = rev(peak_order))

  ct_totals <- tapply(long$value, long$celltype, sum, na.rm = TRUE)
  ct_levels <- names(sort(ct_totals, decreasing = TRUE))
  long$celltype <- factor(long$celltype, levels = ct_levels)
  long <- long %>% dplyr::arrange(peak_label, -value)

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

  # Compute between-locus separator positions on the (flipped) x-axis.
  # peak_label factor levels are rev(peak_order), so the first locus
  # sits at the TOP of the plot. Separator between locus i and i+1 is
  # placed at cumulative peak-count from the bottom + 0.5.
  n_per_locus <- meta_df %>%
    dplyr::group_by(locus_rank) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(locus_rank) %>%
    dplyr::pull(n)
  n_total <- sum(n_per_locus)
  if (length(n_per_locus) > 1) {
    # cumulative count from the BOTTOM (last locus) upwards
    rev_counts <- rev(n_per_locus)
    cum_bottom <- cumsum(rev_counts)
    sep_positions <- head(cum_bottom, -1) + 0.5
  } else {
    sep_positions <- numeric(0)
  }

  n_loci    <- nrow(locus_summary)
  n_multi   <- sum(locus_summary$n_peaks >= 2)
  subtitle_txt <- paste0(
    n_total, " peaks across ", n_loci, " loci (",
    n_multi, " multi-peak); locus separators are dotted lines"
  )

  p <- ggplot(long, aes(x = peak_label, y = value, fill = celltype,
                        group = factor(stack_rank))) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.85) +
    geom_hline(yintercept = cumvar_threshold,
               linetype = "dashed", color = "black", linewidth = 0.3) +
    coord_flip() +
    fill_scale +
    labs(
      title    = "Variant-peaks grouped by locus: cell-type L2^2 decomposition",
      subtitle = subtitle_txt,
      x        = "Peak / CCV (grouped by locus)",
      y        = "Variance explained (L2^2, cumulative = 1.0)",
      fill     = "Cell type"
    ) +
    theme_bw() +
    theme(
      axis.text.y  = element_text(size = 4),
      axis.text.x  = element_text(size = 6),
      legend.text  = element_text(size = 6),
      legend.title = element_text(size = 8),
      plot.title   = element_text(size = 10)
    )

  if (length(sep_positions)) {
    p <- p + geom_vline(xintercept = sep_positions,
                        linetype = "dotted", color = "grey40",
                        linewidth = 0.35)
  }

  pdf_height <- max(6, n_total * 0.14 + 3)
  grDevices::pdf(outfile, width = 12, height = min(pdf_height, 60))
  print(p)
  grDevices::dev.off()

  message("  Saved locus-grouped cell-type stackbar (",
          n_total, " peaks, ", n_loci, " loci): ", outfile)
  invisible(peak_order)
}

# ---------------------------------------------------------------
# Plot: Lineage-level stacked bar chart — grouped by LOCUS (signal)
# Lineage counterpart of plotCelltypeStackbarsByLocus. Same
# locus-grouped ordering rules; values are lineage L2^2 (sums to 1
# per peak).
# ---------------------------------------------------------------
plotLineageStackbarsByLocus <- function(
    se_overlap,
    value_metric     = "raw_l2",
    cumvar_threshold = 0.8,
    outfile,
    cell_lineage,
    multi_peak_only  = FALSE,
    lineage_colours  = NULL) {
  stopifnot(!missing(outfile))

  if (!value_metric %in% assayNames(se_overlap)) {
    warning("Assay '", value_metric, "' not found; skipping locus lineage stackbar.")
    return(invisible(FALSE))
  }

  val_mat <- SummarizedExperiment::assay(se_overlap, value_metric)
  if (inherits(val_mat, "Matrix")) val_mat <- as.matrix(val_mat)
  storage.mode(val_mat) <- "double"

  lineage_mat_sqrt <- getLineageMatrix(val_mat, cell_lineage)
  if (ncol(lineage_mat_sqrt) < 2) {
    warning("Fewer than 2 lineages; skipping locus lineage stackbar.")
    return(invisible(FALSE))
  }
  lineage_mat <- lineage_mat_sqrt^2  # rows sum to 1

  rd <- as.data.frame(rowData(se_overlap))
  if (!"signal" %in% names(rd)) {
    warning("rowData missing 'signal' (locus id); skipping locus lineage stackbar.")
    return(invisible(FALSE))
  }
  locus_id <- as.character(rd$signal)
  locus_id[is.na(locus_id) | !nzchar(locus_id)] <- "unassigned"

  peakid     <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se_overlap)
  ccv_label  <- if ("CCVs" %in% names(rd)) rd$CCVs else ""
  base_label <- ifelse(nzchar(ccv_label), ccv_label, peakid)
  labels <- base_label

  max_val_v <- vapply(seq_len(nrow(lineage_mat)),
                      function(i) max(lineage_mat[i, ]), numeric(1))

  meta_df <- tibble::tibble(
    idx        = seq_len(nrow(lineage_mat)),
    locus      = locus_id,
    max_val    = max_val_v,
    peak_label = labels
  )

  locus_summary <- meta_df %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(n_peaks = dplyr::n(),
                     locus_max = max(max_val),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_peaks), dplyr::desc(locus_max))

  if (multi_peak_only) {
    locus_summary <- dplyr::filter(locus_summary, n_peaks >= 2)
    if (!nrow(locus_summary)) {
      warning("No multi-peak loci; skipping locus lineage stackbar.")
      grDevices::pdf(outfile, width = 10, height = 6)
      plot.new(); text(0.5, 0.5, "No multi-peak loci")
      grDevices::dev.off()
      return(invisible(FALSE))
    }
    meta_df <- dplyr::filter(meta_df, locus %in% locus_summary$locus)
  }

  locus_levels <- locus_summary$locus
  meta_df <- meta_df %>%
    dplyr::mutate(locus_rank = match(locus, locus_levels)) %>%
    dplyr::arrange(locus_rank, -max_val)

  dup_lab <- duplicated(meta_df$peak_label) |
             duplicated(meta_df$peak_label, fromLast = TRUE)
  meta_df$peak_label[dup_lab] <- paste0(meta_df$peak_label[dup_lab],
                                        " [", meta_df$locus[dup_lab], "]")

  peak_order <- meta_df$peak_label

  long_list <- lapply(seq_len(nrow(meta_df)), function(j) {
    i <- meta_df$idx[j]
    row_vals <- lineage_mat[i, ]
    lin_ord  <- order(row_vals, decreasing = TRUE)
    data.frame(
      peak_label = meta_df$peak_label[j],
      locus      = meta_df$locus[j],
      lineage    = colnames(lineage_mat)[lin_ord],
      value      = row_vals[lin_ord],
      stack_rank = seq_along(lin_ord),
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, long_list)
  long$peak_label <- factor(long$peak_label, levels = rev(peak_order))

  lin_totals <- tapply(long$value, long$lineage, sum, na.rm = TRUE)
  lin_levels <- names(sort(lin_totals, decreasing = TRUE))
  long$lineage <- factor(long$lineage, levels = lin_levels)
  long <- long %>% dplyr::arrange(peak_label, -value)

  if (!is.null(lineage_colours) && all(lin_levels %in% names(lineage_colours))) {
    fill_scale <- scale_fill_manual(values = lineage_colours)
  } else {
    fill_scale <- ggsci::scale_fill_npg()
  }

  n_per_locus <- meta_df %>%
    dplyr::group_by(locus_rank) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(locus_rank) %>%
    dplyr::pull(n)
  n_total <- sum(n_per_locus)
  if (length(n_per_locus) > 1) {
    rev_counts <- rev(n_per_locus)
    cum_bottom <- cumsum(rev_counts)
    sep_positions <- head(cum_bottom, -1) + 0.5
  } else {
    sep_positions <- numeric(0)
  }

  n_loci    <- nrow(locus_summary)
  n_multi   <- sum(locus_summary$n_peaks >= 2)
  subtitle_txt <- paste0(
    n_total, " peaks across ", n_loci, " loci (",
    n_multi, " multi-peak); locus separators are dotted lines"
  )

  p <- ggplot(long, aes(x = peak_label, y = value, fill = lineage,
                        group = factor(stack_rank))) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.85) +
    geom_hline(yintercept = cumvar_threshold,
               linetype = "dashed", color = "black", linewidth = 0.3) +
    coord_flip() +
    fill_scale +
    labs(
      title    = "Variant-peaks grouped by locus: lineage L2^2 decomposition",
      subtitle = subtitle_txt,
      x        = "Peak / CCV (grouped by locus)",
      y        = "Lineage L2^2 (sums to 1.0)",
      fill     = "Lineage"
    ) +
    theme_bw() +
    theme(
      axis.text.y  = element_text(size = 4),
      axis.text.x  = element_text(size = 6),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.title   = element_text(size = 10)
    )

  if (length(sep_positions)) {
    p <- p + geom_vline(xintercept = sep_positions,
                        linetype = "dotted", color = "grey40",
                        linewidth = 0.35)
  }

  pdf_height <- max(6, n_total * 0.14 + 3)
  grDevices::pdf(outfile, width = 12, height = min(pdf_height, 60))
  print(p)
  grDevices::dev.off()

  message("  Saved locus-grouped lineage stackbar (",
          n_total, " peaks, ", n_loci, " loci): ", outfile)
  invisible(peak_order)
}

.plotDominanceDodge <- function(genome_mat, gwas_mat,
                                 ratio_thresh = DEFAULT_RATIO_THRESH,
                                 title_label  = "",
                                 lineage_palette = NULL) {
  dom_genome <- addDominance(genome_mat, ratio_thresh = ratio_thresh)
  dom_gwas   <- addDominance(gwas_mat,   ratio_thresh = ratio_thresh)
  
  n_dom_genome <- sum(dom_genome$is_dominant)
  n_dom_gwas   <- sum(dom_gwas$is_dominant)
  
  df <- dplyr::bind_rows(
    dom_genome %>% dplyr::filter(is_dominant) %>% dplyr::count(maximum_ct) %>%
      dplyr::mutate(proportion = n / n_dom_genome, set = "Genome-wide"),
    dom_gwas   %>% dplyr::filter(is_dominant) %>% dplyr::count(maximum_ct) %>%
      dplyr::mutate(proportion = n / n_dom_gwas,   set = "GWAS")
  )
  
  df <- tidyr::complete(df, maximum_ct, set,
                        fill = list(n = 0L, proportion = 0))
  
  # One-sided proportion test: is GWAS proportion > genome-wide?
  all_cts <- unique(df$maximum_ct)
  prop_test_df <- lapply(all_cts, function(ct) {
    n_gwas_ct   <- df$n[df$maximum_ct == ct & df$set == "GWAS"]
    n_genome_ct <- df$n[df$maximum_ct == ct & df$set == "Genome-wide"]
    if (is.na(n_gwas_ct))   n_gwas_ct   <- 0L
    if (is.na(n_genome_ct)) n_genome_ct <- 0L
    pval <- tryCatch({
      prop.test(
        x = c(n_gwas_ct, n_genome_ct),
        n = c(n_dom_gwas, n_dom_genome),
        alternative = "greater"
      )$p.value
    }, error = function(e) NA_real_)
    data.frame(maximum_ct = ct, prop_p = pval, stringsAsFactors = FALSE)
  }) %>% do.call(rbind, .)
  
  df <- dplyr::left_join(df, prop_test_df, by = "maximum_ct")
  df$stars <- ifelse(df$set == "GWAS", .pvalStars(df$prop_p), "")
  
  df$set       <- factor(df$set, levels = c("Genome-wide", "GWAS"))
  df$bar_label <- ifelse(df$n > 0,
                         paste0(df$n, " (", round(df$proportion * 100, 1), "%)", df$stars), "")
  
  ct_order <- df %>% dplyr::filter(set == "GWAS") %>%
    dplyr::arrange(proportion) %>% dplyr::pull(maximum_ct)
  df$maximum_ct <- factor(df$maximum_ct, levels = ct_order)

  # Colour bars by cell type using lineage palette; pattern distinguishes sets
  if (!is.null(lineage_palette) && all(as.character(df$maximum_ct) %in% names(lineage_palette))) {
    fill_scale <- scale_fill_manual(values = lineage_palette)
    p <- ggplot(df, aes(x = maximum_ct, y = proportion, fill = maximum_ct, alpha = set)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = bar_label),
                position = position_dodge(width = 0.8), hjust = -0.05, size = 2.8) +
      coord_flip() +
      fill_scale +
      scale_alpha_manual(values = c("Genome-wide" = 0.4, "GWAS" = 1.0)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
      guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2)) +
      theme_bw()
  } else {
    p <- ggplot(df, aes(x = maximum_ct, y = proportion, fill = set)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = bar_label),
                position = position_dodge(width = 0.8), hjust = -0.05, size = 2.8) +
      coord_flip() +
      scale_fill_manual(values = c("Genome-wide" = "grey60", "GWAS" = "firebrick")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
      theme_bw()
  }
  
  p <- p + labs(
    title    = paste0("Dominant cell types (ratio > ", ratio_thresh, ") ", title_label),
    subtitle = paste0("Genome-wide: ", nrow(genome_mat), " peaks (", n_dom_genome,
                      " dominant);  GWAS: ", nrow(gwas_mat), " peaks (", n_dom_gwas,
                      " dominant)  [* p<0.05  ** p<0.01  *** p<0.001, one-sided proportion test]"),
    x    = "Cell type",
    y    = "Proportion within dominant peaks",
    fill = "Cell type",
    alpha = "Peak set"
  )

  p
}

.plotDominanceLineageDodge <- function(genome_mat, gwas_mat, cell_lineage,
                                         ratio_thresh = DEFAULT_RATIO_THRESH,
                                         title_label  = "",
                                         lineage_colours = NULL) {
  dom_genome <- addDominanceLineage(genome_mat, cell_lineage, ratio_thresh = ratio_thresh)
  dom_gwas <- addDominanceLineage(gwas_mat, cell_lineage, ratio_thresh = ratio_thresh)

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

  # One-sided proportion test: is GWAS proportion > genome-wide?
  all_lin <- unique(df$maximum_lineage)
  prop_test_df <- lapply(all_lin, function(lg) {
    n_gwas_lg   <- df$n[df$maximum_lineage == lg & df$set == "GWAS"]
    n_genome_lg <- df$n[df$maximum_lineage == lg & df$set == "Genome-wide"]
    if (is.na(n_gwas_lg))   n_gwas_lg   <- 0L
    if (is.na(n_genome_lg)) n_genome_lg <- 0L
    pval <- tryCatch({
      prop.test(
        x = c(n_gwas_lg, n_genome_lg),
        n = c(n_dom_gwas, n_dom_genome),
        alternative = "greater"
      )$p.value
    }, error = function(e) NA_real_)
    data.frame(maximum_lineage = lg, prop_p = pval, stringsAsFactors = FALSE)
  }) %>% do.call(rbind, .)

  df <- dplyr::left_join(df, prop_test_df, by = "maximum_lineage")
  df$stars <- ifelse(df$set == "GWAS", .pvalStars(df$prop_p), "")

  df$set <- factor(df$set, levels = c("Genome-wide", "GWAS"))
  df$bar_label <- ifelse(df$n > 0,
                         paste0(df$n, " (", round(df$proportion * 100, 1), "%)", df$stars), "")

  ct_order <- df %>% dplyr::filter(set == "GWAS") %>%
    dplyr::arrange(proportion) %>% dplyr::pull(maximum_lineage)
  df$maximum_lineage <- factor(df$maximum_lineage, levels = ct_order)

  # Colour bars by lineage using lineage colours; alpha distinguishes sets
  if (!is.null(lineage_colours) && all(as.character(df$maximum_lineage) %in% names(lineage_colours))) {
    fill_scale <- scale_fill_manual(values = lineage_colours,
                                     labels = function(x) sapply(x, function(lab) paste(strwrap(gsub(":", ", ", lab), width = 28), collapse = "\n"), USE.NAMES = FALSE))
    p <- ggplot(df, aes(x = maximum_lineage, y = proportion, fill = maximum_lineage, alpha = set)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = bar_label),
                position = position_dodge(width = 0.8), hjust = -0.05, size = 2.8) +
      coord_flip() +
      fill_scale +
      scale_alpha_manual(values = c("Genome-wide" = 0.4, "GWAS" = 1.0)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
      guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2)) +
      theme_bw()
  } else {
    p <- ggplot(df, aes(x = maximum_lineage, y = proportion, fill = set)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = bar_label),
                position = position_dodge(width = 0.8), hjust = -0.05, size = 2.8) +
      coord_flip() +
      scale_fill_manual(values = c("Genome-wide" = "grey60", "GWAS" = "firebrick")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
      theme_bw()
  }

  # Wrap long lineage labels (colon-separated cell types) to prevent overflow
  wrap_lineage <- function(x, width = 28) {
    sapply(x, function(lab) paste(strwrap(gsub(":", ", ", lab), width = width), collapse = "\n"),
           USE.NAMES = FALSE)
  }

  p <- p + labs(
    title = paste0("Dominant lineages (ratio > ", ratio_thresh, ") ", title_label),
    subtitle = paste0("Genome-wide: ", nrow(genome_mat), " peaks (", n_dom_genome,
                      " dominant);  GWAS: ", nrow(gwas_mat), " peaks (", n_dom_gwas,
                      " dominant)  [* p<0.05  ** p<0.01  *** p<0.001, one-sided proportion test]"),
    x = "Lineage",
    y = "Proportion within dominant peaks",
    fill = "Lineage",
    alpha = "Peak set"
  ) +
    scale_x_discrete(labels = wrap_lineage)
  p
}

plotSpecificityBar <- function(specificity_summary_l2, se_name, trait_name, plots_dir) {
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

# ---------------------------------------------------------------
# Plot 1: Density of max L2 specificity per GWAS peak
# Shows the distribution of maximum_score (max L2 across cell types)
# for GWAS peaks vs genome-wide peaks, with threshold lines.
# ---------------------------------------------------------------
plotSpecificityDensity <- function(mat_l2, gwas_mat, se_name, trait_name, plots_dir) {
  if (!nrow(gwas_mat)) return(invisible(NULL))
  n_cells <- ncol(gwas_mat)
  high_thresh <- 1 / sqrt(2)
  mid_thresh  <- 1 / sqrt(n_cells)

  gwas_max  <- apply(gwas_mat, 1, max)
  genome_max <- apply(mat_l2, 1, max)
  df <- rbind(
    data.frame(max_l2 = genome_max, set = "Genome-wide"),
    data.frame(max_l2 = gwas_max,   set = "GWAS")
  )

  p <- ggplot(df, aes(x = max_l2, fill = set, colour = set)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    geom_vline(xintercept = high_thresh, linetype = "dashed", colour = "red") +
    geom_vline(xintercept = mid_thresh,  linetype = "dotted", colour = "orange") +
    annotate("text", x = high_thresh, y = Inf, label = "high", vjust = 2, hjust = -0.1, colour = "red", size = 3) +
    annotate("text", x = mid_thresh,  y = Inf, label = "mid",  vjust = 2, hjust = -0.1, colour = "orange", size = 3) +
    scale_fill_manual(values = c("Genome-wide" = "grey70", "GWAS" = "#2166AC")) +
    scale_colour_manual(values = c("Genome-wide" = "grey40", "GWAS" = "#2166AC")) +
    labs(
      title = paste0(se_name, " x ", trait_name, " — Max L2 specificity per peak"),
      subtitle = paste0("high > 1/sqrt(2) = ", round(high_thresh, 3),
                        ";  mid > 1/sqrt(", n_cells, ") = ", round(mid_thresh, 3)),
      x = "Max L2 specificity (across cell types)",
      y = "Density",
      fill = "Peak set", colour = "Peak set"
    ) +
    theme_bw()

  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_l2_specificity_density.pdf")),
         p, width = 7, height = 5)
  invisible(p)
}

# ---------------------------------------------------------------
# Plot 2: Weighted specificity — max_L2 x weight (mean raw accessibility)
# Scatter of max L2 vs weight for GWAS peaks. Peaks that are both
# specific AND accessible are top-right and most biologically relevant.
# Weight-scaled specificity = max_L2 * weight.
# ---------------------------------------------------------------
plotSpecificityWeighted <- function(specificity_summary_l2, se_name, trait_name, plots_dir,
                                      lineage_palette = NULL) {
  if (!nrow(specificity_summary_l2) || all(is.na(specificity_summary_l2$weight))) return(invisible(NULL))
  n_cells <- max(specificity_summary_l2$n_high + specificity_summary_l2$n_mid + specificity_summary_l2$n_low, na.rm = TRUE)
  high_thresh <- 1 / sqrt(2)

  df <- specificity_summary_l2 %>%
    dplyr::mutate(
      specificity_tier = dplyr::case_when(
        is_dominant ~ "Dominant (ratio >= 2)",
        maximum_score > high_thresh ~ "Highly specific",
        TRUE ~ "Other"
      ),
      specificity_tier = factor(specificity_tier,
                                levels = c("Dominant (ratio >= 2)", "Highly specific", "Other")),
      weighted_specificity = maximum_score * weight
    )

  tier_colours <- c(
    "Dominant (ratio >= 2)" = "#B2182B",
    "Highly specific" = "#2166AC",
    "Other" = "grey60"
  )

  p <- ggplot(df, aes(x = weight, y = maximum_score)) +
    geom_point(aes(colour = specificity_tier), alpha = 0.6, size = 2) +
    geom_hline(yintercept = high_thresh, linetype = "dashed", colour = "red", alpha = 0.5) +
    scale_colour_manual(values = tier_colours) +
    labs(
      title = paste0(se_name, " x ", trait_name, " — Specificity vs Accessibility"),
      subtitle = "Top-right = specific AND accessible peaks (most biologically relevant)",
      x = "Mean raw accessibility (weight = rowMean of raw counts)",
      y = "Max L2 specificity",
      colour = "Specificity tier"
    ) +
    theme_bw()

  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_l2_specificity_vs_weight.pdf")),
         p, width = 8, height = 6)
  invisible(p)
}

# ---------------------------------------------------------------
# Plot: Weight (mean raw accessibility) vs distance to nearest TSS
# Colour encodes TSS bin (core_prom, prox_prom, near_reg, distal, long_range).
# ---------------------------------------------------------------
plotWeightVsTss <- function(specificity_summary_l2, se_name, trait_name, plots_dir) {
  if (!nrow(specificity_summary_l2) ||
      all(is.na(specificity_summary_l2$dist_to_tss)) ||
      all(is.na(specificity_summary_l2$weight))) return(invisible(NULL))

  bin_levels <- c("core_prom", "prox_prom", "near_reg", "distal", "long_range")

  df <- specificity_summary_l2 %>%
    dplyr::filter(!is.na(dist_to_tss), !is.na(weight)) %>%
    dplyr::mutate(
      tss_bin = factor(tss_bin, levels = bin_levels),
      log_dist = log10(dist_to_tss + 1)
    )

  p <- ggplot(df, aes(x = log_dist, y = weight, colour = tss_bin)) +
    geom_point(alpha = 0.6, size = 2.2) +
    ggsci::scale_colour_igv() +
    geom_vline(xintercept = log10(c(1000, 5000, 50000, 200000)),
               linetype = "dotted", colour = "grey50", linewidth = 0.3) +
    labs(
      title = paste0(se_name, " x ", trait_name, " \u2014 Weight vs TSS distance"),
      subtitle = "Mean raw accessibility by distance to nearest TSS",
      x = expression(log[10](distance~to~nearest~TSS + 1)),
      y = "Mean raw accessibility (weight)",
      colour = "TSS bin"
    ) +
    theme_bw() +
    theme(legend.position = "right")

  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_weight_vs_tss.pdf")),
         p, width = 8, height = 6)
  invisible(p)
}

# ---------------------------------------------------------------
# Plot: Lineage-level inter-peak correlation heatmap
# Computes peak-by-lineage matrix (L2^2 summed per lineage group)
# then shows the inter-lineage Pearson correlation.
# ---------------------------------------------------------------
plotLineageCorHeatmap <- function(mat_l2, cell_lineage, label, plots_dir,
                                     lineage_colours = NULL) {
  lineage_mat <- getLineageMatrix(mat_l2, cell_lineage)
  if (ncol(lineage_mat) < 2) return(invisible(NULL))
  cor_lin <- cor(lineage_mat, use = "pairwise.complete.obs")

  # Fixed [-1, 1] symmetric scale anchored at 0 (white) so this heatmap
  # is directly comparable in colour with the cell-type L2 correlation
  # heatmap (plotAssayCorHeatmaps, raw_l2 branch). Actual observed
  # range is reported in the column title for context.
  rng <- range(cor_lin, na.rm = TRUE)
  col_fun <- circlize::colorRamp2(c(-1, 0, 1),
                                  c("#2166AC", "white", "#B2182B"))

  # Row/column colour annotation from lineage palette
  ha <- NULL
  if (!is.null(lineage_colours)) {
    lg_names <- colnames(lineage_mat)
    lg_cols <- lineage_colours[lg_names]
    lg_cols <- lg_cols[!is.na(lg_cols)]
    if (length(lg_cols) == length(lg_names)) {
      ha <- ComplexHeatmap::HeatmapAnnotation(
        Lineage = lg_names,
        col = list(Lineage = setNames(lg_cols, lg_names)),
        show_legend = FALSE
      )
    }
  }

  out <- file.path(plots_dir, label)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  fname <- file.path(out, paste0(label, "_lineage_cor_heatmap.pdf"))
  grDevices::pdf(fname, width = 8, height = 7)
  ht <- ComplexHeatmap::Heatmap(
    cor_lin,
    name = "Pearson r",
    col = col_fun,
    column_title = paste0(label, " — Lineage L2^2 correlation",
                          "  [range: ", round(rng[1], 2), " to ", round(rng[2], 2), "]"),
    top_annotation = ha,
    show_row_dend = TRUE, show_column_dend = TRUE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 7)
  )
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()
  invisible(fname)
}

# ---------------------------------------------------------------
# Plot: Per-lineage L2^2 density distributions
# Shows the distribution of L2^2 values for each lineage group,
# coloured by lineage master hues.
# ---------------------------------------------------------------
plotLineageDensities <- function(mat_l2, cell_lineage, label, plots_dir,
                                    lineage_colours = NULL) {
  lineage_mat <- getLineageMatrix(mat_l2, cell_lineage)
  if (ncol(lineage_mat) < 1) return(invisible(NULL))

  long <- as.data.frame(lineage_mat) %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "lineage", values_to = "L2_sq")

  if (!is.null(lineage_colours) &&
      all(unique(long$lineage) %in% names(lineage_colours))) {
    colour_scale <- scale_colour_manual(values = lineage_colours)
  } else {
    colour_scale <- ggsci::scale_colour_npg()
  }

  p <- ggplot(long, aes(x = L2_sq, colour = lineage)) +
    geom_density(linewidth = 0.7) +
    colour_scale +
    labs(
      title = paste0(label, " — Lineage L2^2 density distributions"),
      subtitle = "L2^2 summed within lineage groups (sums to 1 per peak)",
      x = "Lineage L2^2",
      y = "Density",
      colour = "Lineage"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  out <- file.path(plots_dir, label)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  fname <- file.path(out, paste0(label, "_lineage_densities.pdf"))
  ggsave(fname, p, width = 8, height = 5)
  invisible(fname)
}

plotDominance <- function(mat_l2, gwas_mat, se_name, trait_name, plots_dir, lineage_palette = NULL) {
  p_dom_l2 <- .plotDominanceDodge(mat_l2, gwas_mat, title_label = "L2mat", lineage_palette = lineage_palette)
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_dominance_l2.pdf")), p_dom_l2, width = 8, height = 6)
  invisible(p_dom_l2)
}

plotDominanceLineage <- function(mat_l2, gwas_mat, cell_lineage, se_name, trait_name, plots_dir,
                                    lineage_colours = NULL) {
  p_dom_lineage <- .plotDominanceLineageDodge(mat_l2, gwas_mat, cell_lineage = cell_lineage,
                                                title_label = "lineage-summed L2^2",
                                                lineage_colours = lineage_colours)
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_dominance_lineage_l2.pdf")), p_dom_lineage, width = 8, height = 6)
  invisible(p_dom_lineage)
}

plotCoherence <- function(locus_results, se_name, trait_name, plots_dir) {
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

  # --- Plot 2: Dominance fraction volcano (lineage-level) ---
  # Shows loci by the fraction of peaks that are individually lineage-dominant
  # vs the permutation p-value. Based on the best-performing metric from marker-gene
  # validation (dominant_score_l_sqrt).
  if (all(c("z_score_dom_lin", "p_value_dom_lin", "dominant_frac_lineage", "dominant_category") %in% names(plot_df))) {
    dom_df <- plot_df %>%
      dplyr::filter(!is.na(z_score_dom_lin), !is.na(p_value_dom_lin), p_value_dom_lin > 0) %>%
      dplyr::mutate(
        neg_log10_p_dom = -log10(p_value_dom_lin),
        dominant_category = factor(dominant_category,
          levels = c("Highly_dominant_specific", "Dominant_moderate", "Intermediate", "Non_dominant", "Single_peak"))
      ) %>%
      dplyr::filter(is.finite(neg_log10_p_dom))

    if (nrow(dom_df)) {
      y_cap_dom <- max(as.numeric(stats::quantile(dom_df$neg_log10_p_dom, 0.99, na.rm = TRUE)),
                       -log10(0.05) + 0.2)
      p_dom_volcano <- ggplot(dom_df, aes(x = dominant_frac_lineage, y = neg_log10_p_dom)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.55, color = "gray35") +
        geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.5, color = "gray45") +
        geom_point(aes(size = n_peaks, color = dominant_category, shape = dominant_lineage), alpha = 0.88) +
        scale_color_manual(
          values = c(
            Highly_dominant_specific = "#C43C39",
            Dominant_moderate        = "#E7965A",
            Intermediate             = "#5A5A5A",
            Non_dominant             = "#2F78B7",
            Single_peak              = "#bdbdbd"
          ), drop = FALSE) +
        scale_shape_discrete(name = "Dominant lineage") +
        scale_size_continuous(range = c(2.3, 8.2), breaks = c(2, 5, 10, 20, 30)) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                           expand = expansion(mult = c(0.02, 0.02))) +
        coord_cartesian(ylim = c(0, y_cap_dom)) +
        theme_minimal(base_size = 12) +
        theme(
          legend.position  = "right",
          plot.title       = element_text(face = "bold", size = 13),
          plot.subtitle    = element_text(color = "gray30", size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
          panel.grid.major.y = element_line(color = "gray88", linewidth = 0.3),
          panel.border     = element_rect(fill = NA, color = "gray40", linewidth = 0.45)
        ) +
        labs(
          x        = "Dominant peak fraction (lineage-level)",
          y        = expression(-log[10](p-value)),
          title    = paste0(se_name, " x ", trait_name, " - Locus dominance (lineage-level)"),
          subtitle = paste0("Fraction of peaks with lineage ratio > ", DEFAULT_RATIO_THRESH,
                            "; n = ", nrow(dom_df), " loci (y-axis capped at 99th percentile)"),
          color    = "Dominant category",
          shape    = "Dominant lineage",
          size     = "N peaks"
        )
      ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_locus_dominance_lineage.pdf")),
             p_dom_volcano, width = 8, height = 6)
    }
  }

  # --- Plot 3: Spearman concordance Z vs dominance Z (method comparison) ---
  # Loci that score high on both metrics are the most robustly specific.
  # Discordant loci (high on one, low on other) reveal method-specific signal.
  if (all(c("z_score_dom_lin", "z_score") %in% names(plot_df))) {
    cmp_df <- plot_df %>%
      dplyr::filter(!is.na(z_score_dom_lin), !is.na(z_score)) %>%
      dplyr::mutate(dominant_lineage = as.character(dominant_lineage))

    if (nrow(cmp_df)) {
      p_method_cmp <- ggplot(cmp_df, aes(x = z_score, y = z_score_dom_lin)) +
        geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray60") +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray60") +
        geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.4, color = "gray40") +
        geom_point(aes(size = n_peaks, colour = dominant_lineage), alpha = 0.80) +
        scale_size_continuous(range = c(2.0, 7.5), breaks = c(2, 5, 10, 20, 30)) +
        theme_minimal(base_size = 12) +
        theme(
          legend.position  = "right",
          plot.title       = element_text(face = "bold", size = 13),
          plot.subtitle    = element_text(color = "gray30", size = 10),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(fill = NA, color = "gray40", linewidth = 0.45)
        ) +
        labs(
          x        = "Spearman concordance Z-score",
          y        = "Dominance fraction Z-score (lineage-level)",
          title    = paste0(se_name, " x ", trait_name, " - Concordance vs Dominance"),
          subtitle = "Dotted diagonal = perfect agreement between methods; each point is a locus",
          colour   = "Dominant lineage",
          size     = "N peaks"
        )
      ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_locus_concordance_vs_dominance.pdf")),
             p_method_cmp, width = 8, height = 6)
    }
  }

  invisible(p_locus_cat)
}


# ============================================================
# plotSignalComparison
# Per-signal comparison of all GWAS peaks vs dominant-filtered
# peaks.  Four panels:
#   1. Per-signal peak counts (all / cell-dominant / lineage-dominant)
#   2. Mean weight_max per signal coloured by lineage-dominant fraction
#   3. Concordance scatter: all peaks vs lineage-dominant peaks
#   4. Z-score scatter: all peaks vs lineage-dominant peaks
# ============================================================
plotSignalComparison <- function(signal_summary,
                                   locus_results_dom_lineage = NULL,
                                   se_name, trait_name, plots_dir) {
  if (!nrow(signal_summary)) return(invisible(NULL))

  ss <- signal_summary %>%
    dplyr::filter(!is.na(locus), nzchar(locus))
  if (!nrow(ss)) return(invisible(NULL))

  # ---- Panel 1: Peak counts per signal ----
  count_long <- ss %>%
    dplyr::select(locus, n_peaks, n_dominant_cell, n_dominant_lineage) %>%
    tidyr::pivot_longer(-locus, names_to = "metric", values_to = "n") %>%
    dplyr::mutate(
      metric = factor(metric,
        levels = c("n_peaks", "n_dominant_cell", "n_dominant_lineage"),
        labels = c("All peaks", "Cell-dominant", "Lineage-dominant"))
    )
  locus_order <- ss %>% dplyr::arrange(-n_peaks) %>% dplyr::pull(locus)
  count_long$locus <- factor(count_long$locus, levels = rev(locus_order))

  p1 <- ggplot(count_long, aes(x = locus, y = n, fill = metric)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c(
      "All peaks"        = "#4393C3",
      "Cell-dominant"    = "#D6604D",
      "Lineage-dominant" = "#F4A582")) +
    theme_bw(base_size = 9) +
    labs(
      title   = paste0(se_name, " \u00d7 ", trait_name, " \u2014 Peak counts per signal"),
      x = "Signal / locus", y = "Number of peaks", fill = NULL
    ) +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = max(4, 8 - nrow(ss) %/% 10)))

  # ---- Panel 2: Mean weight_max per signal ----
  p2 <- NULL
  if ("mean_weight_max" %in% names(ss) && !all(is.na(ss$mean_weight_max))) {
    fd_col <- if ("frac_dominant_lineage" %in% names(ss)) "frac_dominant_lineage" else NULL
    ss2 <- ss %>%
      dplyr::filter(!is.na(mean_weight_max)) %>%
      dplyr::arrange(-mean_weight_max) %>%
      dplyr::mutate(locus = factor(locus, levels = rev(locus)))

    fill_aes <- if (!is.null(fd_col)) aes(fill = frac_dominant_lineage) else aes(fill = n_peaks)
    p2 <- ggplot(ss2, aes(x = locus, y = mean_weight_max)) +
      geom_col(fill_aes, width = 0.75) +
      coord_flip() +
      scale_fill_gradient(
        low = "grey82", high = "#D6604D",
        limits = if (!is.null(fd_col)) c(0, 1) else NULL,
        name = if (!is.null(fd_col)) "Frac.\ndom.\n(lineage)" else "N peaks"
      ) +
      theme_bw(base_size = 9) +
      labs(
        title = paste0(se_name, " \u00d7 ", trait_name, " \u2014 Mean max accessibility per signal"),
        x = "Signal / locus", y = "Mean weight_max (max raw accessibility)"
      ) +
      theme(axis.text.y = element_text(size = max(4, 8 - nrow(ss2) %/% 10)))
  }

  # ---- Panels 3 & 4: Concordance comparison (requires dominant lineage results) ----
  p3 <- NULL; p4 <- NULL
  has_concordance <- "concordance" %in% names(ss) && !all(is.na(ss$concordance))
  has_dom         <- !is.null(locus_results_dom_lineage) && nrow(locus_results_dom_lineage) > 0

  if (has_concordance && has_dom) {
    comp_df <- ss %>%
      dplyr::select(locus, concordance_all = concordance, z_all = z_score, n_peaks_all = n_peaks) %>%
      dplyr::inner_join(
        locus_results_dom_lineage %>%
          dplyr::select(locus, concordance_dom = concordance, z_dom = z_score, n_peaks_dom = n_peaks),
        by = "locus"
      ) %>%
      dplyr::filter(!is.na(concordance_all), !is.na(concordance_dom))

    if (nrow(comp_df) >= 2L) {
      lim_conc <- range(c(comp_df$concordance_all, comp_df$concordance_dom), na.rm = TRUE)
      p3 <- ggplot(comp_df, aes(x = concordance_all, y = concordance_dom)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
        geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
        geom_point(aes(size = n_peaks_all, colour = n_peaks_dom), alpha = 0.82) +
        scale_size_continuous(range = c(2, 8), name = "N peaks\n(all)") +
        scale_colour_gradient(low = "grey78", high = "#D6604D", name = "N peaks\n(lin.-dom.)") +
        xlim(lim_conc) + ylim(lim_conc) +
        theme_bw(base_size = 10) +
        labs(
          title    = paste0(se_name, " \u00d7 ", trait_name),
          subtitle = "Concordance: all peaks vs lineage-dominant peaks",
          x = "Spearman concordance (all peaks)",
          y = "Spearman concordance (lineage-dominant peaks)"
        )

      lim_z <- range(c(comp_df$z_all, comp_df$z_dom), na.rm = TRUE)
      p4 <- ggplot(comp_df, aes(x = z_all, y = z_dom)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
        geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
        geom_point(aes(size = n_peaks_dom, colour = n_peaks_all), alpha = 0.82) +
        scale_size_continuous(range = c(2, 8), name = "N peaks\n(lin.-dom.)") +
        scale_colour_gradient(low = "grey78", high = "#4393C3", name = "N peaks\n(all)") +
        xlim(lim_z) + ylim(lim_z) +
        theme_bw(base_size = 10) +
        labs(
          subtitle = "Z-score: all peaks vs lineage-dominant peaks",
          x = "Concordance Z-score (all peaks)",
          y = "Concordance Z-score (lineage-dominant peaks)"
        )
    }
  }

  # ---- Assemble and save ----
  # Page 1: per-signal bar plots (p1, p2) — tall, scales with number of signals.
  # Page 2: concordance scatter plots (p3, p4) — compact, independent of signal
  # count (scatters plot one point per locus, not one row per signal).
  bar_plots     <- Filter(Negate(is.null), list(p1, p2))
  scatter_plots <- Filter(Negate(is.null), list(p3, p4))
  if (length(bar_plots) + length(scatter_plots) == 0L) return(invisible(NULL))

  fname  <- file.path(plots_dir, paste0(se_name, "_", trait_name, "_signal_comparison.pdf"))
  ph_bar <- min(max(5, nrow(ss) * 0.25 + 3), 60)  # bar-plot height scales with signal count
  ph_sc  <- 6                                      # fixed compact height for scatter page

  # Use a single device, one page per plot group. PDF page size is set by the
  # pdf() call, so we size it to the taller page and pad the scatter page so
  # the scatters stay small in absolute terms.
  page_h <- max(ph_bar, ph_sc)
  grDevices::pdf(fname, width = 12, height = page_h, onefile = TRUE)

  # Page 1: bar plots
  if (length(bar_plots) > 0L) {
    bar_grob <- if (length(bar_plots) >= 2L) {
      gridExtra::arrangeGrob(grobs = bar_plots, ncol = 2L)
    } else {
      gridExtra::arrangeGrob(grobs = bar_plots, ncol = 1L)
    }
    # Pad so bar plots always fill their natural height even if page is taller
    pad_bar <- max(0, (page_h - ph_bar) / 2)
    blank   <- grid::rectGrob(gp = grid::gpar(col = NA, fill = NA))
    page1   <- gridExtra::arrangeGrob(
      blank, bar_grob, blank,
      ncol = 1L,
      heights = grid::unit(c(pad_bar, ph_bar, pad_bar), "in")
    )
    grid::grid.newpage()
    grid::grid.draw(page1)
  }

  # Page 2: scatter plots (fixed compact size, padded vertically)
  if (length(scatter_plots) > 0L) {
    scatter_grob <- if (length(scatter_plots) >= 2L) {
      gridExtra::arrangeGrob(grobs = scatter_plots, ncol = 2L)
    } else {
      gridExtra::arrangeGrob(grobs = scatter_plots, ncol = 1L)
    }
    pad_sc <- max(0, (page_h - ph_sc) / 2)
    blank  <- grid::rectGrob(gp = grid::gpar(col = NA, fill = NA))
    page2  <- gridExtra::arrangeGrob(
      blank, scatter_grob, blank,
      ncol = 1L,
      heights = grid::unit(c(pad_sc, ph_sc, pad_sc), "in")
    )
    grid::grid.newpage()
    grid::grid.draw(page2)
  }

  grDevices::dev.off()

  message("  Saved signal comparison plot: ", fname)
  invisible(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}
