# ============================================================
# plotting.R
# All visualisation functions for CLEAR pipeline
# ============================================================

pval_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p)    ~ "",
    p < 0.001   ~ "***",
    p < 0.01    ~ "**",
    p < 0.05    ~ "*",
    TRUE        ~ ""
  )
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

dominance_dodge_plot <- function(genome_mat, gwas_mat,
                                 ratio_thresh = DEFAULT_RATIO_THRESH,
                                 title_label  = "",
                                 lineage_palette = NULL) {
  dom_genome <- compute_dominance(genome_mat, ratio_thresh = ratio_thresh)
  dom_gwas   <- compute_dominance(gwas_mat,   ratio_thresh = ratio_thresh)
  
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
  df$stars <- ifelse(df$set == "GWAS", pval_to_stars(df$prop_p), "")
  
  df$set       <- factor(df$set, levels = c("Genome-wide", "GWAS"))
  df$bar_label <- ifelse(df$n > 0,
                         paste0(df$n, " (", round(df$proportion * 100, 1), "%)", df$stars), "")
  
  ct_order <- df %>% dplyr::filter(set == "GWAS") %>%
    dplyr::arrange(proportion) %>% dplyr::pull(maximum_ct)
  df$maximum_ct <- factor(df$maximum_ct, levels = ct_order)
  
  p <- ggplot(df, aes(x = maximum_ct, y = proportion, fill = set)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = bar_label),
              position = position_dodge(width = 0.8), hjust = -0.05, size = 2.8) +
    coord_flip() +
    scale_fill_manual(values = c("Genome-wide" = "grey60", "GWAS" = "firebrick")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
    theme_bw() +
    labs(
      title    = paste0("Dominant cell types (ratio > ", ratio_thresh, ") ", title_label),
      subtitle = paste0("Genome-wide: ", nrow(genome_mat), " peaks (", n_dom_genome,
                        " dominant);  GWAS: ", nrow(gwas_mat), " peaks (", n_dom_gwas,
                        " dominant)  [* p<0.05  ** p<0.01  *** p<0.001, one-sided proportion test]"),
      x    = "Cell type",
      y    = "Proportion within dominant peaks",
      fill = "Peak set"
    )

  if (!is.null(lineage_palette)) {
    ct_cols <- lineage_palette[levels(df$maximum_ct)]
    ct_cols <- ct_cols[!is.na(ct_cols)]
    if (length(ct_cols)) {
      p <- p + theme(axis.text.y = element_text(colour = ct_cols))
    }
  }
  p
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
  df$stars <- ifelse(df$set == "GWAS", pval_to_stars(df$prop_p), "")

  df$set <- factor(df$set, levels = c("Genome-wide", "GWAS"))
  df$bar_label <- ifelse(df$n > 0,
                         paste0(df$n, " (", round(df$proportion * 100, 1), "%)", df$stars), "")

  ct_order <- df %>% dplyr::filter(set == "GWAS") %>%
    dplyr::arrange(proportion) %>% dplyr::pull(maximum_lineage)
  df$maximum_lineage <- factor(df$maximum_lineage, levels = ct_order)

  ggplot(df, aes(x = maximum_lineage, y = proportion, fill = set)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = bar_label),
              position = position_dodge(width = 0.8), hjust = -0.05, size = 2.8) +
    coord_flip() +
    scale_fill_manual(values = c("Genome-wide" = "grey60", "GWAS" = "firebrick")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
    theme_bw() +
    labs(
      title = paste0("Dominant lineages (ratio > ", ratio_thresh, ") ", title_label),
      subtitle = paste0("Genome-wide: ", nrow(genome_mat), " peaks (", n_dom_genome,
                        " dominant);  GWAS: ", nrow(gwas_mat), " peaks (", n_dom_gwas,
                        " dominant)  [* p<0.05  ** p<0.01  *** p<0.001, one-sided proportion test]"),
      x = "Lineage",
      y = "Proportion within dominant peaks",
      fill = "Peak set"
    )
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
