# ============================================================
# global_enrichment.R
# Trait-level rank enrichment test (CHEERS-style) at cell-type
# or lineage resolution. Tests whether GWAS-overlapping peaks
# have higher specificity ranks than expected by chance.
#
# Method:
#   1. Rank all peaks genome-wide by specificity for each column
#      (cell type or lineage), ascending: rank N = most specific
#   2. Take the mean rank of GWAS-overlapping peaks per column
#   3. Compare to the uniform null expectation (N+1)/2
#   4. Z = (observed - expected) / SE, SE = sqrt((N^2-1)/(12*n))
#   5. One-sided p-value (upper tail), BH FDR across columns
#
# TIES. CLEAR/R/global_enrichment.R ranks with ties.method = "min".
# That breaks the null used on line 3: with a tied block of size m
# starting at position r, all m peaks get rank r, so the mean rank
# over ALL peaks falls below (N+1)/2 while the null stays at
# (N+1)/2. A specificity column has a large tie block at zero (every
# peak with no signal in that cell type), so the bias is real and,
# worse, DIFFERENTIAL: a cell type with many zero peaks is pushed
# towards negative z even for a random peak set, and one with few is
# favoured. Since comparing cell types is the entire point, that is
# not a constant offset that cancels.
#
# "average" is used here instead. Average ranks always sum to
# N(N+1)/2, so the population mean is exactly (N+1)/2 and the null
# on line 3 is correct whatever the tie structure. The variance is
# slightly smaller than (N^2-1)/12 when ties exist, which makes the
# test mildly conservative -- the safe direction. Set
# ties = "min" to reproduce CLEAR's numbers.
#
# The SE on line 4 treats the n GWAS peaks as an iid sample and
# ignores the finite-population correction (N-n)/(N-1); with n in
# the hundreds and N in the hundreds of thousands that is a
# sub-percent effect. What it does NOT handle is that credible-set
# peaks are correlated through LD, which no variance formula here
# can fix.
# ============================================================

addGlobalEnrichment <- function(se, snps, mat_l2,
                                results_dir, plots_dir,
                                se_name, trait_name,
                                level = c("cell", "lineage"),
                                cell_lineage   = NULL,
                                lineage_palette = NULL,
                                lineage_colours = NULL,
                                ties = c("average", "min"),
                                write_z_plot = FALSE) {
  level <- match.arg(level)
  ties  <- match.arg(ties)
  mat_l2 <- .specMatrix(mat_l2)
  if (nrow(mat_l2) != nrow(se) || !identical(rownames(mat_l2), rownames(se))) {
    stop("mat_l2 rows must match se in name and order.")
  }
  if (nrow(mat_l2) < 2L) stop("Global enrichment requires at least two background peaks.")

  if (level == "lineage") {
    if (is.null(cell_lineage) || length(cell_lineage) == 0) {
      warning("cell_lineage missing/empty; skipping lineage global enrichment.")
      return(invisible(data.frame()))
    }
    score_mat <- getLineageMatrix(mat_l2, cell_lineage)
    if (ncol(score_mat) < 2) {
      warning("Fewer than 2 lineages; skipping lineage global enrichment.")
      return(invisible(data.frame()))
    }
    x_label      <- "Lineage"
    id_col       <- "lineage"
    file_suffix  <- "_lineage"
    title_suffix <- " (Lineage)"
    fill_palette <- lineage_colours
    output_file  <- "global_enrichment_lineage.txt"
  } else {
    score_mat    <- mat_l2
    x_label      <- "Cell type"
    id_col       <- "celltype"
    file_suffix  <- ""
    title_suffix <- ""
    fill_palette <- lineage_palette
    output_file  <- "global_enrichment.txt"
  }

  # Column-wise rank across all peaks. See the note on ties above:
  # "average" keeps the population mean rank at (N+1)/2, which is the
  # null this test compares against.
  nmat_r <- vapply(seq_len(ncol(score_mat)),
                   function(j) rank(score_mat[, j], ties.method = ties),
                   numeric(nrow(score_mat)))
  nmat_r <- matrix(nmat_r, nrow(score_mat), ncol(score_mat),
                   dimnames = dimnames(score_mat))
  nse_rank <- SummarizedExperiment(nmat_r, rowRanges = rowRanges(se))

  gwas_rank_idx <- unique(subjectHits(findOverlaps(snps, rowRanges(nse_rank))))
  gwas_rank_mat <- if (length(gwas_rank_idx))
    assay(nse_rank[gwas_rank_idx, , drop = FALSE], 1) else
    nmat_r[0, , drop = FALSE]
  if (!nrow(gwas_rank_mat)) {
    message("No GWAS-overlapping peaks; global enrichment is not evaluated.")
    ge <- data.frame(observed_mean_rank = numeric(), n_peaks = integer(),
                      exp_rank = numeric(), obs.exp = numeric(), z = numeric(),
                      p = numeric(), fdr = numeric(), positive = logical())
    ge[[id_col]] <- character()
    ge$fdr_stars <- character()
    .specWrite(ge, file.path(results_dir, output_file))
    return(invisible(ge))
  }

  N <- nrow(nse_rank)
  n <- nrow(gwas_rank_mat)
  exp_mean <- (N + 1) / 2
  se_z <- sqrt((N^2 - 1) / (12 * n))

  ge <- data.frame(observed_mean_rank = colMeans(gwas_rank_mat),
                   stringsAsFactors = FALSE)
  ge$n_peaks   <- n
  ge$exp_rank  <- exp_mean
  ge$obs.exp   <- ge$observed_mean_rank / exp_mean
  ge$z         <- (ge$observed_mean_rank - exp_mean) / se_z
  ge$p         <- pnorm(ge$z, lower.tail = FALSE)
  ge$fdr       <- p.adjust(ge$p, method = "BH")
  ge$positive  <- ge$z > 0 & ge$fdr < 0.05
  ge[[id_col]] <- rownames(ge)

  ge$fdr_stars <- dplyr::case_when(
    is.na(ge$fdr)  ~ "",
    ge$fdr < 0.001 ~ "***",
    ge$fdr < 0.01  ~ "**",
    ge$fdr < 0.05  ~ "*",
    TRUE           ~ ""
  )

  .specWrite(ge, file.path(results_dir, output_file))
  if (is.null(plots_dir)) return(ge)

  ids <- ge[[id_col]]
  use_palette <- !is.null(fill_palette) && all(ids %in% names(fill_palette))

  title_z   <- paste0(se_name, " x ", trait_name, " - Global GWAS Enrichment", title_suffix)
  title_fdr <- paste0(se_name, " x ", trait_name, " - Global Enrichment", title_suffix, " (FDR)")
  subtitle  <- "* FDR<0.05  ** FDR<0.01  *** FDR<0.001"

  # --- Z-score bar plot ---
  if (use_palette) {
    ge_ord <- ge[order(ge$z), ]
    ge_ord[[id_col]] <- factor(ge_ord[[id_col]], levels = ge_ord[[id_col]])
    p_enrich <- ggplot(ge_ord, aes(x = .data[[id_col]], y = z, fill = .data[[id_col]])) +
      geom_bar(stat = "identity") + coord_flip() + theme_bw() +
      geom_text(aes(label = fdr_stars), hjust = ifelse(ge_ord$z >= 0, -0.3, 1.3), size = 4) +
      scale_fill_manual(values = fill_palette, labels = .wrapCellList) +
      scale_x_discrete(labels = .wrapCellList) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
      labs(x = x_label, y = "Z-score (rank enrichment)",
           title = title_z, subtitle = subtitle, fill = x_label)
  } else {
    p_enrich <- ggplot(ge, aes(x = reorder(.data[[id_col]], z), y = z, fill = positive)) +
      geom_bar(stat = "identity") + coord_flip() + theme_bw() +
      geom_text(aes(label = fdr_stars), hjust = ifelse(ge$z >= 0, -0.3, 1.3), size = 4) +
      scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
      scale_x_discrete(labels = .wrapCellList) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
      labs(x = x_label, y = "Z-score (rank enrichment)",
           title = title_z, subtitle = subtitle, fill = "FDR < 0.05")
  }
  # The Z plot and the FDR plot are the same ranking twice over -- FDR is
  # a monotone function of Z here -- so only the FDR one is written.
  if (write_z_plot) {
    ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name,
                                       "_global_enrichment", file_suffix, ".pdf")),
           p_enrich, width = 8, height = 6)
  }

  # --- FDR bar plot ---
  if (use_palette) {
    ge_ord2 <- ge[order(-log10(ge$fdr)), ]
    ge_ord2[[id_col]] <- factor(ge_ord2[[id_col]], levels = ge_ord2[[id_col]])
    p_fdr <- ggplot(ge_ord2, aes(x = .data[[id_col]], y = -log10(fdr), fill = .data[[id_col]])) +
      geom_bar(stat = "identity") + coord_flip() + theme_bw() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.6) +
      geom_text(aes(label = fdr_stars), hjust = -0.3, size = 4) +
      scale_fill_manual(values = fill_palette, labels = .wrapCellList) +
      scale_x_discrete(labels = .wrapCellList) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
      labs(x = x_label, y = "-log10(FDR)",
           title = title_fdr, subtitle = subtitle, fill = x_label)
  } else {
    p_fdr <- ggplot(ge, aes(x = reorder(.data[[id_col]], -log10(fdr)),
                            y = -log10(fdr), fill = positive)) +
      geom_bar(stat = "identity") + coord_flip() + theme_bw() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.6) +
      geom_text(aes(label = fdr_stars), hjust = -0.3, size = 4) +
      scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
      scale_x_discrete(labels = .wrapCellList) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
      labs(x = x_label, y = "-log10(FDR)",
           title = title_fdr, subtitle = subtitle, fill = "FDR < 0.05")
  }
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name,
                                     "_global_enrichment", file_suffix, "_fdr.pdf")),
         p_fdr, width = 8, height = 6)

  ge
}
