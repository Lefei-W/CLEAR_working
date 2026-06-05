# ============================================================
# specificity.R
# Variant-peak specificity tables (cell + lineage resolution)
# and their comparison.
# ============================================================

# ------------------------------------------------------------
# addSpecificity
# Per variant-peak summary on the L2 (cell-type) matrix:
# top cell type, dominance ratio, tier counts (high/mid/low),
# weighted accessibility, TSS distance, lineage-level mirror
# fields, and GWAS locus / SNP annotation.
# Writes specificity_summary_l2.txt.
# ------------------------------------------------------------
addSpecificity <- function(gwas_mat, snps, high_thresh, mid_thresh, results_dir,
                           cell_lineage = NULL, se_gwas = NULL, gtf_path = NULL) {
  if (!nrow(gwas_mat)) return(data.frame())

  specificity_summary_l2_detailed <- addSpecificitySummary(gwas_mat)

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
    dom_lineage <- addDominanceLineage(gwas_mat, cell_lineage, ratio_thresh = DEFAULT_RATIO_THRESH) %>%
      dplyr::select(
        peak,
        top_lineage_ratio = ratio,
        top_lineage_is_dominant = is_dominant
      )
    specificity_summary_l2 <- specificity_summary_l2 %>%
      left_join(dom_lineage, by = "peak")

    # --- Lineage-level high / mid counts and names ---
    lineage_mat <- getLineageMatrix(gwas_mat, cell_lineage)

    # Same thresholds as cell-type level (lineage_mat is on L2 scale too)
    lin_high_thresh <- high_thresh
    lin_mid_thresh  <- mid_thresh

    specificity_summary_l2$n_high_lineage <- rowSums(lineage_mat > lin_high_thresh)
    specificity_summary_l2$n_mid_lineage  <- rowSums(lineage_mat > lin_mid_thresh & lineage_mat <= lin_high_thresh)
    specificity_summary_l2$n_low_lineage  <- rowSums(lineage_mat <= lin_mid_thresh)

    high_lin <- apply(lineage_mat, 1, function(x) names(x)[x > lin_high_thresh])
    mid_lin  <- apply(lineage_mat, 1, function(x) names(x)[x > lin_mid_thresh & x <= lin_high_thresh])
    specificity_summary_l2$high_lineages <- sapply(high_lin, paste, collapse = ", ")
    specificity_summary_l2$mid_lineages  <- sapply(mid_lin,  paste, collapse = ", ")

    specificity_summary_l2$top_lineage <- colnames(lineage_mat)[max.col(lineage_mat, ties.method = "first")]
    specificity_summary_l2$top_lineage_score <- apply(lineage_mat, 1, max)
  } else {
    specificity_summary_l2$top_lineage_ratio        <- NA_real_
    specificity_summary_l2$top_lineage_is_dominant  <- NA
    specificity_summary_l2$n_high_lineage           <- NA_integer_
    specificity_summary_l2$n_mid_lineage            <- NA_integer_
    specificity_summary_l2$n_low_lineage            <- NA_integer_
    specificity_summary_l2$high_lineages            <- NA_character_
    specificity_summary_l2$mid_lineages             <- NA_character_
    specificity_summary_l2$top_lineage              <- NA_character_
    specificity_summary_l2$top_lineage_score        <- NA_real_
  }

  # --- Weights: mean and max raw accessibility ---
  # weight_mean reflects overall accessibility; weight_max separates
  # specifically accessible peaks from uniformly moderate ones.
  if (!is.null(se_gwas) && "raw" %in% assayNames(se_gwas)) {
    raw_gwas <- assay(se_gwas, "raw")
    if (inherits(raw_gwas, "Matrix")) raw_gwas <- as.matrix(raw_gwas)
    raw_sub <- raw_gwas[rownames(gwas_mat), , drop = FALSE]
    specificity_summary_l2$weight_mean <- rowMeans(raw_sub)
    specificity_summary_l2$weight_max  <- apply(raw_sub, 1, max)
  } else {
    specificity_summary_l2$weight_mean <- NA_real_
    specificity_summary_l2$weight_max  <- NA_real_
  }

  # --- Distance to nearest TSS and TSS bin ---
  if (!is.null(se_gwas) && !is.null(gtf_path)) {
    gtf <- rtracklayer::import(gtf_path)
    tss_gr <- .getTSS(gtf)
    peaks_gr <- rowRanges(se_gwas)
    nearest_hits <- distanceToNearest(peaks_gr, tss_gr)
    specificity_summary_l2$dist_to_tss <- mcols(nearest_hits)$distance[match(
      seq_len(nrow(gwas_mat)), queryHits(nearest_hits)
    )]
    specificity_summary_l2$tss_bin <- as.character(cut(
      specificity_summary_l2$dist_to_tss,
      breaks = c(0, 1000, 5000, 50000, 200000, Inf),
      labels = c("core_prom", "prox_prom", "near_reg", "distal", "long_range"),
      include.lowest = TRUE
    ))
  } else {
    specificity_summary_l2$dist_to_tss <- NA_real_
    specificity_summary_l2$tss_bin <- NA_character_
  }

  specificity_summary_l2 <- .addGWASPeakAnnotation(specificity_summary_l2, snps)
  specificity_summary_l2$locus <- ifelse(
    is.na(specificity_summary_l2$gwas_locus),
    "",
    as.character(specificity_summary_l2$gwas_locus)
  )

  specificity_summary_l2 <- specificity_summary_l2 %>%
    select(
      locus, peak, gwas_snps,
      maximum_cell, maximum_score, dominant_ratio, is_dominant,
      n_high, n_mid, n_low, high_cells, mid_cells, multiple_high,
      top_lineage, top_lineage_score, top_lineage_ratio, top_lineage_is_dominant,
      n_high_lineage, n_mid_lineage, n_low_lineage, high_lineages, mid_lineages,
      dist_to_tss, tss_bin,
      weight_mean, weight_max, k_cells, top_k_cells, coherence
    )

  write.table(specificity_summary_l2, file.path(results_dir, "specificity_summary_l2.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  specificity_summary_l2
}

# ------------------------------------------------------------
# addLineageSpecificity
# Variant-peak level specificity table on the LINEAGE-level L2
# matrix. Mirrors addSpecificity() but with lineages as the unit
# of specificity. Each peak's C-vector is collapsed to a G-vector
# via sqrt(sum_{c in g} s_c^2), preserving unit norm. Thresholds:
#   high : 1/sqrt(2)            (one lineage > 50% squared mass)
#   mid  : 1/sqrt(G)            (G = number of lineages)
# Writes lineage_specificity_summary_l2.txt.
# ------------------------------------------------------------
addLineageSpecificity <- function(gwas_mat, snps, results_dir,
                                  cell_lineage, se_gwas = NULL, gtf_path = NULL,
                                  ratio_thresh = DEFAULT_RATIO_THRESH) {
  if (!nrow(gwas_mat)) return(data.frame())
  if (is.null(cell_lineage)) {
    message("  Skipping lineage specificity: cell_lineage not provided")
    return(data.frame())
  }

  lineage_mat <- getLineageMatrix(gwas_mat, cell_lineage)
  G <- ncol(lineage_mat)
  high_thresh_lin <- 1 / sqrt(2)
  mid_thresh_lin  <- 1 / sqrt(max(G, 2))

  spec_lin_detailed <- addSpecificitySummary(lineage_mat, ratio_thresh = ratio_thresh)
  # Rename to lineage-level vocabulary
  names(spec_lin_detailed)[names(spec_lin_detailed) == "maximum_cell"]   <- "maximum_lineage"
  names(spec_lin_detailed)[names(spec_lin_detailed) == "maximum_score"]  <- "maximum_score_lineage"
  names(spec_lin_detailed)[names(spec_lin_detailed) == "dominant_ratio"] <- "dominant_ratio_lineage"
  names(spec_lin_detailed)[names(spec_lin_detailed) == "is_dominant"]    <- "is_dominant_lineage"
  names(spec_lin_detailed)[names(spec_lin_detailed) == "k_cells"]        <- "k_lineages"
  names(spec_lin_detailed)[names(spec_lin_detailed) == "top_k_cells"]    <- "top_k_lineages"
  names(spec_lin_detailed)[names(spec_lin_detailed) == "coherence"]      <- "coherence_lineage"

  spec_basic <- data.frame(
    peak = rownames(lineage_mat),
    n_high_lineage = rowSums(lineage_mat > high_thresh_lin),
    n_mid_lineage  = rowSums(lineage_mat > mid_thresh_lin & lineage_mat <= high_thresh_lin),
    n_low_lineage  = rowSums(lineage_mat <= mid_thresh_lin),
    stringsAsFactors = FALSE
  )
  high_lin <- apply(lineage_mat, 1, function(x) names(x)[x > high_thresh_lin])
  mid_lin  <- apply(lineage_mat, 1, function(x) names(x)[x > mid_thresh_lin & x <= high_thresh_lin])
  spec_basic$high_lineages         <- sapply(high_lin, paste, collapse = ", ")
  spec_basic$mid_lineages          <- sapply(mid_lin,  paste, collapse = ", ")
  spec_basic$multiple_high_lineage <- lengths(high_lin) > 1

  lineage_summary <- spec_lin_detailed %>% dplyr::left_join(spec_basic, by = "peak")

  if (!is.null(se_gwas) && "raw" %in% assayNames(se_gwas)) {
    raw_gwas <- assay(se_gwas, "raw")
    if (inherits(raw_gwas, "Matrix")) raw_gwas <- as.matrix(raw_gwas)
    raw_sub <- raw_gwas[rownames(lineage_mat), , drop = FALSE]
    lineage_summary$weight_mean <- rowMeans(raw_sub)
    lineage_summary$weight_max  <- apply(raw_sub, 1, max)
  } else {
    lineage_summary$weight_mean <- NA_real_
    lineage_summary$weight_max  <- NA_real_
  }

  if (!is.null(se_gwas) && !is.null(gtf_path) && file.exists(gtf_path)) {
    gtf <- rtracklayer::import(gtf_path)
    tss_gr <- .getTSS(gtf)
    peaks_gr <- rowRanges(se_gwas)
    nearest_hits <- distanceToNearest(peaks_gr, tss_gr)
    lineage_summary$dist_to_tss <- mcols(nearest_hits)$distance[match(
      seq_len(nrow(lineage_mat)), queryHits(nearest_hits)
    )]
    lineage_summary$tss_bin <- as.character(cut(
      lineage_summary$dist_to_tss,
      breaks = c(0, 1000, 5000, 50000, 200000, Inf),
      labels = c("core_prom", "prox_prom", "near_reg", "distal", "long_range"),
      include.lowest = TRUE
    ))
  } else {
    lineage_summary$dist_to_tss <- NA_real_
    lineage_summary$tss_bin <- NA_character_
  }

  lineage_summary <- .addGWASPeakAnnotation(lineage_summary, snps)
  lineage_summary$locus <- ifelse(
    is.na(lineage_summary$gwas_locus),
    "",
    as.character(lineage_summary$gwas_locus)
  )

  lineage_summary <- lineage_summary %>%
    dplyr::select(
      locus, peak, gwas_snps,
      maximum_lineage, maximum_score_lineage,
      dominant_ratio_lineage, is_dominant_lineage,
      n_high_lineage, n_mid_lineage, n_low_lineage,
      high_lineages, mid_lineages, multiple_high_lineage,
      dist_to_tss, tss_bin,
      weight_mean, weight_max,
      k_lineages, top_k_lineages, coherence_lineage
    )

  write.table(lineage_summary,
              file.path(results_dir, "lineage_specificity_summary_l2.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  lineage_summary
}

# ------------------------------------------------------------
# getSpecificityComparison
# Compare per-peak dominance calls between cell-type and lineage
# resolutions:
#   shared       : dominant at both resolutions
#   cell_only    : dominant only at cell-type resolution
#   lineage_only : dominant only at lineage resolution
#   neither      : not dominant at either resolution
# Writes specificity_comparison_cell_vs_lineage.txt.
# ------------------------------------------------------------
getSpecificityComparison <- function(specificity_summary_l2,
                                     lineage_specificity_summary_l2,
                                     results_dir) {
  if (!nrow(specificity_summary_l2) || !nrow(lineage_specificity_summary_l2)) {
    return(data.frame())
  }

  cell_keep <- intersect(
    c("locus", "peak", "gwas_snps",
      "maximum_cell", "maximum_score", "dominant_ratio", "is_dominant",
      "n_high", "n_mid", "n_low", "high_cells", "mid_cells", "multiple_high",
      "dist_to_tss", "tss_bin", "weight_mean", "weight_max",
      "k_cells", "top_k_cells", "coherence"),
    names(specificity_summary_l2)
  )
  lin_keep <- intersect(
    c("peak",
      "maximum_lineage", "maximum_score_lineage",
      "dominant_ratio_lineage", "is_dominant_lineage",
      "n_high_lineage", "n_mid_lineage", "n_low_lineage",
      "high_lineages", "mid_lineages", "multiple_high_lineage",
      "k_lineages", "top_k_lineages", "coherence_lineage"),
    names(lineage_specificity_summary_l2)
  )

  cmp <- dplyr::full_join(
    specificity_summary_l2[, cell_keep, drop = FALSE],
    lineage_specificity_summary_l2[, lin_keep, drop = FALSE],
    by = "peak"
  )

  is_dom_cell <- .isTRUEvec(cmp$is_dominant)
  is_dom_lin  <- .isTRUEvec(cmp$is_dominant_lineage)

  cmp$category <- dplyr::case_when(
     is_dom_cell &  is_dom_lin ~ "shared",
     is_dom_cell & !is_dom_lin ~ "cell_only",
    !is_dom_cell &  is_dom_lin ~ "lineage_only",
    TRUE ~ "neither"
  )

  front <- intersect(c("category", "locus", "peak", "gwas_snps",
                       "dist_to_tss", "tss_bin",
                       "weight_mean", "weight_max"),
                     names(cmp))
  cmp <- cmp[, c(front, setdiff(names(cmp), front)), drop = FALSE]

  write.table(cmp,
              file.path(results_dir, "specificity_comparison_cell_vs_lineage.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  tab <- table(cmp$category, useNA = "ifany")
  message("  Specificity comparison (cell vs lineage):")
  for (nm in names(tab)) message("    ", nm, ": ", tab[[nm]])

  cmp
}

# Internal: coerce mixed logical/character/NA vector to logical, NA -> FALSE
.isTRUEvec <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  res <- suppressWarnings(as.logical(x))
  !is.na(res) & res
}
