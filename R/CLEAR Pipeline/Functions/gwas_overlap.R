# ============================================================
# gwas_overlap.R
# Load inputs, GWAS-peak overlap, input summary, specificity
# ============================================================

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

addSpecificity <- function(gwas_mat, snps, high_thresh, mid_thresh, results_dir,
                           cell_lineage = NULL, se_gwas = NULL, gtf_path = NULL) {
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
        top_lineage_ratio = ratio,
        top_lineage_is_dominant = is_dominant
      )
    specificity_summary_l2 <- specificity_summary_l2 %>%
      left_join(dom_lineage, by = "peak")

    # --- Lineage-level high / mid counts and names ---
    # Sum L2^2 within each lineage (same as compute_dominance_lineage)
    lineage_vec <- cell_lineage[colnames(gwas_mat)]
    lineage_levels <- unique(as.character(lineage_vec))
    mat_sq <- gwas_mat^2
    lineage_mat <- sapply(lineage_levels, function(lg) {
      rowSums(mat_sq[, lineage_vec == lg, drop = FALSE])
    })
    lineage_mat <- as.matrix(lineage_mat)
    rownames(lineage_mat) <- rownames(gwas_mat)

    # Thresholds on lineage_mat (L2^2-summed): a lineage is "high" if its
    # summed L2^2 > high_thresh^2 (equivalent to a single cell at high_thresh)
    lin_high_thresh <- high_thresh^2
    lin_mid_thresh  <- mid_thresh^2

    specificity_summary_l2$n_high_lineage <- rowSums(lineage_mat > lin_high_thresh)
    specificity_summary_l2$n_mid_lineage  <- rowSums(lineage_mat > lin_mid_thresh & lineage_mat <= lin_high_thresh)
    specificity_summary_l2$n_low_lineage  <- rowSums(lineage_mat <= lin_mid_thresh)

    high_lin <- apply(lineage_mat, 1, function(x) names(x)[x > lin_high_thresh])
    mid_lin  <- apply(lineage_mat, 1, function(x) names(x)[x > lin_mid_thresh & x <= lin_high_thresh])
    specificity_summary_l2$high_lineages <- sapply(high_lin, paste, collapse = ", ")
    specificity_summary_l2$mid_lineages  <- sapply(mid_lin,  paste, collapse = ", ")

    # Top lineage name and score (argmax of lineage_mat, always populated)
    specificity_summary_l2$top_lineage <- colnames(lineage_mat)[max.col(lineage_mat, ties.method = "first")]
    specificity_summary_l2$top_lineage_score <- apply(lineage_mat, 1, max)
  } else {
    specificity_summary_l2$top_lineage_ratio <- NA_real_
    specificity_summary_l2$top_lineage_is_dominant <- NA
    specificity_summary_l2$n_high_lineage <- NA_integer_
    specificity_summary_l2$n_mid_lineage  <- NA_integer_
    specificity_summary_l2$n_low_lineage  <- NA_integer_
    specificity_summary_l2$high_lineages  <- NA_character_
    specificity_summary_l2$mid_lineages   <- NA_character_
    specificity_summary_l2$top_lineage    <- NA_character_
    specificity_summary_l2$top_lineage_score <- NA_real_
  }

  # --- Top cell type (argmax L2, same as maximum_cell — kept as alias) ---
  # (maximum_cell already comes from compute_specificity_summary)

  # --- Weight: mean raw accessibility across cell types ---
  if (!is.null(se_gwas) && "raw" %in% assayNames(se_gwas)) {
    raw_gwas <- assay(se_gwas, "raw")
    if (inherits(raw_gwas, "Matrix")) raw_gwas <- as.matrix(raw_gwas)
    specificity_summary_l2$weight <- rowMeans(raw_gwas[rownames(gwas_mat), , drop = FALSE])
  } else {
    # Fallback: cannot compute without raw counts
    specificity_summary_l2$weight <- NA_real_
  }

  # --- Distance to nearest TSS and TSS bin ---
  if (!is.null(se_gwas) && !is.null(gtf_path)) {
    gtf <- rtracklayer::import(gtf_path)
    tss_gr <- get_tss(gtf)
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

  specificity_summary_l2 <- add_gwas_peak_annotation(specificity_summary_l2, snps)
  specificity_summary_l2$locus <- ifelse(
    is.na(specificity_summary_l2$gwas_locus),
    "",
    as.character(specificity_summary_l2$gwas_locus)
  )
  
  specificity_summary_l2 <- specificity_summary_l2 %>%
    select(
      locus, peak, gwas_snps,
      # Cell-type level
      maximum_cell, maximum_score, dominant_ratio, is_dominant,
      n_high, n_mid, n_low, high_cells, mid_cells, multiple_high,
      # Lineage level
      top_lineage, top_lineage_score, top_lineage_ratio, top_lineage_is_dominant,
      n_high_lineage, n_mid_lineage, n_low_lineage, high_lineages, mid_lineages,
      # Genomic context
      dist_to_tss, tss_bin,
      # Weight and detailed metrics
      weight, k_cells, top_k_cells, coherence
    )
  
  write.table(specificity_summary_l2, file.path(results_dir, "specificity_summary_l2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  specificity_summary_l2
}
