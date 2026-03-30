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

addSpecificity <- function(gwas_mat, snps, high_thresh, mid_thresh, results_dir, cell_lineage = NULL) {
  if (!nrow(gwas_mat)) return(data.frame())
  
  specificity_summary_l2_detailed <- compute_specificity_summary(gwas_mat)
  specificity_basic <- data.frame(
    peak = rownames(gwas_mat),
    n_high = rowSums(gwas_mat > high_thresh),
    n_mid = rowSums(gwas_mat > mid_thresh & gwas_mat <= high_thresh),
    n_low = rowSums(gwas_mat <= mid_thresh)
  )
  
  high_cells <- apply(gwas_mat, 1, function(x) names(x)[x > high_thresh])
  specificity_basic$high_cells <- sapply(high_cells, paste, collapse = ", ")
  specificity_basic$multiple_high <- lengths(high_cells) > 1
  
  specificity_summary_l2 <- specificity_summary_l2_detailed %>%
    left_join(specificity_basic, by = "peak")

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
  } else {
    specificity_summary_l2$dominant_lineage <- NA_character_
    specificity_summary_l2$dominant_lineage_ratio <- NA_real_
    specificity_summary_l2$dominant_lineage_is_dominant <- NA
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
      maximum_cell, maximum_score, dominant_ratio, is_dominant,
      dominant_lineage, dominant_lineage_ratio, dominant_lineage_is_dominant,
      n_high, n_mid, n_low, high_cells, multiple_high,
      k_cells, top_k_cells, coherence
    )
  
  write.table(specificity_summary_l2, file.path(results_dir, "specificity_summary_l2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  specificity_summary_l2
}
