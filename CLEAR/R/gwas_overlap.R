# ============================================================
# gwas_overlap.R
# Load CLEAR inputs, GWAS-peak overlap, input summary stats.
# (Per-peak specificity tables live in specificity.R; signal-
# level summary in signal_summary.R; dominant-peak coherence in
# dominant_peaks.R.)
# ============================================================

# ------------------------------------------------------------
# loadCLEARInputs
# Read SE and SNP RDS files and apply prepareSE().
# ------------------------------------------------------------
loadCLEARInputs <- function(se_path, trait_path) {
  se   <- readRDS(se_path)
  snps <- readRDS(trait_path)
  se   <- prepareSE(se)
  list(se = se, snps = snps)
}

# ------------------------------------------------------------
# addGWASOverlap
# Overlap variants with peaks. Returns the GWAS-overlapping SE
# subset (with CCVs / signal annotation), the L2 GWAS matrix,
# and the raw findOverlaps Hits object.
# ------------------------------------------------------------
addGWASOverlap <- function(se, snps) {
  hits_se <- findOverlaps(snps, rowRanges(se))
  gwas_idx <- unique(subjectHits(hits_se))

  if (!length(gwas_idx)) {
    se_gwas <- se[0, , drop = FALSE]
  } else {
    se_gwas <- se[gwas_idx, , drop = FALSE]

    snp_name_col <- if ("names" %in% colnames(mcols(snps))) mcols(snps)$names else names(snps)
    signal_col   <- if ("signal" %in% colnames(mcols(snps))) mcols(snps)$signal else rep(NA_character_, length(snps))

    ann <- data.frame(
      peak_idx = subjectHits(hits_se),
      CCVs     = snp_name_col[queryHits(hits_se)],
      signal   = signal_col[queryHits(hits_se)],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::group_by(peak_idx) %>%
      dplyr::summarise(
        CCVs   = paste(unique(na.omit(CCVs)),   collapse = ","),
        signal = paste(unique(na.omit(signal)), collapse = ","),
        .groups = "drop"
      )

    rd <- as.data.frame(rowData(se_gwas))
    if (!"CCVs"   %in% names(rd)) rd$CCVs   <- ""
    if (!"signal" %in% names(rd)) rd$signal <- ""

    pos <- match(ann$peak_idx, gwas_idx)
    rd$CCVs[pos]   <- ann$CCVs
    rd$signal[pos] <- ann$signal
    rowData(se_gwas) <- S4Vectors::DataFrame(rd, row.names = rownames(se_gwas))
  }

  list(
    hits_se  = hits_se,
    gwas_idx = gwas_idx,
    se_gwas  = se_gwas,
    gwas_mat = assay(se_gwas, "raw_l2")
  )
}

# ------------------------------------------------------------
# getInputSummary
# Cell-type and lineage-level raw-signal summaries; writes:
#   raw_signal_per_celltype.txt
#   raw_signal_per_lineage.txt
#   input_summary.txt
# ------------------------------------------------------------
getInputSummary <- function(se_name, trait_name, se_path, trait_path,
                            se, snps, overlap, results_dir,
                            cell_lineage = NULL) {
  n_unique_signals <- if ("signal" %in% colnames(mcols(snps)))
    length(unique(na.omit(mcols(snps)$signal))) else NA_integer_
  if (length(overlap$gwas_idx) > 0 && "signal" %in% colnames(mcols(snps))) {
    overlapping_signals <- mcols(snps)$signal[queryHits(overlap$hits_se)]
    n_unique_signals_in_peaks <- length(unique(na.omit(overlapping_signals)))
  } else {
    n_unique_signals_in_peaks <- 0
  }

  n_snp_peak_overlaps    <- length(overlap$hits_se)
  n_unique_snps_in_peaks <- if (length(overlap$hits_se)) length(unique(queryHits(overlap$hits_se))) else 0L

  raw_mat <- assay(se, "raw")
  raw_col_sums <- colSums(raw_mat)

  raw_colsum_summary <- data.frame(
    celltype = names(raw_col_sums),
    total_raw_signal = as.numeric(raw_col_sums),
    stringsAsFactors = FALSE
  )
  raw_colsum_summary <- raw_colsum_summary[order(-raw_colsum_summary$total_raw_signal), ]
  write.table(raw_colsum_summary,
              file.path(results_dir, "raw_signal_per_celltype.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (!is.null(cell_lineage)) {
    lineage_vec <- cell_lineage[colnames(raw_mat)]
    missing_lin <- is.na(lineage_vec) | !nzchar(lineage_vec)
    if (any(missing_lin)) lineage_vec[missing_lin] <- colnames(raw_mat)[missing_lin]
    raw_lineage_sums <- tapply(as.numeric(raw_col_sums), lineage_vec, sum)

    raw_lineage_summary <- data.frame(
      lineage = names(raw_lineage_sums),
      total_raw_signal = as.numeric(raw_lineage_sums),
      stringsAsFactors = FALSE
    )
    raw_lineage_summary <- raw_lineage_summary[order(-raw_lineage_summary$total_raw_signal), ]
    write.table(raw_lineage_summary,
                file.path(results_dir, "raw_signal_per_lineage.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    n_lineages           <- length(raw_lineage_sums)
    mean_total_raw_lin   <- mean(raw_lineage_sums)
    var_total_raw_lin    <- var(as.numeric(raw_lineage_sums))
    max_signal_lineage   <- names(which.max(raw_lineage_sums))
    max_signal_lin_value <- max(raw_lineage_sums)
    min_signal_lineage   <- names(which.min(raw_lineage_sums))
    min_signal_lin_value <- min(raw_lineage_sums)
  } else {
    n_lineages           <- NA_integer_
    mean_total_raw_lin   <- NA_real_
    var_total_raw_lin    <- NA_real_
    max_signal_lineage   <- NA_character_
    max_signal_lin_value <- NA_real_
    min_signal_lineage   <- NA_character_
    min_signal_lin_value <- NA_real_
  }

  input_summary <- data.frame(
    se_name = se_name,
    trait_name = trait_name,
    GWAS_trait = basename(trait_path),
    snATAC_dataset = basename(se_path),
    n_peaks = nrow(se),
    n_celltypes = ncol(se),
    n_lineages = n_lineages,
    n_snps = length(snps),
    n_unique_signals = n_unique_signals,
    n_gwas_peaks = nrow(overlap$gwas_mat),
    n_snp_peak_overlaps = n_snp_peak_overlaps,
    n_unique_snps_in_peaks = n_unique_snps_in_peaks,
    n_unique_signals_in_peaks = n_unique_signals_in_peaks,
    mean_total_raw_signal = mean(raw_col_sums),
    var_total_raw_signal  = var(as.numeric(raw_col_sums)),
    max_signal_celltype = names(which.max(raw_col_sums)),
    max_signal_value    = max(raw_col_sums),
    min_signal_celltype = names(which.min(raw_col_sums)),
    min_signal_value    = min(raw_col_sums),
    mean_total_raw_signal_lineage = mean_total_raw_lin,
    var_total_raw_signal_lineage  = var_total_raw_lin,
    max_signal_lineage = max_signal_lineage,
    max_signal_lineage_value = max_signal_lin_value,
    min_signal_lineage = min_signal_lineage,
    min_signal_lineage_value = min_signal_lin_value,
    stringsAsFactors = FALSE
  )
  write.table(input_summary, file.path(results_dir, "input_summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  input_summary
}