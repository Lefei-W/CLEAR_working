# Looping helpers to run compute_metrics_se (with quantile-norm, L2, covariance/correlation weighting) across multiple SummarizedExperiment objects.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Matrix)
  library(preprocessCore)
  library(corpcor)
  library(GenomicRanges)
})

  # Defaults (override at call-site as needed)
  DEFAULT_RATIO_THRESH        <- 2
  DEFAULT_SPEC_THRESHOLD_SQ   <- 0.5
  DEFAULT_L2_HIGH_THRESHOLD   <- NULL
  DEFAULT_L2_MID_THRESHOLD    <- NULL
  DEFAULT_COV_HIGH_QUANTILE   <- 0.9
  DEFAULT_COV_MID_QUANTILE    <- 0.5

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
compute_metrics_se <- function(se, assay_name = NULL, keep_top_prop = NULL) {
  stopifnot(inherits(se, "SummarizedExperiment"))
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

process_se_list_metrics <- function(se_list, assay_name = NULL, keep_top_prop = NULL) {
  lapply(se_list, function(se) compute_metrics_se(se, assay_name = assay_name, keep_top_prop = keep_top_prop, save_processed = FALSE, out_dir = file.path("data", "processed_se")))
}

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

overlap_se_with_snps <- function(se, snps, assay_name) {
  overlap_se_with_snps <- function(se, snps, assay_name, save_overlap = FALSE, tissue = NULL, trait = NULL, out_dir = file.path("data", "GWAS_overlapped_se")) {
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(inherits(snps, "GenomicRanges"))
  if (!assay_name %in% assayNames(se)) stop("Assay ", assay_name, " not found in SE")
  se <- ensure_peakids(se)
  ov <- findOverlaps(se, snps, ignore.strand = TRUE)
  if (!length(ov)) return(list(se_overlap = se[0, ], mat = matrix(0, 0, ncol(se))))
  hit_rows <- sort(unique(queryHits(ov)))
  se_sub <- se[hit_rows, , drop = FALSE]
  # add SNP IDs and signals per overlapped peak
  snp_ids <- mcols(snps)$names
  snp_signals <- mcols(snps)$signal
  snp_list <- tapply(snp_ids[subjectHits(ov)], queryHits(ov), paste, collapse = ",")
  locus_vec <- tapply(snp_signals[subjectHits(ov)], queryHits(ov), function(x) paste(unique(x), collapse = ", "))
  rowData(se_sub)$overlap_snps <- snp_list[as.character(seq_along(se))[hit_rows]]
  rowData(se_sub)$overlap_signals <- locus_vec[as.character(seq_along(se))[hit_rows]]
  mat <- SummarizedExperiment::assay(se_sub, assay_name)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  save_path <- NULL
  if (save_overlap) {
    if (is.null(tissue)) stop("tissue is required when save_overlap = TRUE")
    if (is.null(trait))  stop("trait is required when save_overlap = TRUE")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    n_ct <- ncol(se_sub)
    n_pk <- nrow(se_sub)
    fname <- paste0(tissue, "_", n_ct, "ct_", trait, "_", n_pk, "peaks_overlapped.rds")
    save_path <- file.path(out_dir, fname)
    saveRDS(se_sub, save_path)
  }
  list(se_overlap = se_sub, mat = mat, save_path = save_path)
}

# Pairwise loop: all SEs x all SNP trait sets; returns nested list
process_se_trait_pairs <- function(
    se_list,
    trait_list,
    metric_assays = c("qn_l2"),
    cov_metric_assays = c("cov_qn_l2"),
  dominance_ratio_thresh = DEFAULT_RATIO_THRESH,
  l2_high_threshold = DEFAULT_L2_HIGH_THRESHOLD,
  l2_mid_threshold  = DEFAULT_L2_MID_THRESHOLD,
  cov_high_quantile = DEFAULT_COV_HIGH_QUANTILE,
  cov_mid_quantile  = DEFAULT_COV_MID_QUANTILE,
  specificity_threshold_sq = DEFAULT_SPEC_THRESHOLD_SQ
) {
  if (is.null(names(se_list))) names(se_list) <- paste0("se", seq_along(se_list))
  if (is.null(names(trait_list))) names(trait_list) <- paste0("trait", seq_along(trait_list))

  lapply(names(se_list), function(se_nm) {
    se_obj <- se_list[[se_nm]]
    lapply(names(trait_list), function(trait_nm) {
      snps <- trait_list[[trait_nm]]

      metric_results <- lapply(metric_assays, function(metric_assay) {
        ov <- overlap_se_with_snps(se_obj, snps, assay_name = metric_assay)
        mat_l2 <- ov$mat
        if (nrow(mat_l2) == 0) {
          return(list(
            metric = metric_assay,
            overlap_n = 0,
            spec_basic = NULL,
            dominance  = NULL,
            spec_summary = NULL,
            spec_cov = NULL
          ))
        }

        spec_basic <- specificity_counts_l2(mat_l2, l2_high_thresh = l2_high_threshold, l2_mid_thresh = l2_mid_threshold)
        dominance  <- compute_dominance(mat_l2, ratio_thresh = dominance_ratio_thresh)
        spec_summary <- compute_specificity_summary(mat_l2, ratio_thresh = dominance_ratio_thresh, threshold_sq = specificity_threshold_sq)

        spec_cov <- NULL
        hits_cov <- cov_metric_assays[cov_metric_assays %in% assayNames(ov$se_overlap)]
        if (length(hits_cov)) {
          spec_cov <- lapply(hits_cov, function(cov_assay) {
            mat_cov <- SummarizedExperiment::assay(ov$se_overlap, cov_assay)
            if (inherits(mat_cov, "Matrix")) mat_cov <- as.matrix(mat_cov)
            list(
              cov_assay = cov_assay,
              spec_quantile = specificity_counts_quantile(mat_cov, cov_high_quantile = cov_high_quantile, cov_mid_quantile = cov_mid_quantile)
            )
          })
        }

        list(
          metric = metric_assay,
          overlap_n = nrow(ov$se_overlap),
          spec_basic = spec_basic,
          dominance  = dominance,
          spec_summary = spec_summary,
          spec_cov = spec_cov
        )
      })

      list(se_name = se_nm, trait = trait_nm, metrics = metric_results)
    })
  })
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

# -------- IO helpers for processed and overlapped SEs --------
save_processed_se <- function(se, tissue, out_dir = file.path("data", "processed_se")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_ct <- ncol(se)
  n_pk <- nrow(se)
  fname <- paste0(tissue, "_", n_ct, "ct_", n_pk, "peaks_processed.rds")
  path <- file.path(out_dir, fname)
  saveRDS(se, path)
  path
}

save_overlap_se <- function(se_overlap, tissue, trait, out_dir = file.path("data", "GWAS_overlapped_se")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_ct <- ncol(se_overlap)
  n_pk <- nrow(se_overlap)
  fname <- paste0(tissue, "_", n_ct, "ct_", trait, "_", n_pk, "peaks_overlapped.rds")
  path <- file.path(out_dir, fname)
  saveRDS(se_overlap, path)
  path
}



