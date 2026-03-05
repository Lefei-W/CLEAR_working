safe_l2 <- function(x) { # on the rows (peaks)
  s <- sum(x^2)
  if (s == 0) x else x / sqrt(s)  # l2 euclidean normalization avoiding row sum of 0 (happen with the Human-scATAC-Corpus)
}

add_rank_desc <- function(M) {
  # rank 1 = highest value (descending)
  apply(M, 2, function(x) rank(-x, ties.method = "average"))
}

compute_metrics_se <- function(se, assay_name = NULL, keep_top_prop = NULL) {

  if (is.null(assay_name)) assay_name <- assayNames(se)[1]

  mat <- assay(se, assay_name)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  rn <- rownames(se); if (is.null(rn)) rn <- paste0("peak_", seq_len(nrow(mat)))
  cn <- colnames(se); if (is.null(cn)) cn <- paste0("ct_",   seq_len(ncol(mat)))
  dimnames(mat) <- list(rn, cn)

  # optional filter: keep top proportion by rowSum
  row_sums <- rowSums(mat)
  if (!is.null(keep_top_prop)) {
    thr  <- stats::quantile(row_sums, probs = 1 - keep_top_prop, na.rm = TRUE)
    keep <- row_sums > thr
    se  <- se[keep, , drop = FALSE]
    mat <- mat[keep, , drop = FALSE]
    rn  <- rn[keep]
    row_sums <- row_sums[keep]
    dimnames(mat) <- list(rn, cn)
  }

  # transforms
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

  # save assays (original assay kept as-is)
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