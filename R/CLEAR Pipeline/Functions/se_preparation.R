# ============================================================
# se_preparation.R
# SE loading, coordinate repair, L2 normalisation, and metrics
# ============================================================

prepare_clear_se <- function(se, genome = "hg38") {
  rd <- as.data.frame(rowData(se))
  
  has_rowdata_coords <- all(c("seqnames", "start", "end") %in% names(rd))
  
  if (has_rowdata_coords) {
    gr <- GRanges(
      seqnames = rd$seqnames,
      ranges   = IRanges::IRanges(start = as.integer(rd$start), end = as.integer(rd$end))
    )
    
    extra_cols <- setdiff(names(rd), c("seqnames", "start", "end", "width", "strand"))
    if (length(extra_cols)) mcols(gr) <- rd[, extra_cols, drop = FALSE]
  } else if (inherits(rowRanges(se), "GenomicRanges")) {
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
  
  se <- SummarizedExperiment(
    assays    = as.list(assays(se)),
    rowRanges = gr,
    colData   = colData(se)
  )
  
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

compute_l2_weights <- function(se) {
  raw_mat <- assay(se, "raw")
  if (inherits(raw_mat, "Matrix")) raw_mat <- as.matrix(raw_mat)
  rowMeans(raw_mat)
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
