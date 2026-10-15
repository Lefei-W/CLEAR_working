# ============================================================
# se_preparation.R
# SE loading, coordinate repair, L2 normalisation and metrics
# ============================================================

# ------------------------------------------------------------
# prepareSE
# Repair / standardise an input SummarizedExperiment so that
# downstream CLEAR functions can use it: ensures rowRanges with
# UCSC seqlevels, peakid rownames, and a consistent layout.
# ------------------------------------------------------------
prepareSE <- function(se, genome = "hg38") {
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
        if (!nm %in% names(mcols(gr))) mcols(gr)[[nm]] <- rd[[nm]]
      }
    }
  } else {
    stop("Input SE has no usable genomic coordinates. Provide either rowData seqnames/start/end or a valid rowRanges(se).")
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

# Internal: safe row-wise L2 normalisation (returns x unchanged if sum(x^2)==0)
.safeL2 <- function(x) {
  s <- sum(x^2)
  if (s == 0) x else x / sqrt(s)
}

# Internal: column-wise descending rank
.rankDesc <- function(M) {
  out <- vapply(seq_len(ncol(M)), function(j) rank(-M[, j], ties.method = "average"),
                numeric(nrow(M)))
  matrix(out, nrow(M), ncol(M), dimnames = dimnames(M))
}

# ------------------------------------------------------------
# addSEMetrics
# Computes assays stored in the SE:
#   raw         - original counts (renamed for clarity)
#   raw_l2      - L2-normalised (unit-norm rows, sum of s^2 = 1)
#   raw_l2_rank - column-wise descending ranks of raw_l2
#   norm        - accessibility after optional per-cell/depth normalisation
# Defaults leave the input scale unchanged. raw_l2 is computed from norm.
# ------------------------------------------------------------
addSEMetrics <- function(se, assay_name = NULL, keep_top_prop = NULL,
                         per_cell = c("no", "auto", "yes"),
                         depth = c("none", "colsum", "quantile"), n_cells = NULL) {
  per_cell <- match.arg(per_cell)
  depth <- match.arg(depth)
  if (is.null(assay_name)) assay_name <- assayNames(se)[1]

  mat <- assay(se, assay_name)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  if (!nrow(mat) || !ncol(mat) || any(!is.finite(mat)) || any(mat < 0)) {
    stop("Input accessibility must be a non-empty, finite, non-negative matrix.")
  }

  rn <- rownames(se); if (is.null(rn)) rn <- paste0("peak_", seq_len(nrow(mat)))
  cn <- colnames(se); if (is.null(cn)) cn <- paste0("ct_", seq_len(ncol(mat)))
  dimnames(mat) <- list(rn, cn)
  dimnames(se) <- list(rn, cn)

  row_sums <- rowSums(mat)
  if (!is.null(keep_top_prop)) {
    if (length(keep_top_prop) != 1L || !is.finite(keep_top_prop) ||
        keep_top_prop <= 0 || keep_top_prop > 1) stop("keep_top_prop must be in (0, 1].")
    thr <- stats::quantile(row_sums, probs = 1 - keep_top_prop, na.rm = TRUE)
    keep <- row_sums > thr
    se <- se[keep, , drop = FALSE]
    mat <- mat[keep, , drop = FALSE]
    rn <- rn[keep]
    dimnames(mat) <- list(rn, cn)
    if (!nrow(mat)) stop("keep_top_prop removed every peak.")
  }

  normalised <- .normaliseAccessibility(mat, se, per_cell, depth, n_cells)
  norm <- normalised$norm
  denom <- sqrt(rowSums(norm^2))
  denom[denom == 0] <- 1
  raw_l2 <- norm / denom

  assays(se)[["raw"]]         <- mat
  assays(se)[["norm"]]        <- norm
  assays(se)[["raw_l2"]]      <- raw_l2
  assays(se)[["raw_l2_rank"]] <- .rankDesc(raw_l2)

  metadata(se)$CLEAR <- list(
    source_assay  = assay_name,
    keep_top_prop = keep_top_prop,
    n_peaks       = nrow(se),
    n_celltypes   = ncol(se),
    use_l2_only   = TRUE,
    normalisation = normalised$metadata
  )
  se
}