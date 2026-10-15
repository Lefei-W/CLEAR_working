SPEC_HIGH_THRESH <- 1 / sqrt(2)
DEFAULT_RATIO_THRESH <- 2

# Pooling and L2 division can move an exact boundary by a few machine epsilons.
.aboveThreshold <- function(x, threshold) {
  x > threshold + 8 * .Machine$double.eps * max(1, abs(threshold))
}

.belowThreshold <- function(x, threshold) {
  x < threshold - 8 * .Machine$double.eps * max(1, abs(threshold))
}

.specMatrix <- function(mat) {
  mat <- as.matrix(mat)
  if (!is.numeric(mat) || !ncol(mat) || any(!is.finite(mat)) || any(mat < 0)) {
    stop("Specificity must be a finite, non-negative matrix with at least one column.")
  }
  if (is.null(colnames(mat))) colnames(mat) <- paste0("ct", seq_len(ncol(mat)))
  if (is.null(rownames(mat)) && nrow(mat)) rownames(mat) <- paste0("peak_", seq_len(nrow(mat)))
  mat
}

getSpecTiers <- function(spec_mat) {
  spec_mat <- .specMatrix(spec_mat)
  mass <- rowSums(spec_mat^2)
  if (any(mass > 0 & abs(mass - 1) > 1e-8)) {
    stop("Specificity rows must have unit L2 norm (or be all zero).")
  }
  high <- .aboveThreshold(spec_mat, SPEC_HIGH_THRESH)
  low <- .belowThreshold(spec_mat, 1 / sqrt(max(ncol(spec_mat), 2L)))
  list(high = high, moderate = !high & !low, low = low)
}

.rowTop2 <- function(mat) {
  n <- nrow(mat)
  top_idx <- max.col(mat, ties.method = "first")
  top <- mat[cbind(seq_len(n), top_idx)]
  second_idx <- rep(NA_integer_, n)
  second <- rep(NA_real_, n)
  if (ncol(mat) > 1L) {
    other <- mat
    other[cbind(seq_len(n), top_idx)] <- -Inf
    second_idx <- max.col(other, ties.method = "first")
    second <- other[cbind(seq_len(n), second_idx)]
  }
  ratio <- rep(NA_real_, n)
  comparable <- !is.na(second) & top > 0
  ratio[comparable] <- top[comparable] / second[comparable]
  list(top_idx = top_idx, top = top, second_idx = second_idx,
       second = second, ratio = ratio)
}

getDominanceMask <- function(spec_mat, ratio_thresh = DEFAULT_RATIO_THRESH) {
  spec_mat <- .specMatrix(spec_mat)
  if (length(ratio_thresh) != 1L || !is.finite(ratio_thresh) || ratio_thresh < 1) {
    stop("ratio_thresh must be one finite number >= 1.")
  }
  t2 <- .rowTop2(spec_mat)
  out <- matrix(FALSE, nrow(spec_mat), ncol(spec_mat), dimnames = dimnames(spec_mat))
  hit <- which(!is.na(t2$ratio) & .aboveThreshold(t2$ratio, ratio_thresh))
  out[cbind(hit, t2$top_idx[hit])] <- TRUE
  out
}

.specNamesFromMask <- function(mask) {
  vapply(seq_len(nrow(mask)), function(i)
    paste(colnames(mask)[mask[i, ]], collapse = ", "), character(1))
}

specPeakCore <- function(mat) {
  mat <- .specMatrix(mat)
  tiers <- getSpecTiers(mat)
  t2 <- .rowTop2(mat)
  high <- rowSums(tiers$high) > 0
  dominant <- rowSums(getDominanceMask(mat)) > 0
  data.frame(
    peak = if (nrow(mat)) rownames(mat) else character(0),
    top = colnames(mat)[t2$top_idx],
    top_score = t2$top,
    top_share = t2$top^2,
    second = colnames(mat)[t2$second_idx],
    second_score = t2$second,
    ratio = t2$ratio,
    is_dominant_r20 = dominant,
    is_high = high,
    class = ifelse(high, "specific", "shared"),
    n_high = as.integer(rowSums(tiers$high)),
    n_moderate = as.integer(rowSums(tiers$moderate)),
    n_low = as.integer(rowSums(tiers$low)),
    high_units = .specNamesFromMask(tiers$high),
    moderate_units = .specNamesFromMask(tiers$moderate),
    stringsAsFactors = FALSE
  )
}

.specAgreement <- function(high, dominant) {
  ifelse(high, ifelse(dominant, "both", "high_only"),
         ifelse(dominant, "dominant_only", "neither"))
}

.specWrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

addSpecPeakTable <- function(gwas_mat, snps, se_gwas, results_dir,
                             cell_lineage = NULL, gtf_path = NULL,
                             file = "spec_peak_summary.txt") {
  gwas_mat <- .specMatrix(gwas_mat)
  idx <- match(rownames(gwas_mat), rownames(se_gwas))
  if (anyNA(idx) || anyDuplicated(rownames(gwas_mat))) {
    stop("gwas_mat must have unique peak names present in se_gwas.")
  }
  se_gwas <- se_gwas[idx, , drop = FALSE]
  out <- specPeakCore(gwas_mat)
  names(out)[match(c("top", "second", "high_units", "moderate_units"), names(out))] <-
    c("top_cell", "second_cell", "high_cells", "moderate_cells")
  out$call_agreement <- .specAgreement(out$is_high, out$is_dominant_r20)

  if (!is.null(cell_lineage)) {
    lin <- specPeakCore(getLineageMatrix(gwas_mat, cell_lineage))
    names(lin)[match(c("top", "second", "high_units", "moderate_units"), names(lin))] <-
      c("top_lineage", "second_lineage", "high_lineages", "moderate_lineages")
    suffix <- setdiff(names(lin), c("peak", "top_lineage", "second_lineage",
                                   "high_lineages", "moderate_lineages"))
    names(lin)[match(suffix, names(lin))] <- paste0(suffix, "_lineage")
    out <- cbind(out, lin[, setdiff(names(lin), "peak"), drop = FALSE])
    out$call_agreement_lineage <- .specAgreement(out$is_high_lineage,
                                                out$is_dominant_r20_lineage)
    out$resolution_class <- ifelse(out$is_high, "cell_specific",
                                  ifelse(out$is_high_lineage, "lineage_specific", "shared"))
  } else {
    out$resolution_class <- ifelse(out$is_high, "cell_specific", "shared")
  }

  access_name <- if ("norm" %in% assayNames(se_gwas)) "norm" else "raw"
  if (!access_name %in% assayNames(se_gwas)) {
    stop("se_gwas must contain a norm or raw accessibility assay.")
  }
  access <- as.matrix(assay(se_gwas, access_name))
  out$access_max <- if (nrow(access)) apply(access, 1, max) else numeric(0)
  out$access_mean <- rowMeans(access)
  out <- cbind(out, .peakTssAnnotation(rowRanges(se_gwas), gtf_path))

  map <- getSpecSignalVariantPeakMap(snps, se_gwas)
  collapse <- function(values) {
    by_peak <- split(values, map$peak)
    vapply(out$peak, function(p) paste(unique(by_peak[[p]]), collapse = ", "),
           character(1), USE.NAMES = FALSE)
  }
  out$locus <- collapse(map$signal)
  out$gwas_snps <- collapse(map$variant)
  front <- c("locus", "peak", "gwas_snps", "resolution_class", "call_agreement")
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]
  .specWrite(out, file.path(results_dir, file))
  out
}

# Legacy table and plotting APIs share the same unrounded dominance rule.
addDominance <- function(mat, ratio_thresh = DEFAULT_RATIO_THRESH) {
  mat <- .specMatrix(mat)
  t2 <- .rowTop2(mat)
  data.frame(maximum_ct = colnames(mat)[t2$top_idx], ratio = t2$ratio,
             is_dominant = rowSums(getDominanceMask(mat, ratio_thresh)) > 0,
             peak = if (nrow(mat)) rownames(mat) else character(0), stringsAsFactors = FALSE)
}

addDominanceLineage <- function(mat, cell_lineage, ratio_thresh = DEFAULT_RATIO_THRESH) {
  out <- addDominance(getLineageMatrix(mat, cell_lineage), ratio_thresh)
  names(out)[names(out) == "maximum_ct"] <- "maximum_lineage"
  out
}

addSpecificitySummary <- function(mat_l2, cell_cor = NULL,
                                  ratio_thresh = DEFAULT_RATIO_THRESH,
                                  threshold_sq = 0.8) {
  mat_l2 <- .specMatrix(mat_l2)
  if (length(threshold_sq) != 1L || !is.finite(threshold_sq) ||
      threshold_sq <= 0 || threshold_sq > 1) stop("threshold_sq must be in (0, 1].")
  t2 <- .rowTop2(mat_l2)
  dom <- addDominance(mat_l2, ratio_thresh)
  if (is.null(cell_cor)) {
    cell_cor <- if (nrow(mat_l2) > 1L) cor(mat_l2, use = "pairwise.complete.obs") else
      matrix(NA_real_, ncol(mat_l2), ncol(mat_l2))
  }
  top <- lapply(seq_len(nrow(mat_l2)), function(i) {
    ord <- order(mat_l2[i, ], decreasing = TRUE)
    k <- which(cumsum(mat_l2[i, ord]^2) >= threshold_sq)[1]
    if (is.na(k)) k <- ncol(mat_l2)
    ord[seq_len(k)]
  })
  data.frame(
    maximum_cell = dom$maximum_ct, maximum_score = t2$top,
    dominant_ratio = dom$ratio, is_dominant = dom$is_dominant,
    k_cells = lengths(top),
    top_k_cells = vapply(top, function(i) paste(colnames(mat_l2)[i], collapse = ","), character(1)),
    coherence = vapply(top, function(i) {
      if (length(i) < 2L) return(NA_real_)
      sub <- cell_cor[i, i, drop = FALSE]
      mean(sub[upper.tri(sub)])
    }, numeric(1)),
    peak = if (nrow(mat_l2)) rownames(mat_l2) else character(0), stringsAsFactors = FALSE
  )
}

.wrapCellList <- function(x, width = 28) {
  vapply(x, function(label)
    paste(strwrap(gsub("\\s*[;:,]\\s*", ", ", label), width = width), collapse = "\n"),
    character(1), USE.NAMES = FALSE)
}
