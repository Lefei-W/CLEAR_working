.cellCounts <- function(se, n_cells = NULL) {
  cn <- colnames(se)
  if (is.character(n_cells) && length(n_cells) == 1L) {
    if (!file.exists(n_cells)) stop("Cell-count file does not exist: ", n_cells)
    tab <- utils::read.delim(n_cells, stringsAsFactors = FALSE)
    if (ncol(tab) < 2L) stop("Cell-count TSV must contain cell type and n_cells columns.")
    n_cells <- setNames(tab[[2]], tab[[1]])
  }
  if (is.null(n_cells)) {
    cd <- as.data.frame(colData(se))
    known <- c("ncells", "ncell", "numcells", "cellcounts", "cellcount",
               "ncellsingroup", "ncellspergroup", "cells")
    pos <- which(gsub("[^a-z0-9]", "", tolower(names(cd))) %in% known)
    if (!length(pos)) stop("Cell counts are required; supply n_cells or colData(se)$nCells.")
    n_cells <- setNames(cd[[pos[1]]], cn)
  }
  if (!is.numeric(n_cells) || any(!is.finite(n_cells)) || any(n_cells <= 0)) {
    stop("Cell counts must be finite, positive numbers.")
  }
  if (is.null(names(n_cells))) {
    if (length(n_cells) != length(cn)) stop("Supply one cell count per cell type.")
    names(n_cells) <- cn
  }
  if (anyDuplicated(names(n_cells)) || anyNA(match(cn, names(n_cells)))) {
    stop("Cell counts must uniquely cover every cell-type column.")
  }
  n_cells[cn]
}

.normaliseAccessibility <- function(mat, se, per_cell, depth, n_cells) {
  use_per_cell <- per_cell == "yes" ||
    (per_cell == "auto" && all(abs(mat - round(mat)) < 1e-8))
  cells <- if (use_per_cell) .cellCounts(se, n_cells) else NULL
  norm <- if (use_per_cell) sweep(mat, 2, cells, "/") else mat
  target <- NA_real_
  if (depth == "colsum") {
    totals <- colSums(norm)
    if (!any(totals > 0)) stop("Cannot depth-normalise an all-zero matrix.")
    if (any(totals == 0)) warning("Zero-total columns are left unscaled: ",
                                 paste(colnames(mat)[totals == 0], collapse = ", "))
    target <- stats::median(totals[totals > 0])
    norm <- sweep(norm, 2, ifelse(totals > 0, target / totals, 1), "*")
  } else if (depth == "quantile") {
    sorted <- vapply(seq_len(ncol(norm)), function(j) sort(norm[, j]),
                     numeric(nrow(norm)))
    reference <- rowMeans(matrix(sorted, nrow = nrow(norm), ncol = ncol(norm)))
    for (j in seq_len(ncol(norm))) {
      ord <- order(norm[, j])
      # Average the reference over tied blocks, not only their midpoint.
      norm[ord, j] <- stats::ave(reference, norm[ord, j], FUN = mean)
    }
  }
  list(norm = norm, metadata = list(
    per_cell = use_per_cell, per_cell_mode = per_cell, depth = depth,
    depth_target = target, n_cells = cells))
}

flagAccessiblePeaks <- function(se, min_max = NULL, min_max_quantile = NULL) {
  if (!is.null(min_max) && !is.null(min_max_quantile)) {
    stop("Choose either min_max or min_max_quantile, not both.")
  }
  check <- function(x, name) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 0) {
      stop(name, " must be one finite, non-negative number.")
    }
  }
  if (!is.null(min_max)) check(min_max, "min_max")
  if (!is.null(min_max_quantile)) {
    check(min_max_quantile, "min_max_quantile")
    if (min_max_quantile >= 1) stop("min_max_quantile must be in [0, 1).")
  }
  access <- as.matrix(assay(se, if ("norm" %in% assayNames(se)) "norm" else "raw"))
  if (!nrow(access)) stop("Cannot filter an empty peak universe.")
  mx <- apply(access, 1, max)
  mode <- "none"
  cut <- NA_real_
  pass <- rep(TRUE, length(mx))
  if (!is.null(min_max) || !is.null(min_max_quantile)) {
    mode <- if (!is.null(min_max)) "absolute" else "quantile"
    cut <- if (!is.null(min_max)) min_max else
      unname(stats::quantile(mx, min_max_quantile))
    pass <- if (cut == 0) mx > 0 else mx >= cut
  }
  rowData(se)$max_access <- mx
  rowData(se)$pass_access <- pass
  md <- metadata(se)
  md$CLEAR$access_filter <- list(mode = mode, cut = cut, quantile = min_max_quantile,
                                 n_peaks = length(mx), n_pass = sum(pass),
                                 n_drop = sum(!pass))
  metadata(se) <- md
  message("Accessibility filter (", mode, "): retained ", sum(pass), "/", length(pass),
          " peaks", if (identical(cut, 0)) "; zero-signal peaks excluded" else "")
  se
}

filterAccessiblePeaks <- function(se) {
  pass <- rowData(se)$pass_access
  if (is.null(pass)) stop("Call flagAccessiblePeaks() before filterAccessiblePeaks().")
  if (anyNA(pass) || !is.logical(pass)) stop("pass_access must be a logical vector without NA.")
  if (!any(pass)) stop("Accessibility filtering removed every peak; choose a lower cutoff.")
  se <- se[pass, , drop = FALSE]
  assays(se)[["raw_l2_rank"]] <- .rankDesc(assay(se, "raw_l2"))
  md <- metadata(se)
  md$CLEAR$n_peaks <- nrow(se)
  metadata(se) <- md
  se
}
