# ============================================================
# dominance.R
# Peak-level dominance, consensus, and specificity summaries
# ============================================================

DEFAULT_RATIO_THRESH <- 2

compute_dominance <- function(mat, ratio_thresh = 2) {
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  if (n_rows == 0L) {
    return(data.frame(maximum_ct = character(), ratio = numeric(), is_dominant = logical(), peak = character()))
  }
  if (is.null(colnames(mat))) colnames(mat) <- paste0("ct", seq_len(n_cols))
  
  if (n_cols < 2L) {
    return(data.frame(
      maximum_ct = rep(colnames(mat)[1], n_rows),
      ratio = rep(NA_real_, n_rows),
      is_dominant = rep(FALSE, n_rows),
      peak = rownames(mat)
    ))
  }
  
  cell_names <- colnames(mat)
  max_idx <- max.col(mat, ties.method = "first")
  max_val <- mat[cbind(seq_len(n_rows), max_idx)]
  
  mat2 <- mat
  mat2[cbind(seq_len(n_rows), max_idx)] <- -Inf
  second_val <- mat2[cbind(seq_len(n_rows), max.col(mat2, ties.method = "first"))]
  
  ratio <- ifelse(second_val == 0, Inf, max_val / second_val)
  
  data.frame(
    maximum_ct = cell_names[max_idx],
    ratio = ratio,
    is_dominant = ratio > ratio_thresh,
    peak = rownames(mat),
    stringsAsFactors = FALSE
  )
}

compute_dominance_lineage <- function(mat, cell_lineage, ratio_thresh = 2) {
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  if (n_rows == 0L) {
    return(data.frame(maximum_lineage = character(), ratio = numeric(), is_dominant = logical(), peak = character()))
  }
  if (is.null(colnames(mat))) colnames(mat) <- paste0("ct", seq_len(n_cols))

  lineage_vec <- cell_lineage[colnames(mat)]
  missing_lineage <- is.na(lineage_vec) | !nzchar(lineage_vec)
  if (any(missing_lineage)) {
    lineage_vec[missing_lineage] <- colnames(mat)[missing_lineage]
  }

  lineage_levels <- unique(as.character(lineage_vec))
  mat_sq <- mat^2
  lineage_mat <- sapply(lineage_levels, function(lg) {
    rowSums(mat_sq[, lineage_vec == lg, drop = FALSE])
  })
  lineage_mat <- as.matrix(lineage_mat)
  if (is.null(dim(lineage_mat))) {
    lineage_mat <- matrix(lineage_mat, ncol = 1)
    colnames(lineage_mat) <- lineage_levels[1]
  }
  rownames(lineage_mat) <- rownames(mat)

  if (ncol(lineage_mat) < 2L) {
    return(data.frame(
      maximum_lineage = rep(colnames(lineage_mat)[1], nrow(lineage_mat)),
      ratio = rep(NA_real_, nrow(lineage_mat)),
      is_dominant = rep(FALSE, nrow(lineage_mat)),
      peak = rownames(lineage_mat),
      stringsAsFactors = FALSE
    ))
  }

  max_idx <- max.col(lineage_mat, ties.method = "first")
  max_val <- lineage_mat[cbind(seq_len(nrow(lineage_mat)), max_idx)]

  lineage_mat2 <- lineage_mat
  lineage_mat2[cbind(seq_len(nrow(lineage_mat2)), max_idx)] <- -Inf
  second_val <- lineage_mat2[cbind(seq_len(nrow(lineage_mat2)), max.col(lineage_mat2, ties.method = "first"))]

  ratio <- ifelse(second_val == 0, Inf, max_val / second_val)

  data.frame(
    maximum_lineage = colnames(lineage_mat)[max_idx],
    ratio = ratio,
    is_dominant = ratio > ratio_thresh,
    peak = rownames(lineage_mat),
    stringsAsFactors = FALSE
  )
}

dominance_consensus <- function(specificity_matrix, peak_weights, group_labels) {
  cell_types <- colnames(specificity_matrix)
  groups     <- group_labels[cell_types]

  group_sums <- lapply(unique(groups), function(g) {
    cols <- cell_types[groups == g]
    rowSums(specificity_matrix[, cols, drop = FALSE])
  }) %>%
    setNames(unique(groups)) %>% as.data.frame()

  top_group_per_peak <- apply(group_sums, 1, function(x) names(which.max(x)))
  peak_spec_weight <- apply(group_sums, 1, max)
  peak_weighted_spec <- peak_spec_weight * peak_weights

  group_totals  <- tapply(peak_weighted_spec, top_group_per_peak, sum)
  top_group     <- names(which.max(group_totals))
  consensus     <- group_totals[top_group] / sum(group_totals)

  list(
    consensus    = unname(consensus),
    top_group    = top_group,
    group_totals = sort(group_totals / sum(group_totals), decreasing = TRUE)
  )
}

compute_specificity_summary <- function(mat_l2, cell_cor = NULL, ratio_thresh = 2, threshold_sq = 0.8) {
  if (is.null(cell_cor)) cell_cor <- cor(mat_l2, use = "pairwise.complete.obs")
  results <- apply(mat_l2, 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    x_sorted <- x[ord]
    cells_sorted <- colnames(mat_l2)[ord]
    cum_sq <- cumsum(x_sorted^2)
    k <- which(cum_sq >= threshold_sq)[1]
    if (is.na(k)) k <- length(x_sorted)
    top_cells <- cells_sorted[seq_len(k)]
    
    maximum_cell <- cells_sorted[1]
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
      maximum_cell = maximum_cell,
      maximum_score = as.numeric(maximum_score),
      dominant_ratio = as.numeric(dominant_ratio),
      is_dominant = as.character(dominant_ratio > ratio_thresh),
      k_cells = as.numeric(k),
      top_k_cells = paste(top_cells, collapse = ","),
      coherence = as.numeric(coherence)
    )
  })
  
  out <- as.data.frame(t(results), stringsAsFactors = FALSE)
  out$peak <- rownames(mat_l2)
  out$maximum_score <- as.numeric(out$maximum_score)
  out$dominant_ratio <- as.numeric(out$dominant_ratio)
  out$is_dominant <- as.logical(out$is_dominant)
  out$k_cells <- as.integer(out$k_cells)
  out$coherence <- as.numeric(out$coherence)
  rownames(out) <- NULL
  out
}
