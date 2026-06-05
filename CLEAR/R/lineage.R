# ============================================================
# lineage.R
# Cell-type lineage clustering, palette, and lineage-summed
# specificity matrix.
# ============================================================

# ------------------------------------------------------------
# addLineageMap
# Hierarchical clustering of cell types from the L2 matrix and
# k-group cut; returns the dendrogram, lineage labels, and a
# hierarchical colour palette (one master hue per lineage).
#
# k_lineage: integer number of lineage groups to cut the
# dendrogram into. User-defined; default 4. Must be in
# [2, ncol(mat_l2) - 1].
# ------------------------------------------------------------
addLineageMap <- function(mat_l2, k_lineage = 4L) {
  if (inherits(mat_l2, "Matrix")) mat_l2 <- as.matrix(mat_l2)
  storage.mode(mat_l2) <- "double"
  if (is.null(colnames(mat_l2))) {
    colnames(mat_l2) <- paste0("ct", seq_len(ncol(mat_l2)))
  }

  k_lineage <- as.integer(k_lineage)
  if (is.na(k_lineage) || k_lineage < 2L || k_lineage >= ncol(mat_l2)) {
    stop("k_lineage must be an integer in [2, ncol(mat_l2) - 1]; got ", k_lineage)
  }

  cell_cor <- cor(mat_l2, method = "pearson")
  hc <- hclust(as.dist(1 - cell_cor), method = "average")
  lineage_ids <- cutree(hc, k = k_lineage)
  lineage_labels <- tapply(
    names(lineage_ids), lineage_ids,
    function(x) paste(sort(x), collapse = ":")
  )
  cell_lineage <- setNames(
    as.character(lineage_labels[as.character(lineage_ids)]),
    names(lineage_ids)
  )

  palette <- .buildLineagePalette(lineage_ids)

  master_hues <- c("#D6354A", "#2E8B57", "#2270B5", "#E68A00", "#333333")
  if (length(lineage_labels) > length(master_hues)) {
    master_hues <- c(master_hues,
                     grDevices::hcl.colors(length(lineage_labels) - length(master_hues), "Set2"))
  }
  lineage_colours <- setNames(
    master_hues[seq_along(lineage_labels)],
    as.character(lineage_labels)
  )

  list(
    cell_cor        = cell_cor,
    hc              = hc,
    lineage_ids     = lineage_ids,
    lineage_labels  = lineage_labels,
    cell_lineage    = cell_lineage,
    k_lineage       = k_lineage,
    palette         = palette,
    lineage_colours = lineage_colours
  )
}

# Internal: hierarchical colour palette (master hue per lineage,
# shades per cell type within lineage).
.buildLineagePalette <- function(lineage_ids) {
  master_hues <- c(
    "#D6354A",  # red-pink
    "#2E8B57",  # green
    "#2270B5",  # blue
    "#E68A00",  # amber
    "#333333"   # charcoal
  )

  groups <- split(names(lineage_ids), lineage_ids)
  n_groups <- length(groups)
  if (n_groups > length(master_hues)) {
    master_hues <- c(master_hues,
                     grDevices::hcl.colors(n_groups - length(master_hues), "Set2"))
  }

  pal <- character()
  for (i in seq_along(groups)) {
    cts  <- sort(groups[[i]])
    n_ct <- length(cts)
    if (n_ct == 1L) {
      shades <- master_hues[i]
    } else {
      base_hsv <- grDevices::rgb2hsv(grDevices::col2rgb(master_hues[i]))
      h <- base_hsv[1]; s <- base_hsv[2]
      vals <- seq(0.55, 0.95, length.out = n_ct)
      if (s < 0.1) {
        shades <- vapply(seq_len(n_ct), function(j) grDevices::hsv(0, 0, vals[j]), character(1))
      } else {
        sats <- seq(1.0, 0.35, length.out = n_ct)
        shades <- vapply(seq_len(n_ct), function(j) grDevices::hsv(h, sats[j], vals[j]), character(1))
      }
    }
    names(shades) <- cts
    pal <- c(pal, shades)
  }
  pal
}

# ------------------------------------------------------------
# getLineageMatrix
# Lineage-aggregated specificity matrix. For each peak, sum L2^2
# within each lineage group and take the sqrt; rows preserve
# unit norm because L2 rows of mat_l2 do.
# ------------------------------------------------------------
getLineageMatrix <- function(mat_l2, cell_lineage) {
  lineage_vec <- cell_lineage[colnames(mat_l2)]
  lineage_levels <- unique(as.character(lineage_vec))
  mat_sq <- mat_l2^2
  lineage_mat <- sapply(lineage_levels, function(lg) {
    sqrt(rowSums(mat_sq[, lineage_vec == lg, drop = FALSE]))
  })
  lineage_mat <- as.matrix(lineage_mat)
  rownames(lineage_mat) <- rownames(mat_l2)
  lineage_mat
}