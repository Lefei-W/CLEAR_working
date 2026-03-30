# ============================================================
# lineage.R
# Cell-type lineage clustering and hierarchical colour palette
# ============================================================

build_lineage_map <- function(mat_l2, k_lineage = 4L) {
  if (inherits(mat_l2, "Matrix")) mat_l2 <- as.matrix(mat_l2)
  storage.mode(mat_l2) <- "double"
  if (is.null(colnames(mat_l2))) {
    colnames(mat_l2) <- paste0("ct", seq_len(ncol(mat_l2)))
  }

  cell_cor <- cor(mat_l2, method = "pearson")
  hc <- hclust(as.dist(1 - cell_cor), method = "average")
  lineage_ids <- cutree(hc, k = k_lineage)
  lineage_labels <- tapply(
    names(lineage_ids),
    lineage_ids,
    function(x) paste(sort(x), collapse = ":")
  )
  cell_lineage <- setNames(as.character(lineage_labels[as.character(lineage_ids)]), names(lineage_ids))

  palette <- build_lineage_palette(lineage_ids)

  list(
    cell_cor = cell_cor,
    hc = hc,
    lineage_ids = lineage_ids,
    lineage_labels = lineage_labels,
    cell_lineage = cell_lineage,
    k_lineage = k_lineage,
    palette = palette
  )
}

# Hierarchical colour palette: one master hue per lineage, shades per cell type
# Supports up to 8 lineages x 8 cell types each
build_lineage_palette <- function(lineage_ids) {
  master_hues <- c(
    "#2166AC",  # blue
    "#B2182B",  # red
    "#1B7837",  # green
    "#E08214",  # orange
    "#6A3D9A",  # purple
    "#A6761D",  # brown
    "#E7298A",  # pink
    "#66C2A5"   # teal
  )

  groups <- split(names(lineage_ids), lineage_ids)
  n_groups <- length(groups)
  if (n_groups > length(master_hues)) {
    master_hues <- c(master_hues, grDevices::hcl.colors(n_groups - length(master_hues), "Set2"))
  }

  pal <- character()
  for (i in seq_along(groups)) {
    cts <- sort(groups[[i]])
    n_ct <- length(cts)
    if (n_ct == 1L) {
      shades <- master_hues[i]
    } else {
      base_rgb <- grDevices::col2rgb(master_hues[i]) / 255
      alphas <- seq(1.0, 0.35, length.out = n_ct)
      shades <- vapply(alphas, function(a) {
        blended <- base_rgb * a + 1 * (1 - a)
        grDevices::rgb(blended[1], blended[2], blended[3])
      }, character(1))
    }
    names(shades) <- cts
    pal <- c(pal, shades)
  }

  pal
}
