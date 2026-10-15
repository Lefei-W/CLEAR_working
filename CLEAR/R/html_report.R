# ============================================================
# html_report.R
# A dependency-light, static HTML explorer for one CLEAR run.
#
# The report is written directly from R and contains its CSS, SVG and
# JavaScript inline. It therefore needs no Shiny server, Pandoc, DT or
# Plotly, and can be opened from a shared filesystem with file://.
#
# Main entry point:
#   writeSpecHtmlReport(...)
#
# The data preparation is deliberately separate and testable:
#   prepareSpecHtmlReportData(...)
# ============================================================

.specReportEscape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

.specReportDisplay <- function(x, empty = "\u2014") {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- empty
  .specReportEscape(x)
}

.specReportJsEncodedArray <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  encoded <- vapply(x, utils::URLencode, character(1), reserved = TRUE)
  paste0("[", paste0('"', encoded, '"', collapse = ","), "].map(decodeURIComponent)")
}

.specReportPalette <- function(ids, supplied = NULL) {
  ids <- as.character(ids)
  if (!length(ids)) return(setNames(character(0), character(0)))
  out <- rep(NA_character_, length(ids))
  names(out) <- ids
  if (!is.null(supplied) && length(supplied)) {
    hit <- match(ids, names(supplied))
    out[!is.na(hit)] <- as.character(supplied[hit[!is.na(hit)]])
  }
  missing <- is.na(out) | !nzchar(out)
  if (any(missing)) {
    out[missing] <- grDevices::hcl.colors(sum(missing), "Dark 3")
  }
  out
}

.specReportSnpFields <- function(snps) {
  n <- length(snps)
  if (!n) {
    return(list(signal = character(0), variant = character(0),
                region = character(0), has_signal = FALSE))
  }
  mc <- as.data.frame(S4Vectors::mcols(snps))
  region <- paste0(as.character(GenomicRanges::seqnames(snps)), ":",
                   BiocGenerics::start(snps), "-", BiocGenerics::end(snps))

  variant <- if ("names" %in% names(mc)) as.character(mc$names) else names(snps)
  if (is.null(variant)) variant <- rep("", n)
  bad_variant <- is.na(variant) | !nzchar(trimws(variant))
  variant[bad_variant] <- region[bad_variant]

  has_signal <- "signal" %in% names(mc)
  signal <- if (has_signal) as.character(mc$signal) else rep("unassigned", n)
  signal[is.na(signal) | !nzchar(trimws(signal))] <- "unassigned"
  signal <- trimws(signal)

  list(signal = signal, variant = variant, region = region,
       has_signal = has_signal)
}

# One row per exact signal x variant x retained-peak overlap. This is the
# canonical report join and avoids filtering comma-collapsed signal strings.
getSpecSignalVariantPeakMap <- function(snps, se_gwas) {
  empty <- data.frame(signal = character(0), variant = character(0),
                      variant_region = character(0), peak = character(0),
                      stringsAsFactors = FALSE)
  if (is.null(snps) || !length(snps) || is.null(se_gwas) || !nrow(se_gwas)) {
    return(empty)
  }
  hits <- GenomicRanges::findOverlaps(snps,
                                      SummarizedExperiment::rowRanges(se_gwas))
  if (!length(hits)) return(empty)
  fields <- .specReportSnpFields(snps)
  peak_names <- rownames(se_gwas)
  if (is.null(peak_names)) {
    peak_names <- paste0(as.character(GenomicRanges::seqnames(
      SummarizedExperiment::rowRanges(se_gwas))), ":",
      BiocGenerics::start(SummarizedExperiment::rowRanges(se_gwas)), "-",
      BiocGenerics::end(SummarizedExperiment::rowRanges(se_gwas)))
  }
  out <- data.frame(
    signal         = fields$signal[S4Vectors::queryHits(hits)],
    variant        = fields$variant[S4Vectors::queryHits(hits)],
    variant_region = fields$region[S4Vectors::queryHits(hits)],
    peak            = peak_names[S4Vectors::subjectHits(hits)],
    stringsAsFactors = FALSE
  )
  out <- unique(out)
  out[order(out$signal, out$peak, out$variant), , drop = FALSE]
}

.specReportMatrix <- function(x, row_ids = NULL) {
  if (inherits(x, "Matrix")) x <- as.matrix(x)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!is.null(row_ids)) {
    if (!length(row_ids) && !nrow(x)) {
      rownames(x) <- character(0)
      return(x)
    }
    if (is.null(rownames(x))) stop("Report matrix is missing peak row names.")
    hit <- match(row_ids, rownames(x))
    if (anyNA(hit)) {
      stop("Report matrix is missing ", sum(is.na(hit)),
           " retained variant peak(s).")
    }
    x <- x[hit, , drop = FALSE]
  }
  x
}

.specReportCorrelation <- function(mat) {
  if (!ncol(mat)) return(matrix(numeric(0), 0, 0))
  if (ncol(mat) == 1L) {
    out <- matrix(1, 1, 1, dimnames = list(colnames(mat), colnames(mat)))
    return(out)
  }
  stats::cor(mat, use = "pairwise.complete.obs")
}

.specReportCategory <- function(peak_tbl, peak_ids, high_col, label_col) {
  if (!length(peak_ids) || !nrow(peak_tbl) ||
      !all(c("peak", high_col, label_col) %in% names(peak_tbl))) {
    return(list(n_high = 0L, category = "no_high", sole = NA_character_))
  }
  d <- peak_tbl[match(peak_ids, as.character(peak_tbl$peak)), , drop = FALSE]
  high <- suppressWarnings(as.logical(d[[high_col]]))
  high[is.na(high)] <- FALSE
  labels <- as.character(d[[label_col]])
  labels <- unique(labels[high & !is.na(labels) & nzchar(labels)])
  n_high <- sum(high)
  category <- if (n_high == 0L) "no_high" else if (n_high == 1L) "single" else
    if (length(labels) == 1L) "coherent" else "complex"
  list(n_high = n_high, category = category,
       sole = if (length(labels) == 1L) labels else NA_character_)
}

.specReportTopShare <- function(share_mat, peak_ids) {
  if (!length(peak_ids) || !nrow(share_mat) || !ncol(share_mat)) {
    return(c(unit = NA_character_, share = NA_character_))
  }
  idx <- match(peak_ids, rownames(share_mat))
  idx <- idx[!is.na(idx)]
  if (!length(idx)) return(c(unit = NA_character_, share = NA_character_))
  means <- colMeans(share_mat[idx, , drop = FALSE], na.rm = TRUE)
  means[!is.finite(means)] <- NA_real_
  if (all(is.na(means))) return(c(unit = NA_character_, share = NA_character_))
  j <- which.max(replace(means, is.na(means), -Inf))
  c(unit = names(means)[j], share = as.character(means[j]))
}

.specReportCollapse <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  paste(x, collapse = ", ")
}

# Build the small, report-specific data model. Only the correlation matrices
# and variant-peak contribution shares are retained; the genome-wide matrix is
# not serialized into the HTML.
prepareSpecHtmlReportData <- function(
    se_name, trait_name, mat_spec, gwas_mat, se_gwas, snps,
    peak_cell, locus_categorisation = data.frame(),
    global_enrichment_cell = data.frame(), global_enrichment_lineage = data.frame(),
    lineage_obj, input_summary = data.frame()) {

  mat_spec <- .specReportMatrix(mat_spec)
  if (is.null(colnames(mat_spec)) || !ncol(mat_spec)) {
    stop("mat_spec must have at least one named cell-type column.")
  }

  peak_ids <- if (is.null(se_gwas)) character(0) else rownames(se_gwas)
  if (is.null(peak_ids) && !is.null(se_gwas) && nrow(se_gwas)) {
    stop("se_gwas is missing peak row names.")
  }
  gwas_mat <- .specReportMatrix(gwas_mat, peak_ids)
  if (!identical(colnames(gwas_mat), colnames(mat_spec))) {
    hit <- match(colnames(mat_spec), colnames(gwas_mat))
    if (anyNA(hit)) stop("gwas_mat and mat_spec have different cell-type columns.")
    gwas_mat <- gwas_mat[, hit, drop = FALSE]
  }

  cell_lineage <- lineage_obj$cell_lineage
  if (is.null(cell_lineage) || is.null(names(cell_lineage))) {
    stop("lineage_obj$cell_lineage must be a named cell-type to lineage map.")
  }
  if (anyNA(match(colnames(mat_spec), names(cell_lineage)))) {
    stop("The lineage map does not cover every mat_spec cell type.")
  }
  cell_lineage <- cell_lineage[colnames(mat_spec)]

  cell_share <- gwas_mat^2
  if (nrow(gwas_mat)) {
    lineage_gwas <- getLineageMatrix(gwas_mat, cell_lineage)
    lineage_share <- lineage_gwas^2
  } else {
    lineage_ids <- unique(as.character(cell_lineage))
    lineage_share <- matrix(numeric(0), nrow = 0L, ncol = length(lineage_ids),
                            dimnames = list(character(0), lineage_ids))
  }
  lineage_global <- getLineageMatrix(mat_spec, cell_lineage)

  cell_ids <- colnames(mat_spec)
  lineage_ids <- colnames(lineage_global)
  cell_colours <- .specReportPalette(cell_ids, lineage_obj$palette)
  lineage_colours <- .specReportPalette(lineage_ids,
                                        lineage_obj$lineage_colours)

  map <- getSpecSignalVariantPeakMap(snps, se_gwas)
  snp_fields <- .specReportSnpFields(snps)
  signal_ids <- unique(c(snp_fields$signal, map$signal))
  signal_ids <- signal_ids[!is.na(signal_ids) & nzchar(signal_ids)]

  if (!is.data.frame(peak_cell)) peak_cell <- as.data.frame(peak_cell)
  if (!nrow(peak_cell) && length(peak_ids)) {
    peak_cell <- data.frame(peak = peak_ids, stringsAsFactors = FALSE)
  }
  if (nrow(peak_cell) && !"peak" %in% names(peak_cell)) {
    stop("peak_cell must contain a 'peak' column.")
  }

  # Keep one category row per locus at the canonical s50 threshold.
  cat_tbl <- as.data.frame(locus_categorisation)
  if (nrow(cat_tbl) && "threshold" %in% names(cat_tbl)) {
    keep <- as.character(cat_tbl$threshold) == "s50"
    if (any(keep)) cat_tbl <- cat_tbl[keep, , drop = FALSE]
  }
  if (nrow(cat_tbl) && "locus" %in% names(cat_tbl)) {
    cat_tbl <- cat_tbl[!duplicated(as.character(cat_tbl$locus)), , drop = FALSE]
  }

  map_peaks <- split(map$peak, map$signal)
  map_variants <- split(map$variant, map$signal)
  all_variants <- split(snp_fields$variant, snp_fields$signal)
  gene_by_peak <- if (nrow(peak_cell) && "nearest_gene" %in% names(peak_cell))
    setNames(as.character(peak_cell$nearest_gene), as.character(peak_cell$peak)) else
      setNames(character(0), character(0))

  n_signal <- length(signal_ids)
  summary_rows <- vector("list", n_signal)
  for (i in seq_along(signal_ids)) {
    sig <- signal_ids[i]
    peaks <- unique(as.character(map_peaks[[sig]]))
    peaks <- peaks[!is.na(peaks) & nzchar(peaks)]
    variants_in_peaks <- unique(as.character(map_variants[[sig]]))
    variants_in_peaks <- variants_in_peaks[!is.na(variants_in_peaks) &
                                             nzchar(variants_in_peaks)]
    variants_all <- unique(as.character(all_variants[[sig]]))
    variants_all <- variants_all[!is.na(variants_all) & nzchar(variants_all)]

    cell_fallback <- .specReportCategory(peak_cell, peaks, "is_high", "top_cell")
    lin_fallback <- .specReportCategory(peak_cell, peaks,
                                        "is_high_lineage", "top_lineage")
    cell_top <- .specReportTopShare(cell_share, peaks)
    lin_top <- .specReportTopShare(lineage_share, peaks)

    cat_idx <- if (nrow(cat_tbl) && "locus" %in% names(cat_tbl))
      match(sig, as.character(cat_tbl$locus)) else NA_integer_
    cat_get <- function(column, fallback = NA) {
      if (is.na(cat_idx) || !column %in% names(cat_tbl)) return(fallback)
      value <- cat_tbl[[column]][cat_idx]
      if (length(value) == 0L || is.na(value)) fallback else value
    }
    retained <- length(peaks) > 0L
    genes <- .specReportCollapse(unname(gene_by_peak[peaks]))

    summary_rows[[i]] <- data.frame(
      signal = sig,
      coverage = if (retained) "retained_peaks" else "no_retained_peak",
      n_variants = length(variants_all),
      n_variants_in_peaks = length(variants_in_peaks),
      variant_ids = .specReportCollapse(variants_in_peaks),
      n_peaks = length(peaks),
      n_high_cell = if (retained) as.integer(cat_get(
        "variant_cell_n_high", cell_fallback$n_high)) else 0L,
      cell_category = if (retained) as.character(cat_get(
        "variant_cell_category", cell_fallback$category)) else "not_evaluated",
      sole_cell = if (retained) as.character(cat_get(
        "sole_contributor_cell", cell_fallback$sole)) else NA_character_,
      top_cell = if (retained) cell_top[["unit"]] else NA_character_,
      top_cell_share = if (retained) as.numeric(cell_top[["share"]]) else NA_real_,
      cell_has_enriched = if (retained) as.logical(cat_get(
        "variant_cell_has_enriched", NA)) else NA,
      n_high_lineage = if (retained) as.integer(cat_get(
        "variant_lin_n_high", lin_fallback$n_high)) else 0L,
      lineage_category = if (retained) as.character(cat_get(
        "variant_lin_category", lin_fallback$category)) else "not_evaluated",
      sole_lineage = if (retained) as.character(cat_get(
        "sole_contributor_lineage", lin_fallback$sole)) else NA_character_,
      top_lineage = if (retained) lin_top[["unit"]] else NA_character_,
      top_lineage_share = if (retained) as.numeric(lin_top[["share"]]) else NA_real_,
      lineage_has_enriched = if (retained) as.logical(cat_get(
        "variant_lin_has_enriched", NA)) else NA,
      nearest_genes = genes,
      stringsAsFactors = FALSE
    )
  }
  signal_summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else
    data.frame(signal = character(0), coverage = character(0),
               n_variants = integer(0), n_variants_in_peaks = integer(0),
               variant_ids = character(0), n_peaks = integer(0),
               n_high_cell = integer(0),
               cell_category = character(0), sole_cell = character(0),
               top_cell = character(0), top_cell_share = numeric(0),
               cell_has_enriched = logical(0), n_high_lineage = integer(0),
               lineage_category = character(0), sole_lineage = character(0),
               top_lineage = character(0), top_lineage_share = numeric(0),
               lineage_has_enriched = logical(0), nearest_genes = character(0),
               stringsAsFactors = FALSE)

  if (nrow(signal_summary)) {
    ord <- order(signal_summary$coverage != "retained_peaks",
                 -signal_summary$n_high_cell, -signal_summary$n_peaks,
                 signal_summary$signal)
    signal_summary <- signal_summary[ord, , drop = FALSE]
    rownames(signal_summary) <- NULL
  }
  signal_summary$signal_index <- seq_len(nrow(signal_summary))

  # Reorder the peak summary to the exact gwas matrix rows, then replace its
  # comma-collapsed identifiers with values derived from the canonical map.
  peak_summary <- data.frame(peak = peak_ids, stringsAsFactors = FALSE)
  if (length(peak_ids) && nrow(peak_cell)) {
    idx <- match(peak_ids, as.character(peak_cell$peak))
    for (nm in setdiff(names(peak_cell), "peak")) {
      peak_summary[[nm]] <- peak_cell[[nm]][idx]
    }
  }
  signals_by_peak <- split(map$signal, map$peak)
  variants_by_peak <- split(map$variant, map$peak)
  peak_signals <- lapply(peak_ids, function(p) unique(as.character(signals_by_peak[[p]])))
  peak_variants <- lapply(peak_ids, function(p) unique(as.character(variants_by_peak[[p]])))
  if (length(peak_ids)) {
    peak_summary$locus <- vapply(peak_signals, .specReportCollapse, character(1))
    peak_summary$gwas_snps <- vapply(peak_variants, .specReportCollapse, character(1))
  } else {
    peak_summary$locus <- character(0)
    peak_summary$gwas_snps <- character(0)
  }
  signal_index <- setNames(signal_summary$signal_index, signal_summary$signal)
  peak_signal_indices <- lapply(peak_signals, function(x) {
    out <- unname(signal_index[x])
    as.integer(out[!is.na(out)])
  })

  list(
    se_name = se_name,
    trait_name = trait_name,
    input_summary = as.data.frame(input_summary),
    signal_metadata_present = snp_fields$has_signal,
    map = map,
    signal_summary = signal_summary,
    peak_summary = peak_summary,
    peak_signals = peak_signals,
    peak_signal_indices = peak_signal_indices,
    cell_share = cell_share,
    lineage_share = lineage_share,
    cell_correlation = .specReportCorrelation(mat_spec),
    lineage_correlation = .specReportCorrelation(lineage_global),
    global_enrichment_cell = as.data.frame(global_enrichment_cell),
    global_enrichment_lineage = as.data.frame(global_enrichment_lineage),
    cell_ids = cell_ids,
    lineage_ids = lineage_ids,
    cell_colours = unname(cell_colours[cell_ids]),
    lineage_colours = unname(lineage_colours[lineage_ids]),
    cell_lineage = as.character(cell_lineage),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}

.specReportShort <- function(x, n = 48L) {
  x <- as.character(x)
  ifelse(nchar(x) <= n, x, paste0(substr(x, 1L, n - 1L), "\u2026"))
}

.specReportNum <- function(x, digits = 3L) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "fg", flag = "#"))
}

.specReportCorOrder <- function(cor_mat) {
  n <- nrow(cor_mat)
  if (n <= 1L) return(seq_len(n))
  clean <- cor_mat
  clean[!is.finite(clean)] <- 0
  diag(clean) <- 1
  tryCatch(stats::hclust(stats::dist(clean), method = "ward.D2")$order,
           error = function(e) seq_len(n))
}

.specReportCombinedSvg <- function(cor_mat, enrichment, id_col, colours,
                                   level_label) {
  n <- nrow(cor_mat)
  if (!n) return('<div class="empty-state">Correlation unavailable.</div>')
  ord <- .specReportCorOrder(cor_mat)
  cor_mat <- cor_mat[ord, ord, drop = FALSE]
  ids <- rownames(cor_mat)
  colours <- .specReportPalette(ids, setNames(colours, names(colours)))[ids]

  if (nrow(enrichment) && id_col %in% names(enrichment)) {
    ge <- enrichment[match(ids, as.character(enrichment[[id_col]])), , drop = FALSE]
    z <- if ("z" %in% names(ge)) as.numeric(ge$z) else rep(NA_real_, n)
    fdr <- if ("fdr" %in% names(ge)) as.numeric(ge$fdr) else rep(NA_real_, n)
    stars <- if ("fdr_stars" %in% names(ge)) as.character(ge$fdr_stars) else
      rep("", n)
  } else {
    z <- fdr <- rep(NA_real_, n)
    stars <- rep("", n)
  }

  cell <- max(16, min(25, floor(540 / max(n, 1))))
  heat_x <- 14
  heat_y <- 118
  heat_w <- n * cell
  label_x <- heat_x + heat_w + 18
  label_chars <- min(52, max(18, max(nchar(ids), na.rm = TRUE)))
  label_w <- min(360, max(145, label_chars * 6.4))
  # Keep the two panel titles from colliding on small cell-type panels.
  bar_x <- max(label_x + label_w, heat_x + 390)
  bar_w <- 300
  info_x <- bar_x + bar_w + 12
  # label_w is fractional (label_chars * 6.4) and propagates into info_x, so
  # round up before the viewBox below: sprintf's %d errors on a non-integral
  # double rather than truncating it.
  width <- ceiling(info_x + 155)
  height <- ceiling(heat_y + n * cell + 54)

  finite_z <- z[is.finite(z)]
  zlim <- if (length(finite_z)) max(1, max(abs(finite_z))) else 1
  zero_x <- bar_x + bar_w / 2
  half_scale <- (bar_w / 2 - 8) / zlim
  pal <- grDevices::colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(201)
  gradient_id <- paste0("cor-gradient-",
                        gsub("[^a-z0-9]+", "-", tolower(level_label)))
  cor_colour <- function(v) {
    if (!is.finite(v)) return("#E5E7EB")
    pal[pmax(1L, pmin(201L, round((v + 1) * 100) + 1L))]
  }

  out <- c(
    sprintf('<svg class="combined-svg" viewBox="0 0 %d %d" role="img" aria-label="%s specificity correlation heatmap with global enrichment">',
            width, height, .specReportEscape(level_label)),
    sprintf('<defs><linearGradient id="%s" x1="0%%" x2="100%%"><stop offset="0%%" stop-color="#2166AC"/><stop offset="50%%" stop-color="#FFFFFF"/><stop offset="100%%" stop-color="#B2182B"/></linearGradient></defs>',
            gradient_id),
    sprintf('<text x="%g" y="22" class="svg-title">Genome-wide specificity correlation</text>', heat_x),
    sprintf('<text x="%g" y="22" class="svg-title">Global rank enrichment</text>', bar_x),
    sprintf('<text x="%g" y="43" class="svg-note">Pearson r, clustered symmetrically</text>', heat_x),
    sprintf('<text x="%g" y="43" class="svg-note">Signed z; stars show BH FDR</text>', bar_x),
    sprintf('<rect x="%g" y="57" width="150" height="10" fill="url(#%s)" stroke="#CBD5E1"/>',
            heat_x, gradient_id),
    sprintf('<text x="%g" y="82" class="svg-note">-1</text>', heat_x),
    sprintf('<text x="%g" y="82" class="svg-note" text-anchor="middle">0</text>', heat_x + 75),
    sprintf('<text x="%g" y="82" class="svg-note" text-anchor="end">1</text>', heat_x + 150),
    sprintf('<line x1="%g" x2="%g" y1="%g" y2="%g" stroke="#475569" stroke-width="1"/>',
            zero_x, zero_x, heat_y - 3, heat_y + n * cell)
  )

  for (j in seq_len(n)) {
    x <- heat_x + (j - 0.5) * cell
    label <- .specReportEscape(.specReportShort(ids[j], 34L))
    out <- c(out, sprintf(
      '<text x="%g" y="%g" class="svg-axis" transform="rotate(-55 %g %g)" text-anchor="start"><title>%s</title>%s</text>',
      x, heat_y - 7, x, heat_y - 7, .specReportEscape(ids[j]), label))
  }

  for (i in seq_len(n)) {
    y <- heat_y + (i - 1) * cell
    for (j in seq_len(n)) {
      v <- cor_mat[i, j]
      tip <- if (is.finite(v)) sprintf("%s vs %s: r = %.4f", ids[i], ids[j], v) else
        sprintf("%s vs %s: correlation unavailable", ids[i], ids[j])
      out <- c(out, sprintf(
        '<rect x="%g" y="%g" width="%g" height="%g" fill="%s" stroke="#FFFFFF" stroke-width="0.35"><title>%s</title></rect>',
        heat_x + (j - 1) * cell, y, cell, cell, cor_colour(v),
        .specReportEscape(tip)))
    }

    cy <- y + cell * 0.69
    out <- c(out,
      sprintf('<rect x="%g" y="%g" width="8" height="8" rx="2" fill="%s"/>',
              label_x, y + (cell - 8) / 2, colours[i]),
      sprintf('<text x="%g" y="%g" class="svg-row"><title>%s</title>%s</text>',
              label_x + 13, cy, .specReportEscape(ids[i]),
              .specReportEscape(.specReportShort(ids[i], 48L))))

    if (is.finite(z[i])) {
      end_x <- zero_x + z[i] * half_scale
      rect_x <- min(zero_x, end_x)
      rect_w <- max(1, abs(end_x - zero_x))
      fill <- if (is.finite(fdr[i]) && fdr[i] < 0.05 && z[i] > 0) "#B2182B" else
        if (z[i] < 0) "#2166AC" else "#94A3B8"
      tip <- paste0(ids[i], ": z = ", formatC(z[i], digits = 4, format = "fg"),
                    if (is.finite(fdr[i])) paste0(", FDR = ",
                      formatC(fdr[i], digits = 4, format = "g")) else ", FDR unavailable")
      out <- c(out, sprintf(
        '<rect x="%g" y="%g" width="%g" height="%g" rx="2" fill="%s"><title>%s</title></rect>',
        rect_x, y + cell * 0.18, rect_w, cell * 0.64, fill,
        .specReportEscape(tip)),
        sprintf('<text x="%g" y="%g" class="svg-value">z=%s  %s</text>',
                info_x, cy, .specReportEscape(.specReportNum(z[i], 3L)),
                .specReportEscape(stars[i])))
    } else {
      out <- c(out, sprintf('<text x="%g" y="%g" class="svg-note">no enrichment result</text>',
                            info_x, cy))
    }
  }
  c(out, "</svg>") |> paste(collapse = "\n")
}

.specReportRawSort <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    out <- as.character(x)
    out[is.na(x)] <- ""
    out[is.infinite(x) & x > 0] <- "1e308"
    out[is.infinite(x) & x < 0] <- "-1e308"
    return(out)
  }
  out <- as.character(x)
  out[is.na(out)] <- ""
  out
}

.specReportFormatColumn <- function(x, name = "") {
  if (is.logical(x)) {
    out <- ifelse(is.na(x), "\u2014", ifelse(x, "yes", "no"))
  } else if (is.integer(x)) {
    out <- ifelse(is.na(x), "\u2014", format(x, trim = TRUE, scientific = FALSE))
  } else if (is.numeric(x)) {
    digits <- if (grepl("share|score|ratio|access|n_eff", name)) 4L else 3L
    out <- ifelse(is.na(x), "\u2014",
                  ifelse(is.infinite(x), ifelse(x > 0, "Inf", "-Inf"),
                         formatC(x, digits = digits, format = "fg", flag = "#")))
  } else {
    out <- as.character(x)
    out[is.na(out) | !nzchar(out)] <- "\u2014"
  }
  .specReportEscape(out)
}

.specReportTd <- function(value, name = "", content = NULL, class = NULL) {
  raw <- .specReportRawSort(value)
  if (is.null(content)) content <- .specReportFormatColumn(value, name)
  cls <- if (is.null(class)) "" else paste0(' class="', .specReportEscape(class), '"')
  paste0('<td data-sort="', .specReportEscape(raw), '"', cls, '>', content, '</td>')
}

.specReportLabel <- function(x) {
  labels <- c(
    signal = "Signal", coverage = "Coverage", n_variants = "Credible variants",
    n_variants_in_peaks = "Variants in peaks", n_peaks = "Variant peaks",
    variant_ids = "Variant IDs in peaks",
    n_high_cell = "High peaks (cell)", cell_category = "Cell category",
    sole_cell = "Sole cell contributor", top_cell = "Top cell",
    top_cell_share = "Top cell mean share", cell_has_enriched = "Enriched cell hit",
    n_high_lineage = "High peaks (lineage)", lineage_category = "Lineage category",
    sole_lineage = "Sole lineage contributor", top_lineage = "Top lineage",
    top_lineage_share = "Top lineage mean share",
    lineage_has_enriched = "Enriched lineage hit", nearest_genes = "Nearest genes",
    peak = "Peak", locus = "Signal", gwas_snps = "Variant(s)",
    top_score = "Top cell score", top_share = "Top cell share",
    is_high = "Cell high", class = "Cell class", second_cell = "Second cell",
    second_score = "Second score", ratio = "Top/second ratio",
    is_dominant_r20 = "Dominant (r>2)",
    is_dominant_r20_lineage = "Dominant lineage (r>2)",
    call_agreement = "High / dominant agreement",
    call_agreement_lineage = "High / dominant agreement (lineage)",
    top_score_lineage = "Top lineage score", top_share_lineage = "Top lineage share",
    is_high_lineage = "Lineage high", class_lineage = "Lineage class",
    second_lineage = "Second lineage", ratio_lineage = "Top/second ratio (lineage)",
    resolution_class = "Specificity resolution", access_max = "Max accessibility",
    access_mean = "Mean accessibility", dist_to_tss = "Distance to TSS",
    tss_bin = "TSS bin", nearest_gene = "Nearest gene",
    high_cells = "High cells", moderate_cells = "Moderate cells",
    high_lineages = "High lineages"
  )
  out <- unname(labels[x])
  out[is.na(out)] <- gsub("_", " ", x[is.na(out)], fixed = TRUE)
  out
}

.specReportCategoryBadge <- function(value) {
  val <- as.character(value)
  if (is.na(val) || !nzchar(val)) val <- "not_available"
  paste0('<span class="badge badge-', .specReportEscape(val), '">',
         .specReportEscape(gsub("_", " ", val, fixed = TRUE)), '</span>')
}

.specReportSignalTable <- function(d) {
  cols <- c("signal", "coverage", "n_variants", "n_variants_in_peaks", "n_peaks",
            "variant_ids", "n_high_cell", "cell_category", "sole_cell", "top_cell",
            "top_cell_share", "n_high_lineage", "lineage_category",
            "sole_lineage", "top_lineage", "top_lineage_share", "nearest_genes")
  cols <- intersect(cols, names(d))
  header <- paste0('<th data-type="',
                   ifelse(vapply(d[cols], is.numeric, logical(1)), "number", "text"),
                   '">', .specReportEscape(.specReportLabel(cols)), '</th>',
                   collapse = "")
  rows <- character(nrow(d))
  for (i in seq_len(nrow(d))) {
    cells <- character(length(cols))
    for (j in seq_along(cols)) {
      nm <- cols[j]
      value <- d[[nm]][i]
      if (nm == "signal") {
        content <- paste0('<button type="button" class="signal-link" data-signal-index="',
                          d$signal_index[i], '">', .specReportDisplay(value), '</button>')
      } else if (nm %in% c("cell_category", "lineage_category", "coverage")) {
        content <- .specReportCategoryBadge(value)
      } else {
        content <- NULL
      }
      cells[j] <- .specReportTd(value, nm, content)
    }
    rows[i] <- paste0('<tr data-signal-index="', d$signal_index[i],
                      '" data-signal-name="', .specReportEscape(d$signal[i]),
                      '" data-cell-category="', .specReportEscape(d$cell_category[i]),
                      '">', paste(cells, collapse = ""), '</tr>')
  }
  paste0('<table id="signal-table" class="data-table sortable"><thead><tr>',
         header, '</tr></thead><tbody>', paste(rows, collapse = "\n"),
         '</tbody></table>')
}

.specReportSignalDetails <- function(d) {
  if (!nrow(d)) return('<div class="empty-state">No signal metadata is available.</div>')
  rows <- character(nrow(d))
  for (i in seq_len(nrow(d))) {
    rows[i] <- paste0(
      '<article class="signal-detail-card" data-signal-index="', d$signal_index[i], '" hidden>',
      '<div class="detail-heading"><div><p class="eyebrow">Selected signal</p><h3>',
      .specReportDisplay(d$signal[i]), '</h3></div><div class="badge-pair">',
      .specReportCategoryBadge(d$cell_category[i]),
      .specReportCategoryBadge(d$lineage_category[i]), '</div></div>',
      '<div class="mini-grid">',
      '<div><span>Coverage</span><strong>', .specReportDisplay(d$coverage[i]), '</strong></div>',
      '<div><span>Credible variants</span><strong>', d$n_variants[i], '</strong></div>',
      '<div><span>Variants in retained peaks</span><strong>',
      .specReportDisplay(d$variant_ids[i]), '</strong></div>',
      '<div><span>Retained variant peaks</span><strong>', d$n_peaks[i], '</strong></div>',
      '<div><span>Cell high peaks</span><strong>', d$n_high_cell[i], '</strong></div>',
      '<div><span>Lineage high peaks</span><strong>', d$n_high_lineage[i], '</strong></div>',
      '<div><span>Nearest genes</span><strong>', .specReportDisplay(d$nearest_genes[i]), '</strong></div>',
      '</div></article>')
  }
  paste(rows, collapse = "\n")
}

.specReportPeakColumns <- function(d) {
  wanted <- c("locus", "peak", "gwas_snps", "top_cell", "top_score", "top_share",
              "is_high", "class", "second_cell", "second_score", "ratio",
              "is_dominant_r20", "high_cells", "moderate_cells",
              "top_lineage", "top_score_lineage", "top_share_lineage",
              "is_high_lineage", "class_lineage", "second_lineage",
              "second_score_lineage", "ratio_lineage", "is_dominant_r20_lineage",
              "high_lineages", "moderate_lineages", "call_agreement",
              "call_agreement_lineage", "resolution_class",
              "access_max", "access_mean", "nearest_gene",
              "dist_to_tss", "tss_bin")
  intersect(wanted, names(d))
}

.specReportSignalButtons <- function(signals, indices) {
  if (!length(signals)) return("\u2014")
  paste(vapply(seq_along(signals), function(i) paste0(
    '<button type="button" class="signal-link compact" data-signal-index="',
    indices[i], '">', .specReportEscape(signals[i]), '</button>'), character(1)),
    collapse = " ")
}

.specReportPeakTable <- function(data) {
  d <- data$peak_summary
  cols <- .specReportPeakColumns(d)
  if (!nrow(d)) return('<div class="empty-state">No retained variant peaks.</div>')
  header <- paste0('<th data-type="',
                   ifelse(vapply(d[cols], is.numeric, logical(1)), "number", "text"),
                   '">', .specReportEscape(.specReportLabel(cols)), '</th>',
                   collapse = "")
  cell_share <- data$cell_share
  lin_share <- data$lineage_share
  rows <- character(nrow(d))
  for (i in seq_len(nrow(d))) {
    signals <- data$peak_signals[[i]]
    indices <- data$peak_signal_indices[[i]]
    cells <- character(length(cols))
    for (j in seq_along(cols)) {
      nm <- cols[j]
      value <- d[[nm]][i]
      content <- if (nm == "locus") .specReportSignalButtons(signals, indices) else NULL
      cells[j] <- .specReportTd(value, nm, content)
    }
    key_string <- paste0(",", paste(indices, collapse = ","), ",")
    cell_values <- paste(formatC(cell_share[i, ], digits = 7, format = "fg"),
                         collapse = ",")
    lin_values <- paste(formatC(lin_share[i, ], digits = 7, format = "fg"),
                        collapse = ",")
    rows[i] <- paste0(
      '<tr data-signal-keys="', key_string,
      '" data-peak-label="', .specReportEscape(d$peak[i]),
      '" data-variants="', .specReportEscape(d$gwas_snps[i]),
      '" data-cell-shares="', cell_values,
      '" data-lineage-shares="', lin_values, '">',
      paste(cells, collapse = ""), '</tr>')
  }
  paste0('<table id="peak-table" class="data-table sortable"><thead><tr>',
         header, '</tr></thead><tbody>', paste(rows, collapse = "\n"),
         '</tbody></table>')
}

.specReportMetadataCards <- function(data) {
  s <- data$input_summary
  get_one <- function(name, fallback = NA) {
    if (nrow(s) && name %in% names(s)) s[[name]][1] else fallback
  }
  values <- c(
    "Signals" = nrow(data$signal_summary),
    "Signals with retained peaks" = sum(data$signal_summary$coverage == "retained_peaks"),
    "Retained variant peaks" = nrow(data$peak_summary),
    "Cell types" = length(data$cell_ids),
    "Lineages" = length(data$lineage_ids),
    "Per-cell normalisation" = get_one("per_cell_mode", "not supplied"),
    "Depth normalisation" = get_one("depth", "not supplied"),
    "Accessibility filter" = if (identical(get_one("access_filter_mode"), "none")) "off" else
      paste0(get_one("access_filter_mode", "not supplied"),
             " (", get_one("access_filter_cut", "?"), ")")
  )
  paste0('<div class="summary-grid">', paste(vapply(seq_along(values), function(i) {
    paste0('<div class="summary-card"><span>', .specReportEscape(names(values)[i]),
           '</span><strong>', .specReportDisplay(values[i]), '</strong></div>')
  }, character(1)), collapse = ""), '</div>')
}

.specReportCss <- r"--{
:root{--ink:#172033;--muted:#64748b;--line:#dbe2ea;--paper:#fff;--wash:#f4f7fb;--accent:#315c88;--accent2:#0f766e;--danger:#b2182b;--shadow:0 10px 28px rgba(15,23,42,.08)}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--wash);color:var(--ink);font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;line-height:1.45}
header{background:linear-gradient(125deg,#10243d,#315c88);color:white;padding:36px max(24px,calc((100vw - 1500px)/2));box-shadow:var(--shadow)}header h1{margin:0 0 6px;font-size:clamp(28px,4vw,44px);letter-spacing:-.03em}header p{margin:0;color:#d8e5f2}.report-meta{margin-top:13px;font-size:13px;color:#bcd0e3}
main{max-width:1500px;margin:26px auto;padding:0 22px 60px}section{background:var(--paper);border:1px solid var(--line);border-radius:14px;padding:22px;margin:0 0 22px;box-shadow:0 2px 10px rgba(15,23,42,.035)}h2{margin:0 0 7px;font-size:23px}h3{margin:0}.section-note,.muted{color:var(--muted);font-size:14px}.eyebrow{text-transform:uppercase;letter-spacing:.12em;font-size:11px;font-weight:800;color:var(--accent);margin:0 0 5px}
.summary-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin-top:16px}.summary-card{border:1px solid var(--line);border-radius:10px;padding:14px;background:#fbfdff}.summary-card span{display:block;color:var(--muted);font-size:12px;margin-bottom:4px}.summary-card strong{font-size:20px;overflow-wrap:anywhere}
.tabs{display:flex;gap:8px;margin:17px 0 13px}.tab-button,.control-button{border:1px solid #b7c4d2;background:white;color:var(--ink);border-radius:8px;padding:8px 13px;font-weight:700;cursor:pointer}.tab-button.active,.control-button.primary{background:var(--accent);border-color:var(--accent);color:white}.tab-panel[hidden]{display:none}.plot-scroll{overflow:auto;border:1px solid var(--line);border-radius:10px;background:white}.combined-svg{display:block;min-width:850px;width:100%;height:auto}.svg-title{font-size:15px;font-weight:800;fill:#172033}.svg-note{font-size:10px;fill:#64748b}.svg-axis{font-size:9px;fill:#475569}.svg-row{font-size:10px;fill:#172033}.svg-value{font-size:9px;fill:#334155}
.controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:16px 0 12px}.controls input,.controls select{min-height:38px;border:1px solid #b7c4d2;border-radius:8px;padding:7px 10px;background:white}.controls input{min-width:min(430px,100%);flex:1}.count-label{color:var(--muted);font-size:13px}
.table-wrap{overflow:auto;max-height:720px;border:1px solid var(--line);border-radius:10px}.data-table{border-collapse:separate;border-spacing:0;width:max-content;min-width:100%;font-size:12px}.data-table th{position:sticky;top:0;z-index:2;background:#edf3f8;color:#24364b;text-align:left;white-space:nowrap;padding:10px 9px;border-bottom:1px solid #bfcbd8;cursor:pointer;user-select:none}.data-table th::after{content:" \u2195";color:#94a3b8}.data-table td{padding:8px 9px;border-bottom:1px solid #e7ecf1;vertical-align:top;max-width:310px;overflow-wrap:anywhere}.data-table tbody tr:nth-child(even){background:#fafcff}.data-table tbody tr:hover,.data-table tbody tr.selected{background:#eaf3fb}.signal-link{border:0;background:none;color:#1d5f96;text-decoration:underline;font:inherit;font-weight:750;padding:0;cursor:pointer;text-align:left}.signal-link.compact{display:inline-block;margin:0 4px 3px 0}.badge{display:inline-flex;border-radius:999px;padding:3px 8px;font-size:10px;font-weight:800;white-space:nowrap;background:#e8edf3;color:#475569}.badge-coherent{background:#d8f3ea;color:#08705c}.badge-complex{background:#ffe3dc;color:#a23723}.badge-single{background:#e3ecfa;color:#315c88}.badge-no_high{background:#eceff2;color:#5f6975}.badge-retained_peaks{background:#d8f3ea;color:#08705c}.badge-no_retained_peak,.badge-not_evaluated,.badge-not_available{background:#fff1cc;color:#805b00}
.data-table th::after{content:'\2195'}
.signal-detail-card{border:1px solid var(--line);border-radius:11px;padding:16px;background:#fbfdff}.detail-heading{display:flex;justify-content:space-between;gap:15px;align-items:flex-start}.badge-pair{display:flex;gap:7px;flex-wrap:wrap}.mini-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:9px;margin-top:13px}.mini-grid div{border-left:3px solid #b9cadb;padding:6px 9px}.mini-grid span{display:block;color:var(--muted);font-size:11px}.mini-grid strong{font-size:13px;overflow-wrap:anywhere}
.contribution-toolbar{display:flex;gap:7px;margin:14px 0 10px}.contribution-list{display:grid;gap:8px}.contribution-row{display:grid;grid-template-columns:minmax(180px,310px) minmax(380px,1fr);gap:12px;align-items:center}.contribution-label{font-size:11px;overflow-wrap:anywhere}.contribution-label strong{display:block}.mass-bar{height:25px;display:flex;position:relative;overflow:hidden;border:1px solid #aab7c5;border-radius:5px;background:#f1f5f9}.mass-bar::after{content:"";position:absolute;left:50%;top:0;bottom:0;border-left:1px dashed rgba(15,23,42,.75);pointer-events:none}.mass-segment{height:100%;min-width:0}.legend{display:flex;flex-wrap:wrap;gap:7px 13px;margin-top:13px;font-size:10px;color:#475569}.legend-item{display:inline-flex;align-items:center;gap:5px}.swatch{width:9px;height:9px;border-radius:2px;display:inline-block}.empty-state{padding:18px;border:1px dashed #b8c4d0;border-radius:9px;color:var(--muted);background:#fbfdff}.filter-banner{display:none;align-items:center;gap:8px;background:#eaf3fb;border:1px solid #bad2e7;border-radius:8px;padding:9px 11px;margin:10px 0}.filter-banner.active{display:flex}.methods-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(270px,1fr));gap:15px}.methods-grid h3{font-size:15px;margin-bottom:5px}.methods-grid p{font-size:13px;color:var(--muted);margin:0}
@media(max-width:760px){main{padding:0 10px 35px}section{padding:15px}.contribution-row{grid-template-columns:1fr}.mass-bar{min-width:520px}.contribution-list{overflow:auto}.detail-heading{display:block}.badge-pair{margin-top:9px}}
}--"

.specReportJs <- r"--{
(() => {
  "use strict";
  const unitMeta = window.CLEAR_REPORT_UNITS;
  const state = { signal: null, contributionLevel: "cell" };

  function cellValue(td, type) {
    const raw = td?.dataset.sort ?? "";
    if (type === "number") {
      if (raw === "") return null;
      const value = Number(raw);
      return Number.isNaN(value) ? null : value;
    }
    return raw.toLocaleLowerCase();
  }

  document.querySelectorAll("table.sortable").forEach(table => {
    table.querySelectorAll("thead th").forEach((th, column) => {
      th.addEventListener("click", () => {
        const type = th.dataset.type || "text";
        const ascending = th.dataset.direction !== "asc";
        table.querySelectorAll("thead th").forEach(x => delete x.dataset.direction);
        th.dataset.direction = ascending ? "asc" : "desc";
        const rows = Array.from(table.tBodies[0].rows);
        rows.sort((a, b) => {
          const av = cellValue(a.cells[column], type);
          const bv = cellValue(b.cells[column], type);
          if (av === null || av === "") return bv === null || bv === "" ? 0 : 1;
          if (bv === null || bv === "") return -1;
          if (av < bv) return ascending ? -1 : 1;
          if (av > bv) return ascending ? 1 : -1;
          return 0;
        });
        rows.forEach(row => table.tBodies[0].appendChild(row));
      });
    });
  });

  document.querySelectorAll(".tab-button[data-tab-group]").forEach(button => {
    button.addEventListener("click", () => {
      const group = button.dataset.tabGroup;
      document.querySelectorAll(`.tab-button[data-tab-group="${group}"]`).forEach(x =>
        x.classList.toggle("active", x === button));
      document.querySelectorAll(`.tab-panel[data-tab-group="${group}"]`).forEach(panel =>
        panel.hidden = panel.dataset.tabPanel !== button.dataset.tabTarget);
    });
  });

  function applySignalFilters() {
    const query = (document.getElementById("signal-search")?.value || "").toLocaleLowerCase();
    const category = document.getElementById("signal-category")?.value || "all";
    const rows = Array.from(document.querySelectorAll("#signal-table tbody tr"));
    let visible = 0;
    rows.forEach(row => {
      const okQuery = !query || row.textContent.toLocaleLowerCase().includes(query);
      const okCategory = category === "all" || row.dataset.cellCategory === category;
      row.hidden = !(okQuery && okCategory);
      if (!row.hidden) visible++;
    });
    const label = document.getElementById("signal-count");
    if (label) label.textContent = `${visible} of ${rows.length} signals`;
  }

  function peakMatchesSignal(row) {
    return state.signal === null || row.dataset.signalKeys.includes(`,${state.signal},`);
  }

  function applyPeakFilters() {
    const query = (document.getElementById("peak-search")?.value || "").toLocaleLowerCase();
    const rows = Array.from(document.querySelectorAll("#peak-table tbody tr"));
    let visible = 0;
    rows.forEach(row => {
      const show = peakMatchesSignal(row) &&
        (!query || row.textContent.toLocaleLowerCase().includes(query));
      row.hidden = !show;
      if (show) visible++;
    });
    const label = document.getElementById("peak-count");
    if (label) label.textContent = `${visible} of ${rows.length} variant peaks`;
  }

  function renderLegend(level) {
    const meta = unitMeta[level];
    const legend = document.getElementById("contribution-legend");
    if (!legend) return;
    legend.replaceChildren();
    meta.names.forEach((name, i) => {
      const item = document.createElement("span");
      item.className = "legend-item";
      const swatch = document.createElement("span");
      swatch.className = "swatch";
      swatch.style.backgroundColor = meta.colors[i];
      item.append(swatch, document.createTextNode(name));
      legend.appendChild(item);
    });
  }

  function renderContributions() {
    const target = document.getElementById("contribution-list");
    if (!target) return;
    target.replaceChildren();
    if (state.signal === null) {
      target.innerHTML = '<div class="empty-state">Select a signal to show its per-peak contribution bars.</div>';
      document.getElementById("contribution-legend")?.replaceChildren();
      return;
    }
    const rows = Array.from(document.querySelectorAll("#peak-table tbody tr"))
      .filter(peakMatchesSignal);
    if (!rows.length) {
      target.innerHTML = '<div class="empty-state">This signal has no retained variant peaks in the analysed peak universe.</div>';
      document.getElementById("contribution-legend")?.replaceChildren();
      return;
    }
    const meta = unitMeta[state.contributionLevel];
    rows.forEach(row => {
      const shares = (state.contributionLevel === "cell" ? row.dataset.cellShares :
        row.dataset.lineageShares).split(",").map(Number);
      const order = shares.map((value, index) => ({value, index}))
        .sort((a, b) => b.value - a.value);
      const wrapper = document.createElement("div");
      wrapper.className = "contribution-row";
      const label = document.createElement("div");
      label.className = "contribution-label";
      const strong = document.createElement("strong");
      strong.textContent = row.dataset.peakLabel;
      label.append(strong, document.createTextNode(row.dataset.variants || "No named variant"));
      const bar = document.createElement("div");
      bar.className = "mass-bar";
      order.forEach(({value, index}) => {
        if (!Number.isFinite(value) || value <= 0) return;
        const segment = document.createElement("span");
        segment.className = "mass-segment";
        segment.style.width = `${Math.max(0, value * 100)}%`;
        segment.style.backgroundColor = meta.colors[index];
        segment.title = `${meta.names[index]}: ${(value * 100).toFixed(2)}%`;
        bar.appendChild(segment);
      });
      wrapper.append(label, bar);
      target.appendChild(wrapper);
    });
    renderLegend(state.contributionLevel);
  }

  function showSignal(index, shouldScroll = true) {
    state.signal = Number(index);
    document.querySelectorAll("#signal-table tbody tr").forEach(row =>
      row.classList.toggle("selected", Number(row.dataset.signalIndex) === state.signal));
    document.querySelectorAll(".signal-detail-card").forEach(card =>
      card.hidden = Number(card.dataset.signalIndex) !== state.signal);
    const signalRow = document.querySelector(`#signal-table tbody tr[data-signal-index="${state.signal}"]`);
    const signalName = signalRow?.dataset.signalName || `signal ${state.signal}`;
    const banner = document.getElementById("peak-filter-banner");
    if (banner) {
      banner.classList.add("active");
      banner.querySelector("strong").textContent = signalName;
    }
    applyPeakFilters();
    renderContributions();
    try {
      history.replaceState(null, "", `#signal=${encodeURIComponent(signalName)}`);
    } catch (_) {
      location.hash = `signal=${encodeURIComponent(signalName)}`;
    }
    if (shouldScroll) document.getElementById("signal-detail")?.scrollIntoView({behavior: "smooth"});
  }

  function clearSignal() {
    state.signal = null;
    document.querySelectorAll("#signal-table tbody tr").forEach(row => row.classList.remove("selected"));
    document.querySelectorAll(".signal-detail-card").forEach(card => card.hidden = true);
    document.getElementById("peak-filter-banner")?.classList.remove("active");
    applyPeakFilters();
    renderContributions();
    try {
      history.replaceState(null, "", location.href.split("#")[0]);
    } catch (_) {
      location.hash = "";
    }
  }

  document.addEventListener("click", event => {
    const signal = event.target.closest(".signal-link");
    if (signal) {
      event.preventDefault();
      showSignal(signal.dataset.signalIndex);
    }
  });
  document.getElementById("clear-signal")?.addEventListener("click", clearSignal);
  document.getElementById("signal-search")?.addEventListener("input", applySignalFilters);
  document.getElementById("signal-category")?.addEventListener("change", applySignalFilters);
  document.getElementById("peak-search")?.addEventListener("input", applyPeakFilters);
  document.querySelectorAll("[data-contribution-level]").forEach(button => {
    button.addEventListener("click", () => {
      state.contributionLevel = button.dataset.contributionLevel;
      document.querySelectorAll("[data-contribution-level]").forEach(x =>
        x.classList.toggle("active", x === button));
      renderContributions();
    });
  });

  document.getElementById("download-peaks")?.addEventListener("click", () => {
    const table = document.getElementById("peak-table");
    if (!table) return;
    const visibleRows = Array.from(table.tBodies[0].rows).filter(row => !row.hidden);
    const lines = [Array.from(table.tHead.rows[0].cells).map(x => x.textContent.trim())]
      .concat(visibleRows.map(row => Array.from(row.cells).map(x => x.textContent.trim())));
    const csv = lines.map(row => row.map(value => `"${value.replaceAll('"', '""')}"`).join(",")).join("\n");
    const blob = new Blob([csv], {type: "text/csv;charset=utf-8"});
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = state.signal === null ? "CLEAR_variant_peaks.csv" : `CLEAR_signal_${state.signal}_variant_peaks.csv`;
    link.click();
    URL.revokeObjectURL(link.href);
  });

  applySignalFilters();
  applyPeakFilters();
  renderContributions();
  if (location.hash.startsWith("#signal=")) {
    const name = decodeURIComponent(location.hash.slice(8));
    const row = Array.from(document.querySelectorAll("#signal-table tbody tr"))
      .find(x => x.dataset.signalName === name);
    if (row) showSignal(row.dataset.signalIndex, false);
  }
})();
}--"

writeSpecHtmlReport <- function(
    se_name, trait_name, mat_spec, gwas_mat, se_gwas, snps,
    peak_cell, locus_categorisation = data.frame(),
    global_enrichment_cell = data.frame(), global_enrichment_lineage = data.frame(),
    lineage_obj, input_summary = data.frame(), output_file,
    mapping_file = NULL, signal_file = NULL) {

  if (missing(output_file) || !nzchar(output_file)) {
    stop("output_file is required.")
  }
  data <- prepareSpecHtmlReportData(
    se_name = se_name, trait_name = trait_name, mat_spec = mat_spec,
    gwas_mat = gwas_mat, se_gwas = se_gwas, snps = snps,
    peak_cell = peak_cell, locus_categorisation = locus_categorisation,
    global_enrichment_cell = global_enrichment_cell,
    global_enrichment_lineage = global_enrichment_lineage,
    lineage_obj = lineage_obj, input_summary = input_summary)

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(mapping_file) && nzchar(mapping_file)) {
    dir.create(dirname(mapping_file), recursive = TRUE, showWarnings = FALSE)
    .specWrite(data$map, mapping_file)
  }
  if (!is.null(signal_file) && nzchar(signal_file)) .specWrite(data$signal_summary, signal_file)

  cell_colours_named <- setNames(data$cell_colours, data$cell_ids)
  lineage_colours_named <- setNames(data$lineage_colours, data$lineage_ids)
  cell_svg <- .specReportCombinedSvg(
    data$cell_correlation, data$global_enrichment_cell, "celltype",
    cell_colours_named, "Cell-type")
  lineage_svg <- .specReportCombinedSvg(
    data$lineage_correlation, data$global_enrichment_lineage, "lineage",
    lineage_colours_named, "Lineage")

  units_js <- paste0(
    "window.CLEAR_REPORT_UNITS={cell:{names:",
    .specReportJsEncodedArray(data$cell_ids), ",colors:",
    .specReportJsEncodedArray(data$cell_colours), "},lineage:{names:",
    .specReportJsEncodedArray(data$lineage_ids), ",colors:",
    .specReportJsEncodedArray(data$lineage_colours), "}};")

  signal_note <- if (isTRUE(data$signal_metadata_present)) "" else paste0(
    '<div class="empty-state"><strong>No signal column was supplied.</strong> ',
    'Variant overlaps are grouped as <em>unassigned</em>; signal-level navigation ',
    'requires <code>mcols(snps)$signal</code>.</div>')

  html <- c(
    '<!doctype html>', '<html lang="en">', '<head>',
    '<meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1">',
    paste0('<title>CLEAR specificity report: ', .specReportEscape(trait_name),
           ' \u00d7 ', .specReportEscape(se_name), '</title>'),
    '<style>', .specReportCss, '</style>', '</head>', '<body>',
    '<header>', '<h1>CLEAR specificity report</h1>',
    paste0('<p>', .specReportEscape(trait_name), ' <span aria-hidden="true">\u00d7</span> ',
           .specReportEscape(se_name), '</p>'),
    paste0('<div class="report-meta">Generated ',
           .specReportEscape(data$generated_at),
           ' \u00b7 static offline report \u00b7 high specificity: s &gt; 1/\u221a2 (s\u00b2 &gt; 50%)</div>'),
    '</header>', '<main>',
    '<section id="overview"><p class="eyebrow">Run overview</p><h2>Analysis coverage</h2>',
    '<p class="section-note">Counts describe the analysed CLEAR peak universe; normalisation and filtering settings are shown below. A signal with no retained peak is distinct from a retained signal whose peaks have no high-specificity call.</p>',
    .specReportMetadataCards(data), '</section>',

    '<section id="landscape"><p class="eyebrow">Specificity landscape</p>',
    '<h2>Correlation and global enrichment</h2>',
    '<p class="section-note">The heatmap uses genome-wide specificity across retained peaks. The aligned bars show the signed global rank-enrichment z-score for GWAS-overlapping peaks; hover over cells and bars for exact values.</p>',
    '<div class="tabs"><button type="button" class="tab-button active" data-tab-group="landscape" data-tab-target="cell">Cell type</button>',
    '<button type="button" class="tab-button" data-tab-group="landscape" data-tab-target="lineage">Lineage</button></div>',
    '<div class="tab-panel" data-tab-group="landscape" data-tab-panel="cell"><div class="plot-scroll">',
    cell_svg, '</div></div>',
    '<div class="tab-panel" data-tab-group="landscape" data-tab-panel="lineage" hidden><div class="plot-scroll">',
    lineage_svg, '</div></div>', '</section>',

    '<section id="signals"><p class="eyebrow">Signal explorer</p><h2>Per-signal summary and category</h2>',
    '<p class="section-note">Cell and lineage categories are separate annotations: single = one high peak; coherent = at least two high peaks with one contributor; complex = high peaks with multiple contributors; no high = retained peaks but no contributor above 50% mass.</p>',
    signal_note,
    '<div class="controls"><input id="signal-search" type="search" placeholder="Search signal, variant, gene or contributor..." aria-label="Search signals">',
    '<select id="signal-category" aria-label="Filter cell category"><option value="all">All cell categories</option><option value="coherent">Coherent</option><option value="complex">Complex</option><option value="single">Single</option><option value="no_high">No high</option><option value="not_evaluated">Not evaluated</option></select>',
    '<span id="signal-count" class="count-label"></span></div>',
    '<div class="table-wrap">', .specReportSignalTable(data$signal_summary), '</div>', '</section>',

    '<section id="signal-detail"><p class="eyebrow">Signal detail</p><h2>Variant-peak contribution profile</h2>',
    '<p class="section-note">Each bar is one retained variant peak. Segments are s\u00b2 accessibility-mass shares. Non-zero profiles sum to 100%; zero-signal peaks have no contribution. Segments are ordered largest first and the dashed marker is 50%.</p>',
    .specReportSignalDetails(data$signal_summary),
    '<div class="contribution-toolbar"><button type="button" class="tab-button active" data-contribution-level="cell">Cell type</button>',
    '<button type="button" class="tab-button" data-contribution-level="lineage">Lineage</button></div>',
    '<div id="contribution-list" class="contribution-list"></div><div id="contribution-legend" class="legend"></div>',
    '</section>',

    '<section id="variant-peaks"><p class="eyebrow">Complete result table</p><h2>All retained variant peaks</h2>',
    '<p class="section-note">This table always contains the full peak result set. Selecting a signal applies an exact signal-membership filter; text search can then narrow by variant, peak, gene or any displayed value.</p>',
    '<div id="peak-filter-banner" class="filter-banner">Showing signal: <strong></strong><button id="clear-signal" type="button" class="control-button">Show all peaks</button></div>',
    '<div class="controls"><input id="peak-search" type="search" placeholder="Search variant, peak, signal, gene, cell type..." aria-label="Search variant peaks">',
    '<button id="download-peaks" type="button" class="control-button">Download visible CSV</button><span id="peak-count" class="count-label"></span></div>',
    '<div class="table-wrap">', .specReportPeakTable(data), '</div>', '</section>',

    '<section id="methods"><p class="eyebrow">Interpretation</p><h2>How to read this report</h2>',
    '<div class="methods-grid"><div><h3>Specificity is relative</h3><p>s\u00b2 is the share of a peak\'s accessibility mass assigned to a cell type or lineage. Read it beside max/mean accessibility; a large share is not itself strong accessibility.</p></div>',
    '<div><h3>Categories are descriptive</h3><p>Single, coherent, complex and no-high describe retained variant peaks. They are not locus significance or causal-probability rankings.</p></div>',
    '<div><h3>Cell and lineage calls differ</h3><p>Related cell types can split mass at cell resolution and pool into a lineage-specific call. Both resolutions are therefore retained throughout.</p></div>',
    '<div><h3>Nearest gene is proximity</h3><p>The nearest-gene and TSS fields provide genomic context; they do not identify a regulatory target gene.</p></div></div>',
    '</section>', '</main>',
    '<script>', units_js, '</script>', '<script>', .specReportJs, '</script>',
    '</body>', '</html>'
  )

  tmp <- paste0(output_file, ".", Sys.getpid(), ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(html, tmp, useBytes = TRUE)
  if (!file.rename(tmp, output_file)) {
    if (!file.copy(tmp, output_file, overwrite = TRUE)) {
      stop("Could not write HTML report: ", output_file)
    }
    unlink(tmp)
  }
  path <- normalizePath(output_file, mustWork = FALSE)
  message("  HTML report: ", path)
  path
}
