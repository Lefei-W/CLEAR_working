
# Minimal, fast CLEAR exploration script with stepwise QC, optional grouping, metrics, and GWAS overlap helpers.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(SummarizedExperiment)
  library(GenomicRanges); library(S4Vectors); library(Matrix); library(preprocessCore)
  library(readr)
})

safe_l2 <- function(x) { # on the rows (peaks)
  s <- sum(x^2)
  if (s == 0) x else x / sqrt(s)  # l2 euclidean normalization avoiding row sum of 0 (happen with the Human-scATAC-Corpus)
}

add_rank_desc <- function(M) {
  # rank 1 = highest value (descending)
  apply(M, 2, function(x) rank(-x, ties.method = "average"))
}

ensure_peakids <- function(se) {
  rr <- rowRanges(se)
  if (is.null(rownames(se)) || anyNA(rownames(se))) {
    rownames(se) <- paste0(as.character(seqnames(rr)), ":", start(rr), "-", end(rr))
  }
  if (!"peakid" %in% names(rowData(se))) rowData(se)$peakid <- rownames(se) # Add RowData 
  se
}

compute_metrics_se <- function(se, assay_name = NULL, keep_top_prop = NULL) {
  stopifnot(inherits(se, "SummarizedExperiment"))
  if (is.null(assay_name)) assay_name <- assayNames(se)[1]
  stopifnot(assay_name %in% assayNames(se))

  # --- load matrix ---
  mat <- assay(se, assay_name)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  rn <- rownames(se); if (is.null(rn)) rn <- paste0("peak_", seq_len(nrow(mat)))
  cn <- colnames(se); if (is.null(cn)) cn <- paste0("ct_",   seq_len(ncol(mat)))
  dimnames(mat) <- list(rn, cn)

  # --- optional filtering by row sum ---
  row_sums <- rowSums(mat)
  if (!is.null(keep_top_prop)) {
    stopifnot(is.numeric(keep_top_prop), keep_top_prop > 0, keep_top_prop <= 1)
    thr  <- stats::quantile(row_sums, probs = 1 - keep_top_prop, na.rm = TRUE)
    keep <- row_sums > thr
    if (!all(keep)) {
      se  <- se[keep, , drop = FALSE]
      mat <- mat[keep, , drop = FALSE]
      rn  <- rn[keep]
      row_sums <- row_sums[keep]
      dimnames(mat) <- list(rn, cn)
    }
  }

  # --- transforms ---
  mat_qn <- preprocessCore::normalize.quantiles(mat)
  dimnames(mat_qn) <- list(rn, cn)

  raw_l2 <- t(apply(mat,    1, safe_l2)); dimnames(raw_l2) <- list(rn, cn)
  qn_l2  <- t(apply(mat_qn, 1, safe_l2)); dimnames(qn_l2)  <- list(rn, cn)

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

  covW_raw    <- cov_weight(mat)
  covW_qn     <- cov_weight(mat_qn)
  covW_raw_l2 <- cov_weight(raw_l2)
  covW_qn_l2  <- cov_weight(qn_l2)

  corW_raw    <- cor_weight(mat)
  corW_qn     <- cor_weight(mat_qn)
  corW_raw_l2 <- cor_weight(raw_l2)
  corW_qn_l2  <- cor_weight(qn_l2)

  # --- save assays (keep original intact) ---
  assays(se)[["CLEAR_raw"]]         <- mat
  assays(se)[["CLEAR_qn"]]          <- mat_qn
  assays(se)[["CLEAR_raw_l2"]]      <- raw_l2
  assays(se)[["CLEAR_qn_l2"]]       <- qn_l2

  assays(se)[["CLEAR_raw_l2_rank"]] <- add_rank_desc(raw_l2)
  assays(se)[["CLEAR_qn_l2_rank"]]  <- add_rank_desc(qn_l2)

  assays(se)[["CLEAR_covW_raw"]]    <- covW_raw
  assays(se)[["CLEAR_covW_qn"]]     <- covW_qn
  assays(se)[["CLEAR_covW_raw_l2"]] <- covW_raw_l2
  assays(se)[["CLEAR_covW_qn_l2"]]  <- covW_qn_l2

  assays(se)[["CLEAR_corW_raw"]]    <- corW_raw
  assays(se)[["CLEAR_corW_qn"]]     <- corW_qn
  assays(se)[["CLEAR_corW_raw_l2"]] <- corW_raw_l2
  assays(se)[["CLEAR_corW_qn_l2"]]  <- corW_qn_l2

  metadata(se)$CLEAR <- list(
    source_assay   = assay_name,
    keep_top_prop  = keep_top_prop,
    n_peaks        = nrow(se),
    n_celltypes    = ncol(se),
    rowSum_summary = c(min = min(row_sums), median = stats::median(row_sums), max = max(row_sums))
  )

  se
}

# Group highly correlated cell types into clusters and average their signal
group_celltypes <- function(se, assay_name = NULL, cor_cutoff = 0.9) {
  if (is.null(assay_name)) assay_name <- assayNames(se)[1]
  stopifnot(assay_name %in% assayNames(se))
  mat <- SummarizedExperiment::assay(se, assay_name)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  
  cor_mat <- stats::cor(mat)
  hc      <- stats::hclust(as.dist(1 - cor_mat))
  clusters <- stats::cutree(hc, h = 1 - cor_cutoff)
  
  grouped_mat <- sapply(unique(clusters), function(cl) {
    cols <- names(clusters[clusters == cl])
    if (length(cols) == 1) mat[, cols] else rowMeans(mat[, cols, drop = FALSE])
  })
  grouped_mat <- as.matrix(grouped_mat)
  colnames(grouped_mat) <- sapply(unique(clusters), function(cl) paste(names(clusters[clusters == cl]), collapse = ":"))
  rownames(grouped_mat) <- rownames(mat)
  
  se_grouped <- SummarizedExperiment(
    assays    = setNames(list(grouped_mat), assay_name),
    rowRanges = rowRanges(se),
    rowData   = rowData(se)
  )
  se_grouped <- ensure_peakids(se_grouped)
  
  metadata(se_grouped)$CLEAR <- c(
    metadata(se)$CLEAR,
    list(
      grouping      = TRUE,
      cor_cutoff    = cor_cutoff,
      cluster_sizes = as.integer(table(clusters)),
      clusters      = clusters
    )
  )
  
  list(se = se_grouped, cor_mat = cor_mat, clusters = clusters, hclust = hc)
}


# QC for one assay (correlation heatmap with legend, dendrogram, max-value histogram)
plot_qc_full_matrix <- function(se, tissue, assay, outfile = NULL) {
  stopifnot(assay %in% assayNames(se))
  
  mat <- SummarizedExperiment::assay(se, assay)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  
  cor_mat <- stats::cor(mat)
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  names(cor_df) <- c("celltype1", "celltype2", "cor")
  
  max_vals <- matrixStats::rowMaxs(mat)
  
  hc <- stats::hclust(as.dist(1 - cor_mat))
  cell_order <- hc$labels[hc$order]
  
  if (!is.null(outfile)) grDevices::pdf(outfile, width = 10, height = 8)
  on.exit(if (!is.null(outfile)) grDevices::dev.off(), add = TRUE)
  
  cor_df <- cor_df %>%
    mutate(
      celltype1 = factor(celltype1, levels = cell_order),
      celltype2 = factor(celltype2, levels = cell_order)
    )
  p_cor <- ggplot(cor_df, aes(x = celltype1, y = celltype2, fill = cor)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("blue", "white", "red"),
      values  = scales::rescale(c(0, 0.5, 1)),  # anchor points
      limits  = c(0, 1),                         # enforce mapping
      name    = "corr"
    ) +
    labs(
      title = paste0(tissue, " | Correlation (", assay, ")"),
      x = "", y = ""
    ) +
    coord_fixed() +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7))
  
  print(p_cor)
  
  graphics::plot(hc, main = paste0(tissue, " | Hierarchical clustering (", assay, ")"))
  
  p_hist <- ggplot(data.frame(max_val = max_vals), aes(x = max_val)) +
    geom_histogram(bins = 80, fill = "gray60", color = "white") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = paste0(tissue, " | Max value per peak (", assay, ")"), x = "max value", y = "count") +
    theme_minimal()
  
  print(p_hist)
  
  invisible(list(cor_mat = cor_mat, dendro = hc, max_vals = max_vals))
}

#  GWAS overlap + QC 
overlap_with_credset <- function(se, credset_gr) {
  stopifnot(inherits(se, "SummarizedExperiment"), inherits(credset_gr, "GenomicRanges"))
  se <- ensure_peakids(se)
  ov <- findOverlaps(se, credset_gr, ignore.strand = TRUE)
  hit_rows <- sort(unique(queryHits(ov)))
  se_sub <- se[hit_rows, , drop = FALSE]
  
  # best-effort CCV label
  cs_cols <- names(mcols(credset_gr))
  pick <- cs_cols[cs_cols %in% c("names","SNP","rsid","rsID","variant","marker","ID")]
  rsid <- if (length(pick)) mcols(credset_gr)[[pick[1]]] else rep(NA_character_, length(credset_gr))
  
  kk <- data.frame(orig_row = queryHits(ov), rsid = rsid[subjectHits(ov)], stringsAsFactors = FALSE) |>
    dplyr::group_by(orig_row) |>
    dplyr::summarise(CCVs = paste(unique(na.omit(rsid)), collapse = ","), .groups = "drop") # collapse multiple CCVs per peak into a comma-separated string
  
  rd <- as.data.frame(rowData(se_sub))
  rd$CCVs <- ""
  rd$CCVs[match(kk$orig_row, hit_rows)] <- kk$CCVs
  rowData(se_sub) <- S4Vectors::DataFrame(rd, row.names = rownames(se_sub))
  
  metadata(se_sub)$CLEAR <- c(metadata(se)$CLEAR, list(N_total = nrow(se)))
  se_sub
}

meanrank_from_overlap <- function(se_overlap, rank_assay) {
  if (!rank_assay %in% assayNames(se_overlap)) return(NULL)
  md <- metadata(se_overlap)$CLEAR
  N <- md$N_total
  M <- SummarizedExperiment::assay(se_overlap, rank_assay)
  if (inherits(M, "Matrix")) M <- as.matrix(M)
  n <- nrow(M)
  if (n == 0) return(NULL)
  
  mu    <- (N + 1) / 2
  sigma <- sqrt((N^2 - 1) / (12 * n))
  mr    <- colMeans(M)
  
  out <- data.frame(
    celltype           = colnames(M),
    observed_mean_rank = as.numeric(mr),
    n_peaks            = n,
    exp_rank           = mu,
    stringsAsFactors   = FALSE
  )
  out$obs.exp  <- out$observed_mean_rank / mu
  out$p        <- pnorm(out$observed_mean_rank, mu, sigma, lower.tail = TRUE)
  out$fdr      <- p.adjust(out$p, method = "BH")
  out
}

plot_meanrank <- function(se_overlap, rank_assay, trait, outfile) {
  mr <- meanrank_from_overlap(se_overlap, rank_assay)
  if (is.null(mr)) return(invisible(NULL))
  thresh <- -log10(0.05 / max(1, nrow(mr)))
  ttl <- paste0(trait, " | ", rank_assay, " mean-rank")
  p <- ggplot(mr) +
    geom_hline(yintercept = thresh, linetype = "dashed", color = "gray40") +
    geom_col(aes(x = reorder(celltype, -log10(p)), y = -log10(p), fill = (-log10(p) >= thresh))) +
    scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "firebrick")) +
    coord_flip() +
    labs(title = ttl, x = "Cell type", y = "-log10(p)") +
    theme_minimal() + theme(legend.position = "none")
  if (!missing(outfile) && !is.null(outfile)) {
    grDevices::pdf(outfile, width = 8, height = 6); print(p); grDevices::dev.off()
  }
  invisible(mr)
}

plot_ccv_max_l2 <- function(se_overlap, metric = "qn_l2", trait, outfile) {
  stopifnot(metric %in% assayNames(se_overlap))
  mat <- SummarizedExperiment::assay(se_overlap, metric)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  max_l2 <- matrixStats::rowMaxs(mat)
  p <- ggplot(data.frame(max_l2 = max_l2), aes(x = max_l2)) +
    geom_histogram(bins = 80, fill = "gray60", color = "white") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = paste0(trait, " | Max ", metric, " per CCV peak"), x = "max L2", y = "count") +
    theme_minimal()
  if (!missing(outfile) && !is.null(outfile)) {
    grDevices::pdf(outfile, width = 6, height = 5); print(p); grDevices::dev.off()
  }
  invisible(max_l2)
}


# ---------- IO helpers ----------
save_full_matrix <- function(se, tissue, out_root = "results/full_matrix") {
  md <- metadata(se)$CLEAR
  fname <- paste0(
    tissue, "_", md$n_celltypes, "ct_", md$n_peaks, "peaks_",
    "CLEAR_transformed_l2_rank.se.rds"
  )
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_root, fname)
  saveRDS(se, path)
  path
}

save_overlap_matrix <- function(se_overlap, tissue, trait, out_root = "results/GWAS_overlaps") {
  md     <- metadata(se_overlap)$CLEAR
  n_cred <- length(md$credset_gr)
  ccv_col <- rowData(se_overlap)$CCVs
  n_ccv  <- dplyr::n_distinct(unlist(strsplit(paste(ccv_col, collapse = ","), ",")))
  
  dir1 <- file.path(out_root, paste0(trait, "_", n_cred, "cred_", n_ccv, "ccv"))
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)
  
  fname <- paste0(
    trait, "_", tissue, "_", md$n_celltypes, "ct_",
    nrow(se_overlap), "peaks_CLEAR_transformed_l2_rank.se.rds"
  )
  path <- file.path(dir1, fname)
  saveRDS(se_overlap, path)
  path
}
# ---------- Example usage (uncomment to run) ----------
# setwd("/working/lab_jonathb/lefeiW/projects/CLEAR_breast_demo/")
# source("~/plottheme.R")
# heart <- readRDS("results/snATAC_SE_library/kanemaru2023_heart_ArchR/PeakMatrix/SE_kanemaru2023_heart_ArchR_PeakMatrix.hg38.rds")
# heart <- ensure_peakids(heart)
#
# # QC on raw (choose assay if multiple)
# plot_qc_full_matrix(heart, tissue = "heart", assay = assayNames(heart)[1], outfile = "results/qc_full/heart_raw_qc.pdf")
#
# # Optional: group cell types
# grp <- group_celltypes(heart, cor_cutoff = 0.9)
# heart_grp <- grp$se
#
# # Compute metrics (optionally keep_top_prop)
# heart_met <- compute_metrics_se(heart_grp, keep_top_prop = 0.9)
# save_full_matrix(heart_met, tissue = "heart_grp")
#
# # GWAS overlap + mean-rank + max-L2
# credset <- readRDS("../ATAC_BCAC/data/opentarget_trait/heart_failure_2025_lee_cs165__GCST90455657/heart_failure_2025_lee_cs165__GCST90455657_CLEAR_credible_set_variants.gr.rds")
# GenomeInfoDb::seqlevelsStyle(credset) <- "UCSC"
# ov <- overlap_with_credset(heart_met, credset)
# plot_meanrank(ov, rank_assay = "qn_l2_rank_desc", trait = "HF_cs165", outfile = "results/qc_gwas/HF_cs165_meanrank.pdf")
# plot_ccv_max_l2(ov, metric = "qn_l2", trait = "HF_cs165", outfile = "results/qc_gwas/HF_cs165_maxL2.pdf")
# save_overlap_matrix(ov, tissue = "heart_grp", trait = "HF_cs165")
