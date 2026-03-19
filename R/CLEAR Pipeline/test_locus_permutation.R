#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(dplyr)
  library(rtracklayer)
})

# Set to TRUE to run directly in RStudio with assigned parameters below.
use_rstudio_params <- FALSE

rstudio_params <- list(
  se_name = "",
  se_path = "",
  trait_name = "",
  trait_path = "",
  out_dir = "",
  gtf_path = "/working/lab_jonathb/lefeiW/projects/CLEAR_data/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz",
  n_perm = 100L,
  k_lineage = 4L,
  target_locus = "",
  seed = 1L
)

args <- commandArgs(trailingOnly = TRUE)

if (use_rstudio_params) {
  se_name <- rstudio_params$se_name
  se_path <- rstudio_params$se_path
  trait_name <- rstudio_params$trait_name
  trait_path <- rstudio_params$trait_path
  out_dir <- rstudio_params$out_dir
  gtf_path <- rstudio_params$gtf_path
  n_perm <- as.integer(rstudio_params$n_perm)
  k_lineage <- as.integer(rstudio_params$k_lineage)
  target_locus <- as.character(rstudio_params$target_locus)
  seed <- as.integer(rstudio_params$seed)

  required <- c("se_name", "se_path", "trait_name", "trait_path", "out_dir")
  vals <- c(se_name, se_path, trait_name, trait_path, out_dir)
  if (any(!nzchar(vals))) {
    stop(
      "RStudio mode is enabled. Please fill rstudio_params for: ",
      paste(required[!nzchar(vals)], collapse = ", ")
    )
  }
} else {
  if (length(args) < 5) {
    stop(
      paste0(
        "Usage: Rscript test_locus_permutation.R ",
        "<se_name> <se_path> <trait_name> <trait_path> <out_dir> ",
        "[gtf_path] [n_perm] [k_lineage] [target_locus] [seed]\n",
        "Or set use_rstudio_params <- TRUE and fill rstudio_params."
      )
    )
  }

  se_name <- args[1]
  se_path <- args[2]
  trait_name <- args[3]
  trait_path <- args[4]
  out_dir <- args[5]
  gtf_path <- if (length(args) >= 6) args[6] else "/working/lab_jonathb/lefeiW/projects/CLEAR_data/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz"
  n_perm <- if (length(args) >= 7) as.integer(args[7]) else 100L
  k_lineage <- if (length(args) >= 8) as.integer(args[8]) else 4L
  target_locus <- if (length(args) >= 9) args[9] else ""
  seed <- if (length(args) >= 10) as.integer(args[10]) else 1L
}

set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("\n========================================")
message("Locus permutation test")
message("  SE name:      ", se_name)
message("  Trait name:   ", trait_name)
message("  SE path:      ", se_path)
message("  Trait path:   ", trait_path)
message("  Out dir:      ", out_dir)
message("  GTF path:     ", gtf_path)
message("  n_perm:       ", n_perm)
message("  k_lineage:    ", k_lineage)
message("  target_locus: ", ifelse(nzchar(target_locus), target_locus, "<all>"))
message("  seed:         ", seed)
message("========================================\n")

prepare_clear_se <- function(se, genome = "hg38") {
  rd <- as.data.frame(rowData(se))
  has_rowdata_coords <- all(c("seqnames", "start", "end") %in% names(rd))

  if (has_rowdata_coords) {
    gr <- GRanges(
      seqnames = rd$seqnames,
      ranges = IRanges::IRanges(start = as.integer(rd$start), end = as.integer(rd$end))
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
    stop("Input SE has no usable genomic coordinates in rowData or rowRanges.")
  }

  se <- SummarizedExperiment(
    assays = as.list(assays(se)),
    rowRanges = gr,
    colData = colData(se)
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

compute_raw_l2 <- function(se, assay_name = NULL) {
  if (is.null(assay_name)) assay_name <- assayNames(se)[1]

  mat <- assay(se, assay_name)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  rn <- rownames(se)
  cn <- colnames(se)
  dimnames(mat) <- list(rn, cn)

  raw_l2 <- t(apply(mat, 1, safe_l2))
  dimnames(raw_l2) <- list(rn, cn)
  raw_l2
}

get_tss <- function(gtf) {
  genes <- gtf[gtf$type == "gene"]
  genes <- genes[genes$gene_type %in% c("lncRNA", "protein_coding")]
  tss_gr <- GenomicRanges::promoters(genes, upstream = 0, downstream = 1)
  tss_gr
}

compute_peak_density <- function(peaks_gr, window = 50000) {
  counts <- countOverlaps(peaks_gr, resize(peaks_gr, width = window, fix = "center"))
  counts - 1
}

make_density_bins <- function(density_vec, nbins = 4) {
  cut(
    density_vec,
    breaks = unique(quantile(density_vec, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)),
    include.lowest = TRUE
  )
}

compute_gwas_overlap <- function(se, snps) {
  hits_se <- findOverlaps(snps, rowRanges(se))
  gwas_idx <- unique(subjectHits(hits_se))

  if (!length(gwas_idx)) {
    se_gwas <- se[0, , drop = FALSE]
  } else {
    se_gwas <- se[gwas_idx, , drop = FALSE]

    snp_name_col <- if ("names" %in% colnames(mcols(snps))) mcols(snps)$names else names(snps)
    signal_col <- if ("signal" %in% colnames(mcols(snps))) mcols(snps)$signal else rep(NA_character_, length(snps))

    ann <- data.frame(
      peak_idx = subjectHits(hits_se),
      CCVs = snp_name_col[queryHits(hits_se)],
      signal = signal_col[queryHits(hits_se)],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::group_by(peak_idx) %>%
      dplyr::summarise(
        CCVs = paste(unique(na.omit(CCVs)), collapse = ","),
        signal = paste(unique(na.omit(signal)), collapse = ","),
        .groups = "drop"
      )

    rd <- as.data.frame(rowData(se_gwas))
    if (!"CCVs" %in% names(rd)) rd$CCVs <- ""
    if (!"signal" %in% names(rd)) rd$signal <- ""

    pos <- match(ann$peak_idx, gwas_idx)
    rd$CCVs[pos] <- ann$CCVs
    rd$signal[pos] <- ann$signal
    rowData(se_gwas) <- S4Vectors::DataFrame(rd, row.names = rownames(se_gwas))
  }

  list(
    hits_se = hits_se,
    gwas_idx = gwas_idx,
    se_gwas = se_gwas
  )
}

compute_locus_rank_concordance_verbose <- function(
  gwas_mat,
  full_mat,
  peak_to_locus,
  joint_bin_full,
  cell_lineage,
  n_perm = 100,
  target_locus = ""
) {
  rank_gwas <- t(apply(gwas_mat, 1, rank))
  rank_full <- t(apply(full_mat, 1, rank))

  valid_locus <- peak_to_locus[rownames(gwas_mat)]
  keep <- !is.na(valid_locus) & nzchar(valid_locus)
  if (!any(keep)) {
    return(list(summary = data.frame(), perm_long = data.frame()))
  }

  gwas_peaks <- rownames(gwas_mat)[keep]
  locus_split <- split(gwas_peaks, valid_locus[keep])

  if (nzchar(target_locus)) {
    locus_split <- locus_split[names(locus_split) == target_locus]
    if (!length(locus_split)) {
      stop("target_locus not found among GWAS loci: ", target_locus)
    }
  }

  full_bins <- split(names(joint_bin_full), joint_bin_full)

  summary_list <- list()
  perm_rows <- list()
  row_i <- 1L

  for (locus in names(locus_split)) {
    peaks <- locus_split[[locus]]
    k <- length(peaks)

    if (k < 2) {
      summary_list[[row_i]] <- data.frame(
        locus = locus,
        n_peaks = k,
        concordance = NA_real_,
        perm_mean = NA_real_,
        perm_sd = NA_real_,
        z_score = NA_real_,
        p_value = NA_real_,
        dominant_cell_locus = NA_character_,
        dominant_score_locus = NA_real_,
        dominant_lineage = NA_character_,
        lineage_fraction = NA_real_,
        n_perm_eff = 0L,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
      next
    }

    obs_mat <- rank_gwas[peaks, , drop = FALSE]
    obs_cor <- cor(t(obs_mat), method = "spearman")
    obs_conc <- mean(obs_cor[upper.tri(obs_cor)])

    mean_profile <- colMeans(gwas_mat[peaks, , drop = FALSE])
    dominant_cell_locus <- names(which.max(mean_profile))
    dominant_score_locus <- max(mean_profile)

    lineage_profile <- tapply(mean_profile, cell_lineage[names(mean_profile)], sum)
    dominant_lineage <- names(which.max(lineage_profile))
    lineage_fraction <- max(lineage_profile) / sum(lineage_profile)

    locus_joint_bins <- joint_bin_full[peaks]
    locus_bins <- table(locus_joint_bins)

    if (any(!names(locus_bins) %in% names(full_bins))) {
      summary_list[[row_i]] <- data.frame(
        locus = locus,
        n_peaks = k,
        concordance = obs_conc,
        perm_mean = NA_real_,
        perm_sd = NA_real_,
        z_score = NA_real_,
        p_value = NA_real_,
        dominant_cell_locus = dominant_cell_locus,
        dominant_score_locus = dominant_score_locus,
        dominant_lineage = dominant_lineage,
        lineage_fraction = lineage_fraction,
        n_perm_eff = 0L,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
      next
    }

    sampled_matrix <- matrix(NA_character_, nrow = k, ncol = n_perm)
    row_index <- 1L

    for (bin in names(locus_bins)) {
      n_bin <- as.integer(locus_bins[[bin]])
      pool <- full_bins[[bin]]
      if (is.null(pool) || length(pool) == 0) next

      draws <- replicate(
        n_perm,
        sample(pool, n_bin, replace = length(pool) < n_bin),
        simplify = "matrix"
      )

      sampled_matrix[row_index:(row_index + n_bin - 1), ] <- draws
      row_index <- row_index + n_bin
    }

    perm_vals <- apply(sampled_matrix, 2, function(sampled_peaks) {
      perm_mat <- rank_full[sampled_peaks, , drop = FALSE]
      perm_cor <- cor(t(perm_mat), method = "spearman")
      mean(perm_cor[upper.tri(perm_cor)])
    })

    perm_mean <- mean(perm_vals, na.rm = TRUE)
    perm_sd <- sd(perm_vals, na.rm = TRUE)
    n_perm_eff <- sum(!is.na(perm_vals))
    z_score <- ifelse(is.na(perm_sd) || perm_sd == 0, NA_real_, (obs_conc - perm_mean) / perm_sd)
    p_value <- if (n_perm_eff > 0) {
      (sum(perm_vals >= obs_conc, na.rm = TRUE) + 1) / (n_perm_eff + 1)
    } else {
      NA_real_
    }

    summary_list[[row_i]] <- data.frame(
      locus = locus,
      n_peaks = k,
      concordance = obs_conc,
      perm_mean = perm_mean,
      perm_sd = perm_sd,
      z_score = z_score,
      p_value = p_value,
      dominant_cell_locus = dominant_cell_locus,
      dominant_score_locus = dominant_score_locus,
      dominant_lineage = dominant_lineage,
      lineage_fraction = lineage_fraction,
      n_perm_eff = n_perm_eff,
      stringsAsFactors = FALSE
    )

    perm_rows[[locus]] <- data.frame(
      locus = locus,
      perm_id = seq_along(perm_vals),
      perm_concordance = as.numeric(perm_vals),
      observed_concordance = obs_conc,
      stringsAsFactors = FALSE
    )

    row_i <- row_i + 1L
  }

  summary_df <- dplyr::bind_rows(summary_list)
  perm_df <- dplyr::bind_rows(perm_rows)
  list(summary = summary_df, perm_long = perm_df)
}

if (!file.exists(se_path)) stop("SE file not found: ", se_path)
if (!file.exists(trait_path)) stop("Trait file not found: ", trait_path)
if (!file.exists(gtf_path)) stop("GTF file not found: ", gtf_path)

se <- readRDS(se_path)
snps <- readRDS(trait_path)
se <- prepare_clear_se(se)
mat_l2 <- compute_raw_l2(se)

ov <- compute_gwas_overlap(se, snps)
se_gwas <- ov$se_gwas
if (!nrow(se_gwas)) stop("No GWAS-overlapping peaks found.")

if (!"signal" %in% colnames(rowData(se_gwas))) {
  stop("No signal column found in GWAS-overlap rowData; cannot define loci.")
}

gwas_mat <- mat_l2[rownames(se_gwas), , drop = FALSE]
peak_to_locus <- setNames(as.character(rowData(se_gwas)$signal), rownames(se_gwas))

gtf <- rtracklayer::import(gtf_path)
tss_gr <- get_tss(gtf)

peaks_gr <- rowRanges(se)
nearest_hits <- distanceToNearest(peaks_gr, tss_gr)
dist2tss <- mcols(nearest_hits)$distance

tss_bin <- cut(
  dist2tss,
  breaks = c(0, 1000, 5000, 50000, 200000, Inf),
  labels = c("core_prom", "prox_prom", "near_reg", "distal", "long_range"),
  include.lowest = TRUE
)
names(tss_bin) <- rownames(mat_l2)

density <- compute_peak_density(rowRanges(se))
density_bin <- make_density_bins(density)
names(density_bin) <- rownames(mat_l2)

joint_bin_full <- interaction(tss_bin, density_bin, drop = TRUE)
joint_bin_full <- as.character(joint_bin_full)
names(joint_bin_full) <- names(tss_bin)

cell_cor <- cor(mat_l2, method = "pearson")
hc <- hclust(as.dist(1 - cell_cor), method = "average")
lineage_ids <- cutree(hc, k = k_lineage)
lineage_labels <- tapply(names(lineage_ids), lineage_ids, function(x) paste(sort(x), collapse = ":"))
cell_lineage <- setNames(as.character(lineage_labels[as.character(lineage_ids)]), names(lineage_ids))

res <- compute_locus_rank_concordance_verbose(
  gwas_mat = gwas_mat,
  full_mat = mat_l2,
  peak_to_locus = peak_to_locus,
  joint_bin_full = joint_bin_full,
  cell_lineage = cell_lineage,
  n_perm = n_perm,
  target_locus = target_locus
)

summary_out <- file.path(out_dir, paste0(se_name, "_", trait_name, "_locus_concordance_permutation_test.txt"))
perm_out <- file.path(out_dir, paste0(se_name, "_", trait_name, "_locus_permutation_values_long.txt"))

write.table(res$summary, summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res$perm_long, perm_out, sep = "\t", quote = FALSE, row.names = FALSE)

message("Saved summary: ", summary_out)
message("Saved permutation values: ", perm_out)
message("Done.")
