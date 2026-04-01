# ============================================================
# peak_annotations.R
# TSS distance, peak density, joint binning, locus definition
# ============================================================

get_tss <- function(gtf) {
  genes <- gtf[gtf$type == "gene"]
  genes <- genes[genes$gene_type %in% c("lncRNA", "protein_coding")]
  tss_gr <- GenomicRanges::promoters(genes, upstream = 0, downstream = 1)
  S4Vectors::mcols(tss_gr) <- S4Vectors::mcols(tss_gr)[, c("gene_type", "gene_name")]
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

prepare_peak_annotations <- function(se, gtf, gwas_mat, specificity_summary_GWAS = NULL,
                                     density_window = 50000,
                                     max_locus_distance = 100000) {
  peaks_gr <- rowRanges(se)
  peak_names <- names(peaks_gr)

  tss_gr <- get_tss(gtf)
  nearest_hits <- distanceToNearest(peaks_gr, tss_gr)

  dist_to_tss <- mcols(nearest_hits)$distance
  names(dist_to_tss) <- peak_names[queryHits(nearest_hits)]
  mcols(peaks_gr)$dist_to_tss <- dist_to_tss[peak_names]

  tss_bins <- cut(
    dist_to_tss,
    breaks = c(0, 1000, 5000, 50000, 200000, Inf),
    labels = c("core_prom", "prox_prom", "near_reg", "distal", "long_range"),
    include.lowest = TRUE
  )
  names(tss_bins) <- peak_names
  tss_bins_gwas <- tss_bins[rownames(gwas_mat)]

  density_vec <- compute_peak_density(peaks_gr, window = density_window)
  names(density_vec) <- peak_names

  density_bins <- cut(
    density_vec,
    breaks = unique(quantile(density_vec, probs = seq(0, 1, length.out = 10))),
    include.lowest = TRUE
  )
  names(density_bins) <- peak_names

  joint_bins <- interaction(tss_bins, density_bins, drop = TRUE)
  joint_bins <- as.character(joint_bins)
  names(joint_bins) <- peak_names

  peaks_sorted <- sort(peaks_gr)
  same_chr <- as.character(seqnames(peaks_sorted)[-1]) ==
    as.character(seqnames(peaks_sorted)[-length(peaks_sorted)])
  close_peaks <- diff(start(peaks_sorted)) <= max_locus_distance
  cluster_id <- cumsum(c(1, !(close_peaks & same_chr)))
  locus_ids <- paste0("locus_", cluster_id)
  names(locus_ids) <- names(peaks_sorted)
  locus_ids_full <- locus_ids[peak_names]
  peak_to_locus_full <- setNames(locus_ids_full, peak_names)

  if (!is.null(specificity_summary_GWAS) && "gwas_locus" %in% names(specificity_summary_GWAS)) {
    peak_to_locus_gwas <- specificity_summary_GWAS$gwas_locus
    names(peak_to_locus_gwas) <- specificity_summary_GWAS$peak
  } else if (!is.null(specificity_summary_GWAS) && "locus" %in% names(specificity_summary_GWAS)) {
    peak_to_locus_gwas <- specificity_summary_GWAS$locus
    names(peak_to_locus_gwas) <- specificity_summary_GWAS$peak
  } else {
    peak_to_locus_gwas <- setNames(as.character(rowData(se[rownames(gwas_mat),])$signal), rownames(gwas_mat))
  }

  peak_annot <- data.frame(
    peak = peak_names,
    dist_to_tss = dist_to_tss[peak_names],
    density = density_vec[peak_names],
    tss_bin = as.character(tss_bins[peak_names]),
    density_bin = as.character(density_bins[peak_names]),
    joint_bin = joint_bins[peak_names],
    locus = locus_ids_full,
    stringsAsFactors = FALSE
  )

  list(
    peak_annot = peak_annot,
    peak_to_locus_full = peak_to_locus_full,
    peak_to_locus_gwas = peak_to_locus_gwas,
    tss_bins_full = tss_bins,
    tss_bins_gwas = tss_bins_gwas,
    density_vec = density_vec,
    density_bins = density_bins,
    joint_bins = joint_bins
  )
}

add_gwas_peak_annotation <- function(summary_df, snps) {
  if (!nrow(summary_df)) return(summary_df)
  
  peak_gr <- GRanges(summary_df$peak)
  hits <- findOverlaps(peak_gr, snps)
  snps_names <- if ("names" %in% colnames(mcols(snps))) mcols(snps)$names else names(snps)
  signal_col <- if ("signal" %in% colnames(mcols(snps))) mcols(snps)$signal else rep(NA_character_, length(snps))
  
  snp_list <- tapply(snps_names[subjectHits(hits)], queryHits(hits), paste, collapse = ",")
  locus_vec <- tapply(signal_col[subjectHits(hits)], queryHits(hits), function(x) paste(unique(na.omit(x)), collapse = ", "))
  
  summary_df$gwas_snps <- NA_character_
  summary_df$gwas_locus <- NA_character_
  summary_df$gwas_snps[as.numeric(names(snp_list))] <- snp_list
  summary_df$gwas_locus[as.numeric(names(locus_vec))] <- locus_vec
  summary_df
}
