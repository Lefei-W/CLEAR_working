# ============================================================
# locus_concordance.R
# Locus-level rank concordance, consensus concordance,
# permutation testing, and global enrichment
# ============================================================

compute_locus_rank_concordance_vec <- function(
  gwas_mat,
  full_mat,
  peak_to_locus,
  joint_bin_full,
  cell_lineage,
  n_perm = 100) {
  
  rank_gwas <- t(apply(gwas_mat, 1, rank))
  rank_full <- t(apply(full_mat, 1, rank))
  
  valid_locus <- peak_to_locus[rownames(gwas_mat)]
  keep <- !is.na(valid_locus) & nzchar(valid_locus)
  if (!any(keep)) return(data.frame())
  
  gwas_peaks <- rownames(gwas_mat)[keep]
  locus_split <- split(gwas_peaks, valid_locus[keep])
  full_bins <- split(names(joint_bin_full), joint_bin_full)
  
  results <- lapply(names(locus_split), function(locus) {
    peaks <- locus_split[[locus]]
    k <- length(peaks)
    
    if (k < 2) {
      return(data.frame(
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
        stringsAsFactors = FALSE
      ))
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
      return(data.frame(
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
        stringsAsFactors = FALSE
      ))
    }
    
    sampled_matrix <- matrix(NA_character_, nrow = k, ncol = n_perm)
    row_index <- 1
    
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
    
    data.frame(
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
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, results)
}

compute_locus_consensus_concordance_vec <- function(
    gwas_mat,
    full_mat,
    l2_weights,
    peak_to_locus,
    peak_to_locus_full,
    joint_bin_full,
    cell_lineage,
    n_perm = 100
) {
  locus_to_peaks_full <- split(names(peak_to_locus_full), peak_to_locus_full)
  all_bins <- sort(unique(joint_bin_full))

  message("precomputing bin counts per locus (consensus)...")
  locus_bin_mat <- t(sapply(locus_to_peaks_full, function(peaks) {
    tab <- table(factor(joint_bin_full[peaks], levels = all_bins))
    as.numeric(tab)
  }))
  rownames(locus_bin_mat) <- names(locus_to_peaks_full)
  locus_bin_prop <- locus_bin_mat / rowSums(locus_bin_mat)

  gwas_loci_full <- unique(peak_to_locus_full[rownames(gwas_mat)])
  non_gwas_loci <- setdiff(names(locus_to_peaks_full), gwas_loci_full)

  locus_split <- split(rownames(gwas_mat), peak_to_locus[rownames(gwas_mat)])

  group_labels <- setNames(colnames(gwas_mat), colnames(gwas_mat))

  results <- lapply(names(locus_split), function(locus) {
    message(paste0("calculating consensus coherence for ", locus, "..."))
    peaks <- locus_split[[locus]]
    k <- length(peaks)

    if (k < 2) {
      return(data.frame(
        locus = locus, n_peaks = k,
        consensus_fine = NA_real_, consensus_lineage = NA_real_,
        perm_mean = NA_real_, perm_sd = NA_real_,
        z_score = NA_real_, p_value = NA_real_,
        dominant_cell_consensus = NA_character_,
        dominant_lineage_consensus = NA_character_,
        stringsAsFactors = FALSE))
    }

    fine_obs <- dominance_consensus(gwas_mat[peaks, , drop = FALSE], l2_weights[peaks], group_labels)
    fine_obs_conc <- fine_obs$consensus
    fine_top_cell <- fine_obs$top_group

    lineage_obs <- dominance_consensus(gwas_mat[peaks, , drop = FALSE], l2_weights[peaks], cell_lineage)
    lineage_obs_conc <- lineage_obs$consensus
    lineage_top_grp  <- lineage_obs$top_group

    locus_joint_bins <- joint_bin_full[peaks]

    perm_vals <- replicate(n_perm, {
      valid_loci <- non_gwas_loci[sapply(locus_to_peaks_full[non_gwas_loci], length) >= k]
      if (length(valid_loci) == 0) return(NA)

      lb <- table(factor(locus_joint_bins, levels = all_bins))
      lb <- as.numeric(lb) / sum(lb)
      candidate_mat <- locus_bin_prop[valid_loci, , drop = FALSE]
      dists <- rowSums(abs(candidate_mat - matrix(lb, nrow = nrow(candidate_mat), ncol = length(lb), byrow = TRUE)))

      top_n <- min(20, length(valid_loci))
      best_loci <- valid_loci[order(dists)][seq_len(top_n)]
      sampled_locus <- sample(best_loci, 1)
      sampled_peaks <- sample(locus_to_peaks_full[[sampled_locus]], k)

      m <- full_mat[sampled_peaks, , drop = FALSE]
      dominance_consensus(m, l2_weights[sampled_peaks], group_labels)$consensus
    })

    perm_mean <- mean(perm_vals, na.rm = TRUE)
    perm_sd   <- sd(perm_vals, na.rm = TRUE)
    n_perm_eff <- sum(!is.na(perm_vals))
    z_score <- if (is.na(perm_sd) || perm_sd == 0) NA_real_ else (fine_obs_conc - perm_mean) / perm_sd
    p_value <- if (n_perm_eff > 0) sum(perm_vals >= fine_obs_conc, na.rm = TRUE) / n_perm_eff else NA_real_

    data.frame(
      locus = locus, n_peaks = k,
      consensus_fine = fine_obs_conc,
      consensus_lineage = lineage_obs_conc,
      perm_mean = perm_mean, perm_sd = perm_sd,
      z_score = z_score, p_value = p_value,
      dominant_cell_consensus = fine_top_cell,
      dominant_lineage_consensus = lineage_top_grp,
      stringsAsFactors = FALSE)
  })

  do.call(rbind, results)
}

addLocusCoherence <- function(se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir) {
  if (!nrow(gwas_mat)) return(data.frame())
  if (!file.exists(gtf_path)) {
    message("  Skipping locus analysis: GTF file not found at ", gtf_path)
    return(data.frame())
  }
  
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

  cell_lineage <- lineage_obj$cell_lineage
  
  locus_results <- compute_locus_rank_concordance_vec(
    gwas_mat = gwas_mat,
    full_mat = mat_l2,
    peak_to_locus = peak_to_locus,
    joint_bin_full = joint_bin_full,
    cell_lineage = cell_lineage,
    n_perm = n_perm
  )
  
  if (!nrow(locus_results)) return(data.frame())
  
  locus_results$significant <- locus_results$p_value < 0.05
  locus_results$category <- dplyr::case_when(
    locus_results$n_peaks == 1 ~ "Single_peak",
    locus_results$z_score > 2 & locus_results$p_value < 0.05 & locus_results$dominant_score_locus > 0.5 ~ "Highly_coherent_specific",
    locus_results$z_score > 2 & locus_results$p_value < 0.05 ~ "Coherent_moderate",
    locus_results$z_score < 1 ~ "Incoherent",
    TRUE ~ "Intermediate"
  )
  
  write.table(locus_results, file.path(results_dir, "locus_concordance_permutation.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  locus_results
}

addLocusCoherence_consensus <- function(se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir, specificity_summary_l2 = NULL) {
  if (!nrow(gwas_mat)) return(data.frame())
  if (!file.exists(gtf_path)) {
    message("  Skipping consensus locus analysis: GTF file not found at ", gtf_path)
    return(data.frame())
  }

  gtf <- rtracklayer::import(gtf_path)

  peak_annotation <- prepare_peak_annotations(
    se = se, gtf = gtf, gwas_mat = gwas_mat,
    specificity_summary_GWAS = specificity_summary_l2
  )

  peak_to_locus_gwas <- peak_annotation$peak_to_locus_gwas
  if (is.null(peak_to_locus_gwas) || all(is.na(peak_to_locus_gwas))) {
    peak_to_locus_gwas <- setNames(as.character(rowData(se_gwas)$signal), rownames(se_gwas))
  }

  l2_weights <- compute_l2_weights(se)

  cell_lineage <- lineage_obj$cell_lineage

  locus_results_consensus <- compute_locus_consensus_concordance_vec(
    gwas_mat = gwas_mat,
    full_mat = mat_l2,
    l2_weights = l2_weights,
    peak_to_locus = peak_to_locus_gwas,
    peak_to_locus_full = peak_annotation$peak_to_locus_full,
    joint_bin_full = peak_annotation$joint_bins,
    cell_lineage = cell_lineage,
    n_perm = n_perm
  )

  if (!nrow(locus_results_consensus)) return(data.frame())

  locus_results_consensus$significant <- locus_results_consensus$p_value < 0.05
  locus_results_consensus$category <- dplyr::case_when(
    locus_results_consensus$n_peaks == 1 ~ "Single_peak",
    locus_results_consensus$z_score > 2 & locus_results_consensus$p_value < 0.05 &
      locus_results_consensus$consensus_lineage > 0.5 ~ "Highly_coherent_specific",
    locus_results_consensus$z_score > 2 & locus_results_consensus$p_value < 0.05 ~ "Coherent_moderate",
    locus_results_consensus$z_score < 1 ~ "Incoherent",
    TRUE ~ "Intermediate"
  )

  write.table(locus_results_consensus, file.path(results_dir, "locus_concordance_consensus_permutation.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  locus_results_consensus
}

addGlobalEnrichment <- function(se, snps, mat_l2, results_dir, plots_dir, se_name, trait_name) {
  nmat_r <- apply(mat_l2, 2, function(x) rank(x, ties.method = "min"))
  nse_rank <- SummarizedExperiment(nmat_r, rowRanges = rowRanges(se))
  gwas_rank_idx <- unique(subjectHits(findOverlaps(snps, rowRanges(nse_rank))))
  gwas_rank_mat <- if (length(gwas_rank_idx)) assay(nse_rank[gwas_rank_idx, , drop = FALSE], 1) else nmat_r[0, , drop = FALSE]
  
  if (!nrow(gwas_rank_mat)) return(data.frame())
  
  N <- nrow(nse_rank)
  n <- nrow(gwas_rank_mat)
  exp_mean <- (N + 1) / 2
  se_z <- sqrt((N^2 - 1) / (12 * n))
  
  global_enrichment <- data.frame(observed_mean_rank = colMeans(gwas_rank_mat), stringsAsFactors = FALSE)
  global_enrichment$n_peaks <- n
  global_enrichment$exp_rank <- exp_mean
  global_enrichment$obs.exp <- global_enrichment$observed_mean_rank / exp_mean
  global_enrichment$z <- (global_enrichment$observed_mean_rank - exp_mean) / se_z
  global_enrichment$p <- pnorm(global_enrichment$z, lower.tail = FALSE)
  global_enrichment$fdr <- p.adjust(global_enrichment$p, method = "BH")
  global_enrichment$positive <- global_enrichment$z > 0 & global_enrichment$fdr < 0.05
  global_enrichment$celltype <- rownames(global_enrichment)
  
  # FDR significance stars
  global_enrichment$fdr_stars <- dplyr::case_when(
    is.na(global_enrichment$fdr) ~ "",
    global_enrichment$fdr < 0.001 ~ "***",
    global_enrichment$fdr < 0.01  ~ "**",
    global_enrichment$fdr < 0.05  ~ "*",
    TRUE ~ ""
  )

  write.table(global_enrichment, file.path(results_dir, "global_enrichment.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  p_enrich <- ggplot(global_enrichment, aes(x = reorder(celltype, z), y = z, fill = positive)) +
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    geom_text(aes(label = fdr_stars), hjust = ifelse(global_enrichment$z >= 0, -0.3, 1.3), size = 4) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    labs(
      x = "Cell type", y = "Z-score (rank enrichment)",
      title = paste0(se_name, " x ", trait_name, " - Global GWAS Enrichment"),
      subtitle = "* FDR<0.05  ** FDR<0.01  *** FDR<0.001",
      fill = "FDR < 0.05"
    )
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_global_enrichment.pdf")), p_enrich, width = 8, height = 6)
  
  p_fdr <- ggplot(global_enrichment, aes(x = reorder(celltype, -log10(fdr)), y = -log10(fdr), fill = positive)) +
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.6) +
    geom_text(aes(label = fdr_stars), hjust = -0.3, size = 4) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      x = "Cell type", y = "-log10(FDR)",
      title = paste0(se_name, " x ", trait_name, " - Global Enrichment (FDR)"),
      subtitle = "* FDR<0.05  ** FDR<0.01  *** FDR<0.001",
      fill = "FDR < 0.05"
    )
  ggsave(file.path(plots_dir, paste0(se_name, "_", trait_name, "_global_enrichment_fdr.pdf")), p_fdr, width = 8, height = 6)
  
  global_enrichment
}
