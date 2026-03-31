# ============================================================
# locus_concordance.R
# Locus-level rank concordance, consensus concordance,
# permutation testing, and global enrichment
# ============================================================
#
# This file contains two complementary locus-level coherence methods:
#
# 1. RANK CONCORDANCE (compute_locus_rank_concordance_vec)

# 2. CONSENSUS CONCORDANCE (compute_locus_consensus_concordance_vec)

# Both methods use bin-matched permutation testing to generate null distributions:

#
# Wrapper functions (addLocusCoherence, addLocusCoherence_consensus) handle
# peak annotation, locus assignment, and result categorisation.
# ============================================================


# ============================================================
# RANK CONCORDANCE
# ============================================================
# Measures whether GWAS-overlapping peaks within a locus show
# consistent cell-type specificity RANKINGS. 
#
# Metric: mean pairwise Spearman rho across peaks in the locus
# Null: bin-matched random peaks from the genome-wide matrix
# ============================================================
compute_locus_rank_concordance_vec <- function(
  gwas_mat,            # L2 GWAS-overlapping peaks 
  full_mat,            # genome-wide (used for permutation null)
  peak_to_locus,       # locus ID (for GWAS peaks)
  joint_bin_full,      # joint bin label (TSS distance x density, genome-wide)
  cell_lineage,        # lineage group label
  n_perm = 1000) {     # number of permutations for null 

  # Rank each peak's specificity values across cell types (row-wise ranking)
  # Not the same as global enrichment ranking; this is per-peak to capture relative specificity patterns
  rank_gwas <- t(apply(gwas_mat, 1, rank))
  rank_full <- t(apply(full_mat, 1, rank))

  # this can take advantage of the peak-level summary with the dominancy 
  # Map GWAS peaks to their locus (OpenTargets uses the lead SNP)
  valid_locus <- peak_to_locus[rownames(gwas_mat)]
  keep <- !is.na(valid_locus) & nzchar(valid_locus)
  if (!any(keep)) return(data.frame())

  # Group GWAS peaks by locus for per-locus analysis
  gwas_peaks <- rownames(gwas_mat)[keep]
  locus_split <- split(gwas_peaks, valid_locus[keep])

  # genome-wide peaks by joint bin for efficient permutation sampling
  full_bins <- split(names(joint_bin_full), joint_bin_full)
  
  results <- lapply(names(locus_split), function(locus) {
    peaks <- locus_split[[locus]]
    k <- length(peaks)

    # Single-peak loci NAs
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

    # --- Observed concordance ---
    # Subset the ranked matrix to this locus's peaks, then compute
    # Spearman correlation between all pairs of peaks (transposed so
    # columns = peaks). The mean of the upper triangle is the locus concordance.
    obs_mat <- rank_gwas[peaks, , drop = FALSE]
    obs_cor <- cor(t(obs_mat), method = "spearman")
    obs_conc <- mean(obs_cor[upper.tri(obs_cor)])

    # --- Dominant cell type annotation ---
    # Average the raw L2 specificity across locus peaks to identify
    # which cell type has the strongest mean specificity at this locus
    mean_profile <- colMeans(gwas_mat[peaks, , drop = FALSE])
    dominant_cell_locus <- names(which.max(mean_profile))
    dominant_score_locus <- max(mean_profile)

    # --- Lineage-level aggregation ---
    # Sum cell-type specificities within each lineage group;
    # the dominant lineage is the group capturing the most total specificity
    lineage_profile <- tapply(mean_profile, cell_lineage[names(mean_profile)], sum)
    dominant_lineage <- names(which.max(lineage_profile))
    lineage_fraction <- max(lineage_profile) / sum(lineage_profile)
    
    # --- Bin-matched permutation null ---
    # Each peak belongs to a joint bin (TSS distance x peak density).
    # The permutation preserves this bin composition: for each bin in the
    # observed locus, we sample the same number of peaks from that bin genome-wide.
    locus_joint_bins <- joint_bin_full[peaks]
    locus_bins <- table(locus_joint_bins)

    # If any bins are missing from the genome-wide pool, we cannot build
    # a proper null -- return observed concordance without permutation stats
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

    # Build a (k x n_perm) matrix of sampled peak names.
    # For each bin, draw n_bin peaks from the genome-wide pool of that bin,
    # repeated n_perm times. This preserves the promoter/distal and density
    # configuration of the observed locus in every permutation.
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

    # Compute concordance for each permuted "locus" the same way as observed:
    # rank -> pairwise Spearman -> mean upper triangle
    perm_vals <- apply(sampled_matrix, 2, function(sampled_peaks) {
      perm_mat <- rank_full[sampled_peaks, , drop = FALSE]
      perm_cor <- cor(t(perm_mat), method = "spearman")
      mean(perm_cor[upper.tri(perm_cor)])
    })

    # --- Z-score and empirical p-value ---
    # Z = (observed - null mean) / null SD
    # p = fraction of null values >= observed (one-sided, with +1 correction)
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

# ============================================================
# CONSENSUS CONCORDANCE
# ============================================================
# Measures whether GWAS-overlapping peaks within a locus agree
# on which SINGLE cell type (or lineage) is dominant.
#
# Uses dominance_consensus() (from dominance.R):
#   1. For each peak, sum L2 specificity scores within each group
#   2. Assign each peak to its highest-scoring group
#   3. Weight each peak's "vote" by (max group specificity) * (peak accessibility)
#   4. Tally weighted votes per group; the winner is the plurality group
#   5. Consensus score = winner's share of total weighted votes
#
# A score of 1.0 = all peaks point to the same group
# A score of 1/n_groups = votes spread evenly (no agreement)
#
# Permutation null: sample whole loci from non-GWAS genome-wide loci
# matched by bin composition (L1 distance on bin proportions)
# ============================================================
compute_locus_consensus_concordance_vec <- function(
    gwas_mat,            # L2-normalised specificity matrix for GWAS peaks (peaks x cell types)
    full_mat,            # L2-normalised specificity for ALL peaks genome-wide (for permutation)
    l2_weights,          # named vector: peak accessibility weights (rowMeans of raw matrix)
    peak_to_locus,       # named vector: GWAS peak name -> locus ID
    peak_to_locus_full,  # named vector: ALL peak names -> locus ID (genome-wide locus structure)
    joint_bin_full,      # named vector: peak name -> joint bin label (TSS distance x density)
    cell_lineage,        # named vector: cell type name -> lineage group label
    n_perm = 100         # number of permutations for null distribution
) {
  # Build a lookup from locus -> all peak names (genome-wide)
  locus_to_peaks_full <- split(names(peak_to_locus_full), peak_to_locus_full)
  all_bins <- sort(unique(joint_bin_full))

  # Precompute the bin composition (as a proportion vector) for every genome-wide locus.
  # This is used later to find loci with similar promoter/distal + density configurations.
  message("precomputing bin counts per locus (consensus)...")
  locus_bin_mat <- t(sapply(locus_to_peaks_full, function(peaks) {
    tab <- table(factor(joint_bin_full[peaks], levels = all_bins))
    as.numeric(tab)
  }))
  rownames(locus_bin_mat) <- names(locus_to_peaks_full)
  locus_bin_prop <- locus_bin_mat / rowSums(locus_bin_mat)

  # Exclude GWAS-overlapping loci from the permutation pool to avoid circularity
  gwas_loci_full <- unique(peak_to_locus_full[rownames(gwas_mat)])
  non_gwas_loci <- setdiff(names(locus_to_peaks_full), gwas_loci_full)

  # Group GWAS peaks by their locus assignment
  locus_split <- split(rownames(gwas_mat), peak_to_locus[rownames(gwas_mat)])

  # Fine-level group labels: each cell type is its own group (identity mapping)
  # This means consensus at the cell-type level, not the lineage level
  group_labels <- setNames(colnames(gwas_mat), colnames(gwas_mat))

  results <- lapply(names(locus_split), function(locus) {
    message(paste0("calculating consensus coherence for ", locus, "..."))
    peaks <- locus_split[[locus]]
    k <- length(peaks)

    # Single-peak loci: consensus is undefined (need >=2 peaks to measure agreement)
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

    # --- Observed consensus (cell-type level) ---
    # Each peak votes for its top cell type, weighted by specificity * accessibility.
    # consensus_fine = fraction of weighted votes going to the winning cell type.
    fine_obs <- dominance_consensus(gwas_mat[peaks, , drop = FALSE], l2_weights[peaks], group_labels)
    fine_obs_conc <- fine_obs$consensus
    fine_top_cell <- fine_obs$top_group

    # --- Observed consensus (lineage level) ---
    # Same logic but groups are lineage clusters instead of individual cell types.
    # This is a coarser measure - easier to achieve high scores.
    lineage_obs <- dominance_consensus(gwas_mat[peaks, , drop = FALSE], l2_weights[peaks], cell_lineage)
    lineage_obs_conc <- lineage_obs$consensus
    lineage_top_grp  <- lineage_obs$top_group

    # --- Locus-matched permutation null ---
    # Unlike rank concordance (which samples individual peaks by bin), consensus
    # permutation samples WHOLE LOCI from the genome, matched by bin composition.
    # Steps per permutation:
    #   1. Find non-GWAS loci with >= k peaks
    #   2. Compute L1 distance between observed locus bin proportions and each candidate
    #   3. Pick one of the top-20 closest-matching loci at random
    #   4. Sample k peaks from that locus
    #   5. Compute consensus on the sampled peaks using the genome-wide L2 matrix
    locus_joint_bins <- joint_bin_full[peaks]

    perm_vals <- replicate(n_perm, {
      # Filter to loci with enough peaks to match the observed locus size
      valid_loci <- non_gwas_loci[sapply(locus_to_peaks_full[non_gwas_loci], length) >= k]
      if (length(valid_loci) == 0) return(NA)

      # Compute bin proportion vector for the observed locus
      lb <- table(factor(locus_joint_bins, levels = all_bins))
      lb <- as.numeric(lb) / sum(lb)

      # L1 distance between observed and each candidate locus's bin composition
      candidate_mat <- locus_bin_prop[valid_loci, , drop = FALSE]
      dists <- rowSums(abs(candidate_mat - matrix(lb, nrow = nrow(candidate_mat), ncol = length(lb), byrow = TRUE)))

      # Select from the top-20 best-matching loci (introduces some stochasticity
      # while still controlling for promoter/distal + density configuration)
      top_n <- min(20, length(valid_loci))
      best_loci <- valid_loci[order(dists)][seq_len(top_n)]
      sampled_locus <- sample(best_loci, 1)
      sampled_peaks <- sample(locus_to_peaks_full[[sampled_locus]], k)

      # Compute consensus on the genome-wide L2 matrix (not the GWAS subset)
      m <- full_mat[sampled_peaks, , drop = FALSE]
      dominance_consensus(m, l2_weights[sampled_peaks], group_labels)$consensus
    })

    # --- Z-score and empirical p-value ---
    # Z = (observed consensus - null mean) / null SD
    # p = fraction of null values >= observed (one-sided)
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

# ============================================================
# WRAPPER: addLocusCoherence (rank concordance)
# ============================================================
# High-level function that prepares peak annotations inline
# (TSS distance bins, peak density bins, joint bins) and calls
# compute_locus_rank_concordance_vec. Categorises results and
# writes output to results_dir.
# ============================================================
addLocusCoherence <- function(se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir) {
  if (!nrow(gwas_mat)) return(data.frame())
  if (!file.exists(gtf_path)) {
    message("  Skipping locus analysis: GTF file not found at ", gtf_path)
    return(data.frame())
  }

  # Locus assignment from GWAS SNP metadata (signal column = locus name)
  peak_to_locus <- setNames(as.character(rowData(se_gwas)$signal), rownames(se_gwas))

  # --- Build joint bins for permutation matching ---
  # TSS distance bins: classify each peak by distance to nearest gene TSS
  #   core_prom (0-1kb), prox_prom (1-5kb), near_reg (5-50kb),
  #   distal (50-200kb), long_range (>200kb)
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

  # Density bins: classify peaks by local peak density (number of peaks within 50kb window)
  density <- compute_peak_density(rowRanges(se))
  density_bin <- make_density_bins(density)
  names(density_bin) <- rownames(mat_l2)

  # Joint bin = interaction of TSS bin x density bin
  # This controls for the confounding effect of genomic context in the permutation null
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

  # --- Categorise loci ---
  # Highly_coherent_specific: significant z>2, p<0.05, and the dominant cell type is strong (>0.5)
  # Coherent_moderate: significant but dominant score is not exceptionally high
  # Incoherent: z < 1, peaks show no more agreement than random
  # Intermediate: between coherent and incoherent
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

# ============================================================
# WRAPPER: addLocusCoherence_consensus (consensus concordance)
# ============================================================
# High-level function that delegates peak annotation to
# prepare_peak_annotations() (from peak_annotations.R) and
# calls compute_locus_consensus_concordance_vec.
# Unlike addLocusCoherence which builds bins inline, this
# uses the shared annotation infrastructure and also needs
# genome-wide locus structure (peak_to_locus_full) for the
# whole-locus permutation strategy.
# ============================================================
addLocusCoherence_consensus <- function(se, se_gwas, mat_l2, gwas_mat, gtf_path, lineage_obj, n_perm, results_dir, specificity_summary_l2 = NULL) {
  if (!nrow(gwas_mat)) return(data.frame())
  if (!file.exists(gtf_path)) {
    message("  Skipping consensus locus analysis: GTF file not found at ", gtf_path)
    return(data.frame())
  }

  # Use prepare_peak_annotations() to compute TSS bins, density bins,
  # joint bins, and locus assignments in one shot
  gtf <- rtracklayer::import(gtf_path)

  peak_annotation <- prepare_peak_annotations(
    se = se, gtf = gtf, gwas_mat = gwas_mat,
    specificity_summary_GWAS = specificity_summary_l2
  )

  # Prefer GWAS locus labels from peak annotation; fall back to SE metadata
  peak_to_locus_gwas <- peak_annotation$peak_to_locus_gwas
  if (is.null(peak_to_locus_gwas) || all(is.na(peak_to_locus_gwas))) {
    peak_to_locus_gwas <- setNames(as.character(rowData(se_gwas)$signal), rownames(se_gwas))
  }

  # L2 weights = mean raw accessibility per peak (rowMeans of raw assay)
  # These weight each peak's "vote" so that more accessible peaks have more influence
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

  # --- Categorise loci ---
  # Same thresholds as rank concordance, but uses consensus_lineage (the lineage-level
  # consensus score) instead of dominant_score_locus for the "highly coherent" threshold.
  # This means a locus needs >50% of weighted votes going to one lineage to qualify.
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

# ============================================================
# GLOBAL ENRICHMENT
# ============================================================
# Tests whether GWAS-overlapping peaks have higher specificity
# RANKS for a given cell type than expected by chance (genome-wide).
#
# This is a rank-sum style test (related to Wilcoxon/Mann-Whitney):
#   1. Rank all peaks genome-wide by L2 specificity for each cell type
#   2. Take the mean rank of GWAS-overlapping peaks for each cell type
#   3. Compare to the expected mean rank under the null (uniform: (N+1)/2)
#   4. Z = (observed mean rank - expected) / SE, where SE = sqrt((N^2-1)/(12*n))
#   5. One-sided p-value (upper tail: are GWAS peaks ranked higher?)
#   6. BH FDR correction across cell types
#
# This is a GLOBAL (trait-level) metric, not locus-level.
# It asks: "across all GWAS loci combined, is this cell type enriched?"
# ============================================================
addGlobalEnrichment <- function(se, snps, mat_l2, results_dir, plots_dir, se_name, trait_name) {
  # Rank each cell type's L2 specificity column across all peaks (column-wise ranking)
  nmat_r <- apply(mat_l2, 2, function(x) rank(x, ties.method = "min"))
  nse_rank <- SummarizedExperiment(nmat_r, rowRanges = rowRanges(se))

  # Identify peaks that overlap GWAS SNPs
  gwas_rank_idx <- unique(subjectHits(findOverlaps(snps, rowRanges(nse_rank))))
  gwas_rank_mat <- if (length(gwas_rank_idx)) assay(nse_rank[gwas_rank_idx, , drop = FALSE], 1) else nmat_r[0, , drop = FALSE]

  if (!nrow(gwas_rank_mat)) return(data.frame())

  # Null expectation: if GWAS peaks were random, their mean rank = (N+1)/2
  # SE derived from the variance of the uniform distribution of ranks
  N <- nrow(nse_rank)   # total number of peaks genome-wide
  n <- nrow(gwas_rank_mat)  # number of GWAS-overlapping peaks
  exp_mean <- (N + 1) / 2
  se_z <- sqrt((N^2 - 1) / (12 * n))

  global_enrichment <- data.frame(observed_mean_rank = colMeans(gwas_rank_mat), stringsAsFactors = FALSE)
  global_enrichment$n_peaks <- n
  global_enrichment$exp_rank <- exp_mean
  global_enrichment$obs.exp <- global_enrichment$observed_mean_rank / exp_mean  # obs/exp ratio
  global_enrichment$z <- (global_enrichment$observed_mean_rank - exp_mean) / se_z
  global_enrichment$p <- pnorm(global_enrichment$z, lower.tail = FALSE)  # one-sided: higher ranks
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
