# CLEAR Summary Column Definitions

This document explains the two primary summary outputs from run_clear_pair_locus.R:

1. specificity_summary_l2.txt (per peak summary)
2. locus_concordance_permutation.txt (per signal/locus summary)

It includes:
- exact output columns
- how each column is computed
- key code snippets used in the current implementation

## Source Script

Main script:
- run_clear_pair_locus.R

Primary functions used:
- addSpecificity
- compute_specificity_summary
- add_gwas_peak_annotation
- addLocusCoherence
- compute_locus_rank_concordance_vec

----------------------------------------

## 1) specificity_summary_l2.txt

### What this table represents

One row per GWAS-overlapping peak (using raw_l2 values for that peak across cell types).

The final table is built in addSpecificity and written with:

~~~r
write.table(
  specificity_summary_l2,
  file.path(results_dir, "specificity_summary_l2.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
~~~

### Column definitions and calculations

| Column | Type | Definition | How calculated |
|---|---|---|---|
| peak | character | Peak identifier for row | rownames(gwas_mat) |
| gwas_snps | character | Comma-separated GWAS SNP names overlapping this peak | add_gwas_peak_annotation via overlap peak_gr vs snps |
| gwas_locus | character | Comma-separated unique signal labels overlapping this peak | add_gwas_peak_annotation from snps signal metadata |
| maximum_cell | character | Cell type with highest raw_l2 value for this peak | In compute_specificity_summary: sort decreasing and take first |
| maximum_score | numeric | Highest raw_l2 value for this peak | First value after decreasing sort |
| dominant_ratio | numeric | Top to second-best raw_l2 ratio | x_sorted[1] / x_sorted[2], NA if only one cell type |
| is_dominant | logical | Whether dominant_ratio > ratio_thresh | ratio_thresh default is 2 |
| n_high | integer | Number of cell types with high specificity | rowSums(gwas_mat > high_thresh), high_thresh = 1/sqrt(3) |
| n_mid | integer | Number of cell types with intermediate specificity | rowSums(gwas_mat > mid_thresh and <= high_thresh), mid_thresh = 1/sqrt(n_celltypes) |
| n_low | integer | Number of cell types with low specificity | rowSums(gwas_mat <= mid_thresh) |
| high_cells | character | Comma-separated cell types above high_thresh | apply row-wise: names(x)[x > high_thresh] |
| multiple_high | logical | Whether more than one cell type is high | lengths(high_cells) > 1 |
| k_cells | integer | Number of top-ranked cell types needed to reach cumulative L2^2 threshold | smallest k where cumsum(sorted(x)^2) >= threshold_sq (default 0.8) |
| top_k_cells | character | Comma-separated top k cell types from k_cells | top cells from descending raw_l2 |
| coherence | numeric | Mean pairwise correlation among top_k_cells based on cell correlation matrix | mean upper triangle of submatrix from cell_cor |

### Key code excerpts

#### Final selected columns

~~~r
specificity_summary_l2 <- specificity_summary_l2 %>%
  select(
    peak, gwas_snps, gwas_locus,
    maximum_cell, maximum_score, dominant_ratio, is_dominant,
    n_high, n_mid, n_low, high_cells, multiple_high,
    k_cells, top_k_cells, coherence
  )
~~~

#### High, mid, low thresholds and counts

~~~r
high_thresh <- 1 / sqrt(3)
mid_thresh <- 1 / sqrt(ncol(mat_l2))

specificity_basic <- data.frame(
  peak = rownames(gwas_mat),
  n_high = rowSums(gwas_mat > high_thresh),
  n_mid = rowSums(gwas_mat > mid_thresh & gwas_mat <= high_thresh),
  n_low = rowSums(gwas_mat <= mid_thresh)
)
~~~

#### Dominance and k_cells logic

~~~r
dominant_ratio <- if (length(x_sorted) > 1) x_sorted[1] / x_sorted[2] else NA
is_dominant = as.character(dominant_ratio > ratio_thresh)

cum_sq <- cumsum(x_sorted^2)
k <- which(cum_sq >= threshold_sq)[1]
if (is.na(k)) k <- length(x_sorted)
~~~

----------------------------------------

## 2) locus_concordance_permutation.txt

### What this table represents

One row per GWAS signal/locus (derived from rowData(se_gwas)$signal), summarizing whether peaks in that locus show stronger cell type rank concordance than matched random peaks.

The final table is built in addLocusCoherence and written with:

~~~r
write.table(
  locus_results,
  file.path(results_dir, "locus_concordance_permutation.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
~~~

### Column definitions and calculations

| Column | Type | Definition | How calculated |
|---|---|---|---|
| locus | character | Locus/signal label | split key from peak_to_locus |
| n_peaks | integer | Number of GWAS-overlapping peaks in this locus | length(peaks in locus_split[[locus]]) |
| concordance | numeric | Observed mean pairwise peak-to-peak rank concordance within locus | Spearman correlation among row ranks, then mean upper triangle |
| perm_mean | numeric | Mean null concordance from matched permutations | mean(perm_vals, na.rm = TRUE) |
| perm_sd | numeric | SD of null concordance from matched permutations | sd(perm_vals, na.rm = TRUE) |
| z_score | numeric | Standardized observed concordance | (concordance - perm_mean) / perm_sd if perm_sd > 0 else NA |
| p_value | numeric | Empirical one-sided p-value for concordance enrichment | (count(perm_vals >= concordance) + 1) / (n_perm_eff + 1) |
| dominant_cell_locus | character | Cell type with highest mean raw_l2 across locus peaks | which.max(colMeans(gwas_mat[peaks, ])) |
| dominant_score_locus | numeric | Highest mean raw_l2 among cell types for locus | max(colMeans(gwas_mat[peaks, ])) |
| dominant_lineage | character | Dominant lineage among cell types for locus | sum mean_profile by lineage and take max |
| lineage_fraction | numeric | Fraction of locus mean signal in dominant lineage | max(lineage_profile) / sum(lineage_profile) |
| significant | logical | Nominal significance flag | p_value < 0.05 |
| category | character | Heuristic locus category | case_when using n_peaks, z_score, p_value, dominant_score_locus |

### Matching strategy used for permutations

Matched null peaks are sampled from full_mat using joint bins:
- TSS distance bin: core_prom, prox_prom, near_reg, distal, long_range
- Peak density bin: quantile bins from local peak density
- Joint bin: interaction(tss_bin, density_bin)

For each locus, the same number of peaks per joint bin is sampled from genome-wide pools for each permutation.

### Key code excerpts

#### Observed concordance

~~~r
rank_gwas <- t(apply(gwas_mat, 1, rank))
obs_mat <- rank_gwas[peaks, , drop = FALSE]
obs_cor <- cor(t(obs_mat), method = "spearman")
obs_conc <- mean(obs_cor[upper.tri(obs_cor)])
~~~

#### Matched permutation null

~~~r
locus_joint_bins <- joint_bin_full[peaks]
locus_bins <- table(locus_joint_bins)
full_bins <- split(names(joint_bin_full), joint_bin_full)

draws <- replicate(
  n_perm,
  sample(pool, n_bin, replace = length(pool) < n_bin),
  simplify = "matrix"
)
~~~

#### Empirical p-value and z-score

~~~r
perm_mean <- mean(perm_vals, na.rm = TRUE)
perm_sd <- sd(perm_vals, na.rm = TRUE)
n_perm_eff <- sum(!is.na(perm_vals))
z_score <- ifelse(is.na(perm_sd) || perm_sd == 0, NA_real_, (obs_conc - perm_mean) / perm_sd)
p_value <- if (n_perm_eff > 0) {
  (sum(perm_vals >= obs_conc, na.rm = TRUE) + 1) / (n_perm_eff + 1)
} else {
  NA_real_
}
~~~

#### Category assignment

~~~r
locus_results$category <- dplyr::case_when(
  locus_results$n_peaks == 1 ~ "Single_peak",
  locus_results$z_score > 2 & locus_results$p_value < 0.05 & locus_results$dominant_score_locus > 0.5 ~ "Highly_coherent_specific",
  locus_results$z_score > 2 & locus_results$p_value < 0.05 ~ "Coherent_moderate",
  locus_results$z_score < 1 ~ "Incoherent",
  TRUE ~ "Intermediate"
)
~~~

----------------------------------------

## Notes on NA and edge behavior

- specificity_summary_l2.txt:
  - dominant_ratio can be NA if only one cell type exists.
  - coherence can be NA if top_k_cells has only one cell type.

- locus_concordance_permutation.txt:
  - n_peaks < 2 returns NA for concordance/null stats.
  - If required joint bins are unavailable in full_bins, permutation stats are returned as NA.
  - If permutation variance is zero, z_score is NA.
  - p_value uses +1 correction, so exact 0 is avoided when permutations are valid.

----------------------------------------

## Output file names in this pipeline

- specificity summary:
  - results/specificity_summary_l2.txt

- locus coherence summary:
  - results/locus_concordance_permutation.txt
