# CLEAR workflow and outputs (20-03 WIP: Step 7 & 8)

This document the CLEAR running STEPS and the details of the results tables and plots produced.

## 1) Main workflow (step-by-step)

The script runs one snATAC dataset x one GWAS trait pair.

### Inputs
- snATAC SummarizedExperiment RDS
- GWAS SNP/CCV GRanges RDS
- Permutation: GTF for TSS-distance matching (used in locus coherence analysis)

### Step 1. Load and prepare data
- Load GWAS SNP object.
- Load processed SE and process it.
- Processing includes:
  - Rebuilding row genomic ranges if needed. RSE
  - L2-normalizing each peak across cell types to create raw_l2.
  - Creating rank matrix raw_l2_rank.

### Step 2. Define core matrix (ALL on L2)
- Use raw_l2 as the main matrix for specificity, dominance, locus coherence, and global enrichment.

### Step 3. Find GWAS-overlapping peaks (snATAC matrix annotated with GWAS vector)
- Overlap SNP/CCV ranges with peak ranges.
- Create a GWAS-only SE subset (se_gwas).
- Aggregate annotation per overlapping peak:
  - CCVs
  - signal (locus label, OpenTargets lead SNP) rowData used for locus coherence. 

### Step 4. QC visualization
- Correlation heatmaps per assay, check lineages.
- Per-cell-type density plots for both the raw and L2-transformed matrices.
- Done for genome-wide SE and GWAS-overlap SE (when GWAS peak count > 1).

### Step 5. Input summary
- Save two summary files that capture data size and baseline signal structure.
- File: input_summary.txt (single-row run metadata):
  - se_name, trait_name, GWAS_trait, snATAC_dataset
  - n_peaks, n_celltypes, n_snps
  - n_unique_signals (all trait signals), n_gwas_peaks (overlapping peaks), n_unique_signals_in_peaks
  - mean_total_raw_signal, var_total_raw_signal across cell types
  - max_signal_celltype / max_signal_value and min_signal_celltype / min_signal_value
- File: raw_signal_per_celltype.txt (one row per cell type):
  - celltype
  - total_raw_signal (column sum of raw assay)
- Purpose:
  - Check assay scale and coverage before interpretation.
  - Confirm overlap size and locus diversity are reasonable for downstream tests.
  - Identify globally high- or low-signal cell types that may influence interpretation.

### Step 6. L2 specificity summary
- For each GWAS peak, compute:
  - max cell type and score
  - dominant ratio (top1/top2) Definition: >2 as dominant
  - number of high/mid/low cell types using fixed thresholds
  - top {k} cell set needed to explain 80% of L2 squared mass
  - coherence_l2 among those top-k cell types
- Threshold details used in Step 6:
  - Dominance threshold: dominant_ratio > 2 (is_dominant = TRUE)
  - High threshold: high_thresh = 1/sqrt(3) (about 0.577)
  - Mid threshold: mid_thresh = 1/sqrt(number_of_cell_types)
  - Category counting:
    - n_high counts cell types with L2 > high_thresh
    - n_mid counts cell types with mid_thresh < L2 <= high_thresh
    - n_low counts cell types with L2 <= mid_thresh
  - Top-k threshold: choose smallest k where cumulative L2^2 >= 0.8
- Coherence_l2 definition:
  - Matrix used: GWAS peak x cell-type raw_l2 matrix (not raw).
  - Build cell-type correlation matrix from that L2 matrix.
  - For each peak, take top_k_cells and compute mean pairwise correlation among those top cells.
- Interpretation intent:
  - low k_cells and high dominant_ratio suggest concentrated, cell-type-specific signal
  - higher k_cells and multiple_high suggest broader or shared programs
- Save specificity_summary_l2.txt.
- Plot summary bar and stackbar decomposition at cumulative 50% and 80% variance.

### TEST: Step 7. Dominance profile
## Run a test on the ratio (GWAS peaks dominant proportion over the genome-wide domiant) for each cell type NOTE they are dependent on each other, not independent tests.
- Compare dominant-cell-type composition between genome-wide peaks and GWAS peaks.
- Add per-cell-type dominance enrichment test using dominant peaks only:
  - Dominant peak definition: dominant_ratio > 2.


### Step 8. Locus-level permutation concordance
## WIP: 
## Currently it's taking properties from each signal, and grabbing peaks genome-wide, thinking of to do one permutation across the snATAC within a GWAS-ish region (500k window) and bin the permutations based on genes, peaks density ... IDEAS: 1. use 1000 random GWAS loci from different traits from OpenTargets 2. Check s-LDSC for baseline annotations 3. Might need to consider the rank correlation calculation, the top cell types or lineage (hypothesis!) might stay in similar ranks, but bottom noisy cell types could be crap and affect correlation, thinking of putting more weights on the top cell types like L2?  
- For each GWAS locus (signal), test whether GWAS peaks in that locus have higher concordant within-peak cell-type ranks.
- Build matched permutation null using joint bins:
  - TSS-distance bin
  - local peak-density bin
- How permutation works (per locus):
  - For each peak, convert cell-type profile to within-peak ranks. (within peaks (row-wise), global enrichment was done within cell types(column-wise))
  - For a locus with k peaks, compute observed concordance as mean pairwise correlation among the k rank vectors.
  - Create matched null by sampling k peaks from genome-wide pool with the same joint-bin composition as the locus.
  - Compute perm_mean, perm_sd, z_score, and p_value.
- Bins used for matching:
  - TSS-distance bins:
    - core_prom: 0-1 kb
    - prox_prom: 1-5 kb
    - near_reg: 5-50 kb
    - distal: 50-200 kb
    - long_range: >200 kb
  - Local density bins:
    - 4 quantile bins from number of neighboring peaks in a 50 kb window
  - Joint bin = interaction(TSS bin, density bin)
- Compute locus z-score and empirical p-value.
- Categorize loci by coherence and specificity rule set.
- Category rules used: # change the dominant naming, it was defined in the peak-level categorisation
  - Single_peak: n_peaks == 1 (no concordence)
  - Highly_coherent_specific: z_score > 2 AND p_value < 0.05 AND dominant_score_locus > 0.5
  - Coherent_moderate: z_score > 2 AND p_value < 0.05 AND dominant_score_locus <= 0.5
  - Incoherent: z_score < 1
  - Intermediate: all remaining loci
- Save locus_concordance_permutation.txt.

### Step 9. Global rank-based enrichment
## Thinking of putting the domiance test of GWAS peak over genome-wide here as well
- Rank all peaks per cell type genome-wide.
- Evaluate whether GWAS peaks have higher-than-expected mean rank per cell type.
- Compute z, p, FDR and save global_enrichment.txt.


## 2) Two key summary tables

### A) specificity_summary_l2.txt

One row per GWAS-overlapping peak.

| Column | Definition |
|---|---|
| locus | GWAS signal string |
| peak | Peak genomic ranges (chr:start-end). |
| gwas_snps | SNPs in the peak, comma-separated if multiple |
| maximum_cell | Cell type with the highest L2 value for this peak. |
| maximum_score | Highest L2 value for this peak. |
| dominant_ratio | Ratio of top L2 to second-highest L2 (top1/top2). |
| is_dominant | TRUE if dominant_ratio > 2. |
| n_high | Number of cell types with L2 > high_thresh. |
| n_mid | Number of cell types with mid_thresh < L2 <= high_thresh. |
| n_low | Number of cell types with L2 <= mid_thresh. |
| high_cells | Comma-separated names of cell types above high_thresh. |
| multiple_high | TRUE if more than one cell type is above high_thresh. |
| k_cells | Minimum number of top-ranked cell types needed so cumulative L2^2 >= 0.8. |
| top_k_cells | Comma-separated cell types included in that top-k set. |
| coherence_l2 | Mean pairwise correlation among top_k_cells in the cell-type correlation matrix computed from GWAS raw_l2 (not raw). |

**Key thresholds**

| Threshold | Value / rule |
|---|---|
| Dominance | dominant_ratio > 2 (DEFAULT_RATIO_THRESH = 2) |
| High-specificity | high_thresh = 1/sqrt(3) (~0.577) |
| Mid | mid_thresh = 1/sqrt(number_of_cell_types) |
| Top-k explanatory | cumulative L2^2 >= 0.8 for k_cells/top_k_cells |
### B) locus_concordance_permutation.txt

One row per GWAS locus (signal value).

| Column | Definition |
|---|---|
| locus | GWAS signal string |
| n_peaks | Number of GWAS-overlapping peaks in this locus. |
| concordance | Observed mean pairwise Spearman correlation of within-peak cell-type rank profiles across peaks in the locus. |
| perm_mean | Mean concordance from matched permutations. |
| perm_sd | SD of concordance from matched permutations. |
| z_score | (concordance - perm_mean) / perm_sd. |
| p_value | Empirical one-sided p-value, using +1 correction: (count(perm >= obs)+1)/(n_perm_eff+1). |
| dominant_cell_locus | Cell type with highest locus-average L2 signal. |
| dominant_score_locus | That maximum locus-average L2 value. |
| dominant_lineage | Dominant lineage cluster label from hierarchical clustering of cell types. |
| lineage_fraction | Fraction of locus-average signal captured by the dominant lineage. |
| significant | TRUE if p_value < 0.05. |
| category | Rule-based locus class. |

**Key thresholds and rules**

| Threshold / rule | Value |
|---|---|
| significant | p_value < 0.05 |
| Single_peak | n_peaks == 1 |
| Highly_coherent_specific | z_score > 2 AND p_value < 0.05 AND dominant_score_locus > 0.5 |
| Coherent_moderate | z_score > 2 AND p_value < 0.05 AND dominant_score_locus <= 0.5 |
| Incoherent | z_score < 1 |
| Intermediate | all remaining loci |

### global_enrichment.txt 

One row per cell type.

| Column | Definition |
|---|---|
| observed_mean_rank | Mean GWAS-peak rank for this cell type. |
| n_peaks | Number of GWAS peaks used in rank enrichment. |
| exp_rank | Null expected mean rank = (N+1)/2. |
| obs.exp | observed_mean_rank / exp_rank. |
| z | Z-score for rank enrichment. |
| p | One-sided p-value for rank enrichment. |
| fdr | BH-adjusted p-value for rank enrichment. |
| positive | TRUE if z > 0 and fdr < 0.05. |
| celltype | Cell type name. |


## 3) Plots and what each plot shows

## Correlation and density QC
- tissue_assay_cor_heatmap.pdf:
  - Cell-type x cell-type correlation heatmap within an assay (raw or raw_l2).
  - Shows similarity structure among cell types.
- tissue_assay_densities.pdf:
  - Density curves of values per cell type, faceted by assay.
  - Shows distributional shifts and dynamic range across assays/cell types.

## Specificity-focused plots
- se_trait_l2_specificity_bar.pdf:
  - Stacked bar of GWAS peaks by n_high (how many highly specific cell types per peak).
  - Labels show count and proportion.
- se_trait_stackbars_l2_cum50.pdf and se_trait_stackbars_l2_cum80.pdf:
  - Horizontal stacked bars per selected GWAS peak.
  - Segment height is L2^2 contribution per cell type (sum to 1).
  - Dashed line marks cumulative threshold (0.5 or 0.8).
  - Highlights whether few or many cell types explain each peak.

## Dominance plot
- se_trait_dominance_l2.pdf:
  - Side-by-side proportions of dominant peaks (ratio > 2) by maximum cell type.
  - Compares genome-wide set vs GWAS set.

## Locus coherence plot
- se_trait_locus_concordance_category.pdf:
  - Each point is a locus (n_peaks > 1).
  - x-axis: observed concordance.
  - y-axis: -log10(p-value) (capped at 99th percentile for display).
  - Color: locus category.
  - Shape: dominant lineage.
  - Size: number of peaks in locus.

## Global enrichment plots
- se_trait_global_enrichment.pdf:
  - Barplot of per-cell-type z-score for GWAS rank enrichment.
  - Positive significant cell types are highlighted.
- se_trait_global_enrichment_fdr.pdf:
  - Barplot of -log10(FDR) by cell type.
  - Dashed line at FDR = 0.05 threshold.


## Notes
- This pipeline intentionally uses raw + L2 + rank-based analyses and does not use quantile normalization or covariance-weighted matrices.
- n_perm controls permutation depth in locus analysis (default 100).