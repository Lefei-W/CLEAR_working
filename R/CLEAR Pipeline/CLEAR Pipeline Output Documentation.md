# CLEAR Pipeline Output Documentation

## Analysis Structure

Pair-wise CLEAR analysis of snATAC-seq datasets across GWAS traits generates outputs in the following structure:

---

## Results Files

### 1. input_summary.txt

**Purpose:** Summary of input data dimensions and quality control metrics

**Columns:**
- `se_name` – snATAC dataset name
- `trait_name` – GWAS trait name  
- `GWAS_trait` – GWAS file basename
- `snATAC_dataset` – snATAC file basename
- `n_peaks` – Total peaks in snATAC dataset
- `n_celltypes` – Number of cell types
- `n_snps` – Total GWAS SNPs
- `n_unique_signals` – Number of unique GWAS loci/signals
- `n_gwas_peaks` – Peaks overlapping GWAS SNPs

---

### 2. specificity_summary_l2.txt

**Purpose:** Comprehensive per-peak specificity analysis combining basic and detailed metrics. Shows which cell types are accessible at each GWAS peak and whether specificity is concentrated or diffuse.

**Columns:**

*Peak identification:*
- `peak` – Peak coordinate (chr:start-end)
- `gwas_snps` – Comma-separated rsIDs overlapping this peak
- `gwas_locus` – Comma-separated signal/locus names

*Dominant cell type:*
- `maximum_cell` – Cell type with highest L2 value
- `maximum_score` – Highest L2 value
- `dominant_ratio` – Ratio of 1st/2nd highest L2 values
- `is_dominant` – TRUE if ratio > 2

*Specificity classification:*
- `n_high` – Number of cell types with L2 > 1/√3 (highly specific)
- `n_mid` – Number of cell types with mid-level specificity
- `n_low` – Number of cell types with low specificity
- `high_cells` – Comma-separated list of highly specific cell types
- `multiple_high` – TRUE if >1 highly specific cell type

*Detailed metrics:*
- `k_cells` – Number of top cell types explaining 50% of variance
- `top_k_cells` – Comma-separated list of those k cell types
- `coherence` – Mean correlation among top k cell types

---

### 3. specificity_summary_covweighted.txt

**Purpose:** Specificity analysis after accounting for cell-type correlations (e.g., T cells and NK cells often show co-accessibility)

**Structure:** Same columns as specificity_summary_l2.txt but computed using covariance-weighted matrix to decorrelate cell types

---

### 4. global_enrichment.txt

**Purpose:** Primary enrichment test. Tests if GWAS peaks rank higher than expected in each cell type's accessibility profile using the mean rank test.

**Columns:**
- `celltype` – Cell type name
- `observed_mean_rank` – Mean rank of GWAS peaks in this cell type
- `n_peaks` – Number of GWAS peaks
- `exp_rank` – Expected mean rank under null hypothesis: (N+1)/2
- `obs.exp` – Ratio of observed/expected (>1 = enriched)
- `z` – Z-score for rank enrichment test
- `p` – One-sided p-value (upper tail)
- `fdr` – Benjamini-Hochberg adjusted p-value
- `positive` – TRUE if FDR < 0.05 and Z > 0

**Interpretation:**
- Higher `z` values indicate stronger enrichment
- `fdr < 0.05` indicates statistically significant enrichment
- `obs.exp > 1` indicates GWAS peaks rank higher than random expectation

---

### 5. results/qn/specificity_summary_qn_l2.txt

**Purpose:** Sensitivity analysis - specificity after quantile normalization to make cell type distributions comparable

**Structure:** Same columns as specificity_summary_l2.txt

---

### 6. results/qn/global_enrichment_qn_l2.txt

**Purpose:** Mean rank test on quantile-normalized data (sensitivity analysis)

**Structure:** Same columns as global_enrichment.txt

---

## Plot Files

### QC and Correlation Plots

#### 1. {se}_{trait}_all_assay_cor_panel.pdf
- **Type:** Multi-panel correlation heatmap (12 panels)
- **Shows:** Cell-type correlations across all matrix transformations:
  - raw, quantile-normalized (qn)
  - L2-normalized (raw_l2, qn_l2)
  - Covariance-weighted (raw_covW, qn_covW)
  - Correlation-weighted (raw_corW, qn_corW)
  - Combined transformations (covW_l2, corW_l2, qn_covW_l2, qn_corW_l2)
- **Purpose:** Quality control to verify transformations preserve biological structure

#### 2. {se}_{trait}_GWAS_all_assay_cor_panel.pdf
- **Type:** Same as above, GWAS peaks only
- **Purpose:** Check if GWAS peaks have different correlation patterns than genome-wide peaks

---

### Specificity Plots (L2-normalized)

#### 3. {se}_{trait}_l2_specificity_bar.pdf
- **Type:** Stacked bar chart
- **Shows:** Proportion of GWAS peaks with 0, 1, 2, 3+ highly specific cell types
- **Purpose:** Overview of peak specificity distribution

#### 4. {se}_{trait}_stackbars_l2_0.5.pdf (multi-page)
- **Type:** Stacked bar plots (one page per cell type)
- **Shows:** For each cell type, displays GWAS peaks where that cell type has L2 ≥ 0.5
  - Height of each bar = L2 value for that cell type
  - Stack composition = contributions from all cell types
- **Purpose:** Distinguish true cell-type-specific peaks from shared accessibility
- **Threshold:** Moderate (L2 ≥ 0.5)

#### 5. {se}_{trait}_stackbars_l2_0.8.pdf (multi-page)
- **Type:** Same as above
- **Threshold:** Stringent (L2 ≥ 0.8)
- **Purpose:** Focus on highly specific peaks

---

### Dominance Plots

#### 6. {se}_{trait}_dominance_l2.pdf
- **Type:** Dodged bar plot
- **Shows:** Two bars per cell type:
  - Proportion of all peaks where cell type is dominant (ratio > 2)
  - Proportion of GWAS peaks where cell type is dominant
- **Purpose:** Identify which cell types are enriched for GWAS signals
- **Interpretation:** GWAS bar higher than genome-wide = enrichment

#### 7. {se}_{trait}_dominance_covweighted.pdf
- **Type:** Same as dominance_l2.pdf
- **Difference:** Uses covariance-weighted matrix (decorrelated cell types)
- **Purpose:** Dominance after accounting for cell-type correlations

#### 8. {se}_{trait}_gwas_covweighted_cor.pdf
- **Type:** Correlation heatmap (cell types × cell types)
- **Shows:** Correlations within GWAS peaks using cov-weighted values
- **Purpose:** Identify groups of co-accessible cell types at GWAS loci

---

### Summary Panels (Covariance-Weighted)

#### 9. {se}_{trait}_summary_panels_ABC.pdf
- **Type:** Three-panel figure
- **Panel A:** Pie chart of overall specificity distribution
  - High: L2 > 1/√3
  - Mid: 1/√n < L2 ≤ 1/√3
  - Low: L2 ≤ 1/√n
- **Panel B:** Histogram of number of highly specific cell types per peak
- **Panel C:** Boxplot of dominant ratio stratified by number of high-specific cells
- **Purpose:** Comprehensive overview of GWAS peak specificity patterns

#### 10. {se}_{trait}_summary_panel_D.pdf
- **Type:** Binary heatmap (peaks × cell types)
- **Shows:** Black = cell type highly specific for this peak, White = not
- **Purpose:** Visual identification of peak-celltype associations

#### 11. {se}_{trait}_summary_panel_E.pdf
- **Type:** Continuous heatmap (peaks × cell types)
- **Shows:** Full covariance-weighted matrix for GWAS peaks
- **Color scale:** Red = high accessibility, Blue = low accessibility
- **Purpose:** Detailed view of accessibility patterns

---

### Enrichment Plots

#### 12. {se}_{trait}_global_enrichment.pdf
- **Type:** Horizontal bar plot
- **Shows:** Z-scores for mean rank test, ordered by Z-score
- **Color:** Gradient based on Z-score (blue = depleted, red = enriched)
- **Purpose:** **Primary enrichment visualization** - which cell types show significant GWAS enrichment

#### 13. {se}_{trait}_global_enrichment_fdr.pdf
- **Type:** Horizontal bar plot
- **Shows:** -log10(FDR) values, ordered by significance
- **Reference line:** Dashed line at FDR = 0.05 (-log10(0.05) ≈ 1.3)
- **Purpose:** Emphasize statistical significance of enrichments

---

### Quantile-Normalized Plots (plots/qn/)

The following plots replicate the L2 and enrichment analyses using quantile-normalized data:

#### 14. {se}_{trait}_qn_l2_specificity_bar.pdf
- QN version of plot #3

#### 15. {se}_{trait}_stackbars_qn_l2_0.5.pdf
- QN version of plot #4

#### 16. {se}_{trait}_stackbars_qn_l2_0.8.pdf
- QN version of plot #5

#### 17. {se}_{trait}_dominance_qn_l2.pdf
- QN version of plot #6

#### 18. {se}_{trait}_global_enrichment_qn_l2.pdf
- QN version of plot #12

#### 19. {se}_{trait}_global_enrichment_qn_l2_fdr.pdf
- QN version of plot #13

**Purpose of QN plots:** Sensitivity analysis to verify results are robust to differences in cell type accessibility distributions

---

## Key Concepts

### Matrix Transformations

**L2 Normalization:**
- Each peak becomes a unit vector (sum of squares = 1)
- Makes peaks comparable regardless of total accessibility
- Emphasizes relative accessibility across cell types

**Quantile Normalization:**
- Makes cell type distributions identical
- Removes distributional differences between cell types
- Tests if results depend on these differences

**Covariance Weighting:**
- Decorrelates cell types based on genome-wide covariance structure
- Accounts for shared biology (e.g., T cells and NK cells)
- Identifies unique contributions of each cell type

### Specificity Metrics

**High/Mid/Low Classification:**
- High: L2 > 1/√3 ≈ 0.577 (highly specific)
- Mid: 1/√n < L2 ≤ 1/√3 (moderate specificity)
- Low: L2 ≤ 1/√n (low specificity, n = number of cell types)

**Dominance:**
- Dominant ratio > 2: Most accessible cell type is 2× higher than second-highest
- Indicates strong cell-type specificity

**Top-k Cells:**
- Smallest set of cell types explaining ≥50% of peak variance
- Low k + high coherence = concentrated specificity
- High k + low coherence = diffuse accessibility

### Enrichment Testing

**Mean Rank Test:**
- Ranks all peaks by accessibility in each cell type (descending)
- Computes mean rank of GWAS peaks
- Tests if GWAS peaks rank higher than expected by chance
- Z-score quantifies enrichment strength
- FDR corrects for multiple testing across cell types

---

## Typical Analysis Workflow

1. **Input QC:** Check `input_summary.txt` for reasonable overlap between GWAS and snATAC peaks

2. **Correlation QC:** Examine `*_all_assay_cor_panel.pdf` to verify cell-type relationships are preserved

3. **Global Enrichment:** Identify significant cell types in `global_enrichment.txt` and `*_global_enrichment.pdf`

4. **Specificity Analysis:** 
   - Use `specificity_summary_l2.txt` to identify which GWAS peaks are specific to enriched cell types
   - Visualize with `*_stackbars_l2_*.pdf` and `*_dominance_l2.pdf`

5. **Sensitivity Check:** Compare raw vs. QN results to assess robustness

6. **Signal-level Analysis:** Use `gwas_locus` column to aggregate results by GWAS signal/locus

---

## File Naming Convention

All files follow the pattern: `{se_name}_{trait_name}_{analysis_type}.{ext}`

Example: `breast_x_breastcancer_global_enrichment.pdf`

Where:
- `se_name` = snATAC dataset identifier
- `trait_name` = GWAS trait identifier  
- `analysis_type` = specific analysis or plot type
- `ext` = txt, pdf, or other file extension
