# Locus-level coherence summary

## The idea

For a GWAS trait we have fine-mapped **credible-set SNPs** grouped into **loci**
(a group of SNPs in LD with a index variant grouped into a signal). For the matched tissue we have a
single-cell ATAC peak matrix where every peak carries a **specificity profile**
across cell types (the `raw_l2` assay): a unit-norm vector measuring cell type contribution to the peak's overall accessibility.

We define a peak is **dominant** when one cell type ranked first (specificity) carries far more of its
specificity than the second contributor, i.e. the top/2nd ratio by a threshold (1.5 or 2). (The L2-based threshold of 1/sqrt(2) ≈ 0.707 corresponds to the top cell contributing > 50% of the specificity could also be used. And the same idea applies to lineages)

The question: **do GWAS loci concentrate regulatory activity into a single cell
type, and what predicts that?** We score each locus two ways:

- **variant-peak subset** — only peaks that physically overlap a credible SNP.
- **GWAS-window subset** — *all* peaks inside the locus window (SNP span ± 100 kb),
  the local genomic characteristics for context.

For each subset we count dominant peaks and, compute the **fraction of
peaks that are dominant**. That fraction is the main metric (below); the rest
of the columns are variables (*potential confounders*) to plot against it.

Currently controlling for locus size (n_snp, span, peak counts, TSS density). Coherence could be driven by other potential confounders. A locus can look 'coherent' with few dominant peaks, we should condition on number of variant peaks and dominant peaks (less stringent dominance threshold). 
We can also construct a tissue-level baseline (permutation). 

## Column dictionary

### `locus_summary_<tissue>.csv` — one row per gwas × locus
Geometry / context:
| column | meaning |
|---|---|
| `locus`, `gwas`, `tissue` | identifiers |
| `matched` | is this the tissue matched to the trait (BCAC↔breast, HF↔heart) |
| `pass_filter` | `variant_n_peaks > 1` (i.e. ≥ 2 SNP-overlapping peaks) — locus has enough variant peaks |
| `seqnames` | chromosome |
| `n_snp` | credible SNPs in the locus |
| `snp_span_start/end`, `snp_span_width_bp` | extent of the credible set = the **signal-wise span** (first SNP → last SNP) |
| `window_start/end`, `window_width_bp` | the ± 100 kb window |
| `n_tss_in_window` | TSS count in the window |
| `gene_density_per_kb` | `n_tss_in_window / (window_width_bp / 1000)` |

Dominance summary, prefixed `variant_` (variant / SNP-overlapping subset), `signal_`
(signal / SNP-span subset) and `window_` (± 100 kb window subset). The three subsets are
nested (variant ⊆ signal ⊆ window), so e.g. `variant_n_peaks ≤ signal_n_peaks ≤ window_n_peaks`:
| column (drop prefix) | meaning |
|---|---|
| `n_peaks` | peaks in the subset |
| `n_dom_cell_r15`, `n_dom_cell_r20` | dominant peaks at cell level (ratio > 1.5 / 2) |
| `n_dom_lin_r15`, `n_dom_lin_r20` | dominant peaks at lineage level |
| **`frac_dom_cell_r15`, `frac_dom_cell_r20`** | **fraction of peaks that are cell-dominant — the key metric** |
| **`frac_dom_lin_r15`, `frac_dom_lin_r20`** | **fraction lineage-dominant** |
| `top_dom_cell_r15`, `top_dom_cell_r20` | most frequent cell type among the dominant peaks (ratio > 1.5 / 2) — the cell type that "wins" the locus |
| **`frac_topcell_of_dom_r15`, `frac_topcell_of_dom_r20`** | **its share: (dominant peaks of that top cell type) / (all dominant peaks). Denominator = dominant peaks, not all peaks. Locus coherence** |
| `top_dom_lin_r15`, `top_dom_lin_r20` | same at lineage level |
| **`frac_toplin_of_dom_r15`, `frac_toplin_of_dom_r20`** | **lineage share among dominant peaks** |
| `dom_cell_list` | `;`-joined dominant cell labels, strongest first (repeats kept) |
| `dom_lin_list` | same at lineage level |

### `signal_peaks_<tissue>.csv` — one row per peak inside the SNP span
All peaks overlapping `[snp_span_start, snp_span_end]` (first SNP → last SNP, **no**
± 100 kb expansion). Same metric columns as the variant table, plus `is_variant_peak`
(TRUE if the peak also overlaps a credible SNP). No `gwas_snps` / `n_snps`.

### `window_peaks_<tissue>.csv` — one row per peak in the ± 100 kb window
Same metric columns as the variant table, plus `is_variant_peak` (TRUE if the peak
also overlaps a credible SNP). No `gwas_snps` / `n_snps`.

### `locus_snp_level.csv` — one row per credible SNP
| column | meaning |
|---|---|
| `locus` | association signal id (index variant for OpenTargets traits) |
| `gwas` | trait name (BCAC / HF) |
| `snp` | SNP id |
| `seqnames`, `start`, `end` | SNP position |
| `dist_to_tss_bp` | to nearest protein-coding/lncRNA TSS |
| `nearest_tss_gene` | gene's name |
| `tss_bin` | distance bin: core_prom / prox_prom / near_reg / distal / long_range `c(0, 1000, 5000, 50000, 200000, Inf)` |

### `variant_peaks_<tissue>.csv` — one row per peak *overlapping* a credible SNP
| column | meaning |
|---|---|
| `peak` | peak id |
| `gwas_snps` | `;`-joined credible SNPs this peak overlaps |
| `n_snps` | how many credible SNPs overlap the peak |
| `maximum_cell` | top cell type in the peak's L2 profile |
| `maximum_score` | that cell type's L2 |
| `dominant_ratio` | top / 2nd cell-type L2 (the dominance strength) |
| `is_dominant_r15`, `is_dominant_r20` | dominant_ratio > 1.5 / > 2 |
| `top_lineage` | top lineage (lineage-summed L2) |
| `top_lineage_score` | that lineage's L2 |
| `top_lineage_ratio` | top / 2nd lineage ratio |
| `top_lineage_is_dominant_r15/r20` | top_lineage_ratio > 1.5 / > 2 |
| `dist_to_tss`, `nearest_tss_gene`, `tss_bin` | nearest-TSS annotation of the peak |



