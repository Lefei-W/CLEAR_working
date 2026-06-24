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

## 23/06/2026 Discussion with Jonno:

The three are **nested**: variant ⊆ signal ⊆ window. The variant overlapping peaks and dominant peaks are the key focus of the analysis, signal and window are considered as background context. We are interested in the fraction of dominant peaks among the variant peaks, and how it differs compared to loci and window fractions.

This can be turned into a 2X3 contigency table (dominant / non-dominant × variant / signal / window) and test with Fisher's exact test, or Wilcoxon test, or proportion z-test. Noticed that the variant overlapping peaks are small numbers and the fraction of dominant peaks are concentrated at 0 0.5 and 1, so the statistical test may not be powerful. We can separate the loci into complex and simple loci, for complex loci with more variant peaks, the coherence metric is more meaningful, and for simple loci with limited variant peaks or dominant peaks (categorising by number of dominant peaks could be biased towards artificial coherence). 
- 151 loci with at least 1 variant peak, 100 loci with at least 1 dominant variant peak,
- 32 loci with more than 5 variant peaks, and 8 loci with more than 10 variant peaks, and the coherence metric is more meaningful for these complex loci.
- 19 loci with more than 3 dominant peaks.
- 16 loci with more than 5 variant peaks and more than 3 dominant peaks at dominance ratio > 2.
- 24 loci with more than 5 variant peaks and more than 3 dominant peaks at dominance ratio > 1.5.

Hypothesis: since we were able to identify global enrichment patterns in the mature luminal cells and mammary epithelial lineage for breast cancer, we expect to see coherence loci to those globally enriched cell types / lineages.

For complex loci we can annotate with coherence metric, and for simple loci we can annotate with at variant level, the top cell type and its dominance ratio. e.g. mammographic density with small sample size, only 33 variant peaks, but we can still annotate the top cell type and its dominance ratio, and identify significant global enrichment in adipocytes and fibroblasts in breast tissue. 

Jonno suggested that although we can use the dominance idea (ratio > 2) to identify globally enriched cell types and lineages. 168/461 variant peaks are dominant across 151 loci. The number of dominant peaks per locus is small. We not nessarily expect to see a null and observed distribution at locus level, but we can still annotate the top cell type lineage and its dominance ratio and the genomic background charateristics for all loci. 

From the locus_summary_breast table, we want the genomic background characteristics (n_snp, span, TSS density), and the number of dominant peaks to each cell type and lineage. Also the (dominant) peaks in variant peaks. 
We should also consider the global enrichment pattern, I would expect many of the loci with only one dominant peak is likely to be dominant in the globally enriched cell type, even those with multiple dominant peaks may still have all of them pointing to the same globally enriched cell type. **BUT** we are interested in annotating all loci, those with dominant peak(s) to the other cell types, they are the loci that post-GWAS studies have been struggling to interpret, and may be pointing to interesting novel biology and therapeutic targets.

To do this, we can use the locus_summary_breast table and the window_peaks_breast table to annotate for each loci, at each cell type and lineage.

## Variables to consider per locus
- locus size: n_snp, span, n_peaks
- TSS density
- number of dominant peaks (variant / signal / window)
- fraction of dominant peaks — the key metric
- **dominant cell type and lineage**

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


## Output files 

- 'locus_snp_level.csv' — one row per CCV, with TSS distance and bin annotation
- variant_peaks_<tissue>.csv — one row per peak overlapping a CCV, with dominance metrics and TSS annotation
- 'signal_' / 'window_' / 'variant_' peaks_<tissue>.csv — one row per peak in the respective subset, with dominance metrics and `is_variant_peak` flag, signal and window are background context of all peaks within the region, variant overlapping is the key focus
- 'locus_summary_<tissue>.csv' — *KEY summary table* one row per locus, with the key metric `frac_dom_cell_<r15 / r20>` and `frac_topcell_of_dom_r15 / r20` (the top cell type's share among dominant peaks) for the variant subset, and the same metrics for signal and window subsets as background context, plus other variables to plot against (locus size, TSS density, etc.)
- variant_peak_count_distribution.csv — loci by variant peak count and dominant peak count