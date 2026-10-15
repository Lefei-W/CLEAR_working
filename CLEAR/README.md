# CLEAR

**C**ell-type-resolved **L**ocus-specific **E**nrichment **A**nalysis by **R**anking

> Quantifies the cell-type and lineage specificity of chromatin accessibility
> at GWAS variant loci, with combined peak annotations and an offline HTML
> explorer of cell-type, lineage, and signal-level results.

## Status

Early development (v0.90). 

The public core includes high specificity **and** dominance at cell and lineage
resolution. Signal categories are descriptive; locus-level permutation testing
is not part of this workflow.

## Installation

```r
install.packages(c("remotes", "BiocManager"),
                 repos = "https://cloud.r-project.org")

remotes::install_github(
  repo         = "Lefei-W/CLEAR_working",
  subdir       = "CLEAR",
  dependencies = TRUE,     # CRAN + Bioconductor (auto via biocViews:)
  upgrade      = "never"
)
```

`BiocManager` must be installed beforehand; `remotes` then resolves the
Bioconductor `Imports` automatically thanks to the `biocViews:` field
in `DESCRIPTION`.

## Quick start

See [vignettes/CLEAR_tutorial_step_by_step.R](vignettes/CLEAR_tutorial_step_by_step.R)
for a block-by-block walk-through, including the combined table and HTML report.

For building a GWAS `GRanges` input (incl. Open Targets and PLINK LD expansion examples), see
[vignettes/gwas_input_tutorial.md](vignettes/gwas_input_tutorial.md).

## Method overview

- **L2 specificity:** each non-zero peak has a unit-norm vector; its squared
  entries sum to one.
- **High specificity:** `s > 1/sqrt(2)`, meaning more than half of the squared
  L2 mass belongs to one cell type or lineage.
- **Dominance:** top score / runner-up score `> 2`. This is a different rule:
  a peak can be high-only, dominant-only, both, or neither.
- **Tiers:** high as above; low is `s < 1/sqrt(max(N, 2))`; moderate includes
  both boundaries. Uniform profiles are therefore moderate, not low.

Calls use full precision, with only machine-roundoff tolerance at boundaries
(eight machine epsilons, scaled to the threshold). Zero-signal peaks have no high or dominant call;
with only one unit, dominance is undefined. Lineage scores pool squared
cell-type scores and take the square root, preserving L2 geometry.

Global enrichment uses ascending specificity ranks, average ranks for ties,
and BH-adjusted one-sided p-values. `ties = "min"` reproduces the older
development rank convention. These peak-level statistics do not account for
within-locus dependence and are not a substitute for locus permutations.

## Combined table and offline report

The runner now writes **one** `results/spec_peak_summary.txt`, containing:

- top and runner-up cell/lineage, scores, squared shares, and dominance ratios;
- high and dominant flags and their agreement at both resolutions;
- cell-specific, lineage-specific, or shared resolution;
- high/moderate/low counts, accessibility, variant/signal IDs, and optional TSS
  distance and nearest-gene annotation.

`report/index.html` contains the correlation/enrichment landscape, searchable
signal and peak tables, linked contribution bars, and visible-row CSV export.
All CSS, SVG and JavaScript are embedded; no server, internet connection,
Shiny, Pandoc, DT or Plotly is needed to view it.

Signal navigation uses exact variant-to-peak overlaps, including peaks shared
by several signals. `results/signal_variant_peak_map.txt` records that join and
`results/signal_summary.txt` records the report's signal table. A signal with
no retained peak is not evaluated, rather than classified as `no_high`.
`single`, `coherent`, `complex`, and `no_high` refer to high-specificity peaks;
they are not significance or causal-probability calls. Nearest genes indicate
proximity, not regulatory targets.

Existing `addSpecificity()`, `addLineageSpecificity()`,
`getSpecificityComparison()` and plotting functions remain callable. Their
separate tables are no longer written by the default runner. Legacy table
functions retain their original low/mid boundary convention.
Legacy per-locus PDFs also remain callable, but the default signal view is now
the HTML report's exact overlap mapping rather than comma-grouped locus plots.

## Optional normalisation and filtering

**Defaults are unchanged scientifically:** `addSEMetrics()` does not divide by
cell counts or depth-match columns, and the runner retains every peak.
The original input is stored as `raw`, comparable accessibility as `norm`,
and its L2 scores/ranks as `raw_l2` / `raw_l2_rank`.

```r
se <- addSEMetrics(prepareSE(readRDS(se_path)))  # original scale, no filtering

# Optional: divide by cell counts, then match column totals.
# Counts come from colData(se)$nCells or a named vector / two-column TSV.
se <- addSEMetrics(prepareSE(readRDS(se_path)),
                   per_cell = "yes", depth = "colsum")

# Optional: remove the low-accessibility tail using max(norm), not row sums.
se_all <- flagAccessiblePeaks(se, min_max_quantile = 0.1)
se <- filterAccessiblePeaks(se_all)
overlap <- addGWASOverlap(se, snps, se_all = se_all)
```

`per_cell` accepts `no`, `yes`, or `auto`; auto uses integer-valued input as a
heuristic for summed counts, so use an explicit setting when provenance is
known. Missing/invalid cell counts are errors when division is requested.
`depth` accepts `none`, `colsum`, or `quantile`. These choices can change
specificity, lineage clustering, **and enrichment ranks**.

The filter also accepts an absolute `min_max` cutoff. Neither cutoff is enabled
by default. Ties can make a quantile remove a different fraction than requested;
a zero cutoff excludes zero-signal peaks. Dropped GWAS-peak counts and processing
settings are recorded in `results/input_summary.txt` and displayed in the report.

## Runner

From a separate output working directory:

```text
Rscript master_runner.R tissue snATAC.rds trait GWAS.rds [gtf_path] [n_perm] [k_lineage] [processed_dir]
```

Use the full path to the runner when it is not in the working directory.
GTF is optional; omit it in the R API or pass an empty fifth CLI argument.
The sixth CLI argument is retained for compatibility but no permutations run.

Environment variables:

| Variable | Default | Choices |
|---|---|---|
| `CLEAR_PER_CELL` | `no` | `no`, `yes`, `auto` |
| `CLEAR_DEPTH` | `none` | `none`, `colsum`, `quantile` |
| `CLEAR_NCELLS_TSV` | unset | Two-column cell-type/count TSV |
| `CLEAR_ACCESS_FILTER` | `none` | `none`, `q0.1`, or an absolute cutoff |
| `CLEAR_WRITE_HTML` | `1` | `1` or `0` |

Interactively, source the runner, edit its input variables, and call
`runCLEAR()`. The same options can be passed as arguments; `write_plots = FALSE`
omits PDFs but still produces tables and the HTML report.

Processed caches are separate from the input RDS and keyed by normalisation.
Source-file changes and changed cell-count inputs invalidate the cache.
Filtering is reapplied per run. Use separate output working directories for
different traits or sensitivity settings so results are not overwritten.

## Tests

```r
testthat::test_local("CLEAR")
```

## Licence

© 2026 Lefei Wang - Genome Variation and Regulation in Disease - QIMR Berghofer - Cancer Research Program 
