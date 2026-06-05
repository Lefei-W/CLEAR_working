# CLEAR

**C**ell-type-resolved and **L**ocus-specific **E**nrichment **A**nalysis by **R**anking

> Quantifies the cell-type and lineage specificity of chromatin accessibility
> at GWAS variant loci, and tests whether candidate causal variants at the same
> locus act on the same regulatory cellular context.

## Status

Early development (v0.90). 

## Installation

```r
install.packages(c("remotes", "BiocManager"),
                 repos = "https://cloud.r-project.org")

remotes::install_github(
  repo         = "Lefei-W/CLEAR_personal",
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
for a block-by-block walk-through (load SE -> pick `k_lineage` from the
dendrogram preview -> variant-level analysis + plots & summaries).

For building a GWAS `GRanges` input (incl. Open Targets and PLINK LD expansion examples), see
[vignettes/gwas_input_tutorial.md](vignettes/gwas_input_tutorial.md).

## Method overview

CLEAR scores ATAC peaks at GWAS loci on three axes:

- **L2 specificity** — per-peak L2-normalised cell-type vector.
- **Dominance** — fraction of peaks whose top cell type beats runner-up by a tunable ratio (default 2x).

- **Coherence** — *upcoming* for locus-level categorisation. For GWAS locus / signal / region with at least 2 variant-peaks, coherence measures whether they converge in one cell-type / lineage group.

All three are computed at both cell-type and lineage resolution; lineage
groups are derived from hierarchical clustering of the cell-type L2 matrix.

## Licence

© 2026 Lefei Wang - Genome Variation and Regulation in Disease - QIMR Berghofer - Cancer Research Program 
