\name{CLEAR-specificity}
\alias{specPeakCore}
\alias{getSpecTiers}
\alias{getDominanceMask}
\alias{addSpecPeakTable}
\alias{addDominance}
\alias{addDominanceLineage}
\alias{addSpecificitySummary}
\title{Combined high-specificity and dominance annotations}
\description{
Annotate variant-overlapping peaks at cell-type and lineage resolution without
rounding scores before threshold decisions. The combined table is the default
public output; legacy dominance and specificity-summary helpers remain available.
}
\usage{
specPeakCore(mat)
getSpecTiers(spec_mat)
getDominanceMask(spec_mat, ratio_thresh = DEFAULT_RATIO_THRESH)
addSpecPeakTable(gwas_mat, snps, se_gwas, results_dir,
                cell_lineage = NULL, gtf_path = NULL,
                file = "spec_peak_summary.txt")
addDominance(mat, ratio_thresh = DEFAULT_RATIO_THRESH)
addDominanceLineage(mat, cell_lineage, ratio_thresh = DEFAULT_RATIO_THRESH)
addSpecificitySummary(mat_l2, cell_cor = NULL,
                      ratio_thresh = DEFAULT_RATIO_THRESH, threshold_sq = 0.8)
}
\arguments{
\item{mat, spec_mat, gwas_mat, mat_l2}{Named numeric peak-by-unit matrices.
Specificity rows must have unit L2 norm or be all zero.}
\item{ratio_thresh}{Top/runner-up ratio threshold; the default is 2.}
\item{snps}{A variant \code{GRanges}, with optional \code{names} and \code{signal}
metadata columns. Missing variant names fall back to genomic coordinates.}
\item{se_gwas}{The overlapping \code{SummarizedExperiment}, containing a
\code{norm} or \code{raw} accessibility assay. Rows are matched by peak name.}
\item{results_dir}{Directory for result tables, created if needed.}
\item{cell_lineage}{Named character vector mapping each cell type to a lineage.}
\item{gtf_path}{Optional local GTF path. A supplied but missing path is an error.
Protein-coding and lncRNA TSSs are used.}
\item{file}{Combined table filename.}
\item{cell_cor}{Optional correlation matrix for the legacy top-k coherence summary.}
\item{threshold_sq}{Cumulative squared-mass threshold for the legacy top-k summary.}
}
\details{
High specificity is strictly greater than \eqn{1/\sqrt{2}}. Low specificity is
strictly less than \eqn{1/\sqrt{\max(N,2)}}; moderate includes the boundaries.
Dominance is a separate rule: top/runner-up strictly greater than 2.
Comparisons allow only machine-roundoff tolerance at boundaries (eight machine
epsilons, scaled to the threshold), not rounding to displayed precision. Zero-signal
peaks and single-unit panels are not dominant. Positive peaks with a zero
runner-up have infinite ratios.

The combined table includes both resolutions, agreement between high and
dominant calls, accessibility, exact overlap-derived variant/signal identifiers,
and optional nearest-TSS context. Nearest genes do not identify regulatory targets.
}
\value{
\code{getSpecTiers} returns three logical matrices named \code{high},
\code{moderate}, and \code{low}. \code{getDominanceMask} returns one logical
matrix. Other functions return data frames; \code{addSpecPeakTable} also writes
a tab-separated table, including headers when there are no overlapping peaks.
}
