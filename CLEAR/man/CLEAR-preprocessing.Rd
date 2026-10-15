\name{CLEAR-preprocessing}
\alias{addSEMetrics}
\alias{flagAccessiblePeaks}
\alias{filterAccessiblePeaks}
\title{Optional accessibility normalisation and filtering}
\description{
Compute L2 specificity while keeping per-cell normalisation, depth matching,
and accessibility filtering disabled by default.
}
\usage{
addSEMetrics(se, assay_name = NULL, keep_top_prop = NULL,
             per_cell = c("no", "auto", "yes"),
             depth = c("none", "colsum", "quantile"), n_cells = NULL)
flagAccessiblePeaks(se, min_max = NULL, min_max_quantile = NULL)
filterAccessiblePeaks(se)
}
\arguments{
\item{se}{A peak-by-cell-type \code{SummarizedExperiment}.}
\item{assay_name}{Source assay; defaults to the first assay.}
\item{keep_top_prop}{Optional legacy row-sum quantile filter. The newer
maximum-accessibility filter is separate and also disabled by default.}
\item{per_cell}{Whether to divide by cell numbers. \code{auto} uses integer-valued
input as a heuristic for summed counts; explicit provenance is preferable.}
\item{depth}{No depth matching, matching column sums to their positive median,
or rank-average quantile normalisation.}
\item{n_cells}{Positive cell counts, as a named numeric vector or a two-column
TSV path. Otherwise detected from cell-count columns in \code{colData(se)}.
Counts are required only when division is requested.}
\item{min_max}{Optional absolute minimum of a peak's maximum comparable accessibility.}
\item{min_max_quantile}{Optional quantile in [0,1) of maximum comparable accessibility.
Cannot be supplied together with \code{min_max}.}
}
\details{
\code{raw} retains the input assay; \code{norm} stores comparable accessibility;
\code{raw_l2} and \code{raw_l2_rank} store specificity and descending column ranks.
Defaults preserve the original L2 values. Normalisation choices can change
specificity, lineage clustering and enrichment ranks.

\code{flagAccessiblePeaks} records \code{pass_access} in row metadata without
subsetting. With neither cutoff supplied, every peak is retained. A requested
zero cutoff excludes zero-signal peaks. Quantile ties can change the fraction
retained. \code{filterAccessiblePeaks} applies the flags and recalculates ranks
on the retained universe; removing every peak is an error.

Settings and counts are stored in \code{metadata(se)$CLEAR}. Keep the flagged,
unfiltered object to pass as \code{se_all} to \code{addGWASOverlap}, so filter
losses can be reported.
}
\value{An updated \code{SummarizedExperiment}.}
