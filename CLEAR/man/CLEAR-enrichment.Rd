\name{CLEAR-enrichment}
\alias{addGlobalEnrichment}
\title{Global rank enrichment of variant-overlapping peaks}
\description{Compare mean specificity ranks of variant peaks to a genome-wide
uniform-rank expectation at cell or lineage resolution.}
\usage{
addGlobalEnrichment(se, snps, mat_l2, results_dir, plots_dir,
                    se_name, trait_name, level = c("cell", "lineage"),
                    cell_lineage = NULL, lineage_palette = NULL,
                    lineage_colours = NULL, ties = c("average", "min"),
                    write_z_plot = FALSE)
}
\arguments{
\item{se}{Experiment defining the retained peak universe.}
\item{snps}{Variant \code{GRanges}.}
\item{mat_l2}{Specificity matrix in the same row order as \code{se}.}
\item{results_dir}{Table output directory.}
\item{plots_dir}{PDF directory, or \code{NULL} to omit plots.}
\item{se_name, trait_name}{Display names.}
\item{level}{Cell-type or lineage resolution.}
\item{cell_lineage}{Named cell-to-lineage mapping, required for lineage resolution.}
\item{lineage_palette, lineage_colours}{Named cell-type or lineage colours.}
\item{ties}{Rank tie convention. Average ranks preserve the null mean;
\code{min} reproduces the older development convention.}
\item{write_z_plot}{Also write the signed z-score PDF. The HTML report always
displays signed z-scores when the enrichment table is supplied.}
}
\details{
The table contains z-scores, one-sided p-values and BH-adjusted FDR.
Correlated peaks within loci violate independent-sampling assumptions; these
statistics do not replace locus-level permutations.
}
\value{An enrichment data frame, or an empty data frame when no GWAS peaks overlap.}
