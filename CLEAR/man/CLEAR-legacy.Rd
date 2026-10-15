\name{CLEAR-legacy}
\alias{prepareSE}
\alias{loadCLEARInputs}
\alias{addGWASOverlap}
\alias{getInputSummary}
\alias{addLineageMap}
\alias{getLineageMatrix}
\alias{addSpecificity}
\alias{addLineageSpecificity}
\alias{getSpecificityComparison}
\alias{plotAssayCorHeatmaps}
\alias{plotAssayDensities}
\alias{plotLineageDendrogram}
\alias{plotCelltypeStackbars}
\alias{plotCelltypeStackbarsAll}
\alias{plotCelltypeStackbarsByLocus}
\alias{plotLineageStackbars}
\alias{plotLineageStackbarsAll}
\alias{plotLineageStackbarsByLocus}
\alias{plotSpecificityBar}
\alias{plotSpecificityDensity}
\alias{plotSpecificityWeighted}
\alias{plotWeightVsTss}
\alias{plotLineageCorHeatmap}
\alias{plotLineageDensities}
\alias{plotDominance}
\alias{plotDominanceLineage}
\alias{plotCoherence}
\alias{plotSignalComparison}
\alias{PLOT_ASSAYS}
\title{Input, lineage and compatibility APIs}
\description{
Existing CLEAR input, lineage, separate-table and plotting functions remain
available alongside the combined-table and offline-report APIs.
}
\details{
\code{prepareSE} standardises genomic coordinates and peak names;
\code{loadCLEARInputs} reads experiment and variant RDS files.
\code{addGWASOverlap(se, snps, se_all = NULL)} returns the overlapping experiment,
specificity matrix and overlap indices. Supplying the unfiltered \code{se_all}
also records the number of variant peaks removed by accessibility filtering.
\code{getInputSummary} writes run metadata and raw-signal summaries.

\code{addLineageMap} uses Pearson correlation distance, average-linkage
clustering and a user-chosen integer \code{k_lineage}. \code{getLineageMatrix}
pools squared cell-type specificity within each lineage and takes square roots,
including for zero or one peak.

\code{addSpecificity}, \code{addLineageSpecificity}, and
\code{getSpecificityComparison} retain their separate-table interfaces.
Use \code{addSpecPeakTable} for the combined output written by the default
runner. Legacy low/mid tiers keep their original equality convention.

The plot functions accept the corresponding legacy tables or specificity
matrices. \code{PLOT_ASSAYS} names the input and specificity assays used for QC.
Coherence and signal-comparison plots require externally supplied experimental
results; their presence does not enable permutation analysis in the public runner.
}
\seealso{\code{\link{addSEMetrics}}, \code{\link{addSpecPeakTable}},
\code{\link{writeSpecHtmlReport}}}
