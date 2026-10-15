\name{CLEAR-report}
\alias{getSpecSignalVariantPeakMap}
\alias{prepareSpecHtmlReportData}
\alias{writeSpecHtmlReport}
\title{Standalone offline specificity report}
\description{
Build an offline HTML explorer with cell and lineage landscapes, exact signal
navigation, peak contribution bars, sortable/searchable tables and CSV export.
}
\usage{
getSpecSignalVariantPeakMap(snps, se_gwas)
prepareSpecHtmlReportData(
  se_name, trait_name, mat_spec, gwas_mat, se_gwas, snps,
  peak_cell, locus_categorisation = data.frame(),
  global_enrichment_cell = data.frame(),
  global_enrichment_lineage = data.frame(),
  lineage_obj, input_summary = data.frame())
writeSpecHtmlReport(
  se_name, trait_name, mat_spec, gwas_mat, se_gwas, snps,
  peak_cell, locus_categorisation = data.frame(),
  global_enrichment_cell = data.frame(),
  global_enrichment_lineage = data.frame(),
  lineage_obj, input_summary = data.frame(), output_file,
  mapping_file = NULL, signal_file = NULL)
}
\arguments{
\item{se_name, trait_name}{Dataset and trait display names.}
\item{mat_spec}{Genome-wide specificity matrix on the retained peak universe,
usually the \code{raw_l2} assay.}
\item{gwas_mat}{Specificity matrix for variant-overlapping peaks.}
\item{se_gwas}{Overlapping experiment, with peak names matching \code{gwas_mat}.}
\item{snps}{Variant \code{GRanges}; \code{signal} metadata enables signal navigation.}
\item{peak_cell}{Combined table returned by \code{addSpecPeakTable}.}
\item{locus_categorisation}{Optional precomputed high-specificity categories.
By default categories are derived from the exact signal-to-peak map.}
\item{global_enrichment_cell, global_enrichment_lineage}{Optional tables returned
by \code{addGlobalEnrichment}. Empty tables are explicitly displayed as unavailable.}
\item{lineage_obj}{Lineage map and palettes returned by \code{addLineageMap}.}
\item{input_summary}{Optional run metadata from \code{getInputSummary}.}
\item{output_file}{Required HTML output path.}
\item{mapping_file}{Optional path for the exact signal/variant/peak overlap TSV.}
\item{signal_file}{Optional path for the report's per-signal summary TSV.}
}
\details{
All CSS, SVG and JavaScript are embedded. No web server or network connection is
needed to view the report. Shared peaks are joined to every overlapping signal
using exact genomic overlaps rather than comma-string matching.

Signals without retained peaks remain visible as not evaluated. This differs
from \code{no_high}, where retained peaks exist but none is high. Signal
categories are descriptive, not significance tests or causal probabilities.
}
\value{
The mapping function returns a data frame. Data preparation returns a list of
matrices, tables, palettes and metadata. The writer returns the HTML path.
}
