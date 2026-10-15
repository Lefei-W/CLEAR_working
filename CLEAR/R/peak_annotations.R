.getTSS <- function(gtf) {
  fields <- names(mcols(gtf))
  if (!all(c("type", "gene_name") %in% fields)) {
    stop("GTF must contain type and gene_name annotations.")
  }
  type_col <- intersect(c("gene_type", "gene_biotype"), fields)[1]
  if (is.na(type_col)) stop("GTF must contain gene_type or gene_biotype.")
  genes <- gtf[gtf$type == "gene" &
                 mcols(gtf)[[type_col]] %in% c("lncRNA", "protein_coding")]
  GenomicRanges::promoters(genes, upstream = 0, downstream = 1)
}

.peakTssAnnotation <- function(peaks, gtf_path) {
  out <- data.frame(dist_to_tss = rep(NA_real_, length(peaks)),
                    tss_bin = rep(NA_character_, length(peaks)),
                    nearest_gene = rep(NA_character_, length(peaks)))
  if (is.null(gtf_path) || identical(gtf_path, "")) return(out)
  if (length(gtf_path) != 1L || is.na(gtf_path) || !file.exists(gtf_path)) {
    stop("GTF file does not exist: ", paste(gtf_path, collapse = ", "))
  }
  if (!length(peaks)) return(out)
  tss <- .getTSS(rtracklayer::import(gtf_path))
  if (!length(tss)) {
    warning("GTF contains no protein-coding or lncRNA genes; TSS annotations are unavailable.")
    return(out)
  }
  hits <- distanceToNearest(peaks, tss)
  pos <- queryHits(hits)
  out$dist_to_tss[pos] <- mcols(hits)$distance
  out$nearest_gene[pos] <- as.character(tss$gene_name[subjectHits(hits)])
  out$tss_bin <- as.character(cut(
    out$dist_to_tss, c(0, 1000, 5000, 50000, 200000, Inf),
    labels = c("core_prom", "prox_prom", "near_reg", "distal", "long_range"),
    include.lowest = TRUE))
  out
}

.addGWASPeakAnnotation <- function(summary_df, snps) {
  se <- SummarizedExperiment(rowRanges = GRanges(summary_df$peak))
  rownames(se) <- summary_df$peak
  map <- getSpecSignalVariantPeakMap(snps, se)
  collapse <- function(values) {
    by_peak <- split(values, map$peak)
    vapply(summary_df$peak, function(p)
      paste(unique(by_peak[[p]]), collapse = ", "), character(1), USE.NAMES = FALSE)
  }
  summary_df$gwas_snps <- collapse(map$variant)
  summary_df$gwas_locus <- collapse(map$signal)
  summary_df
}
