clear_fixture <- function() {
  mat <- rbind(
    c(3, 2.5, .1, .1, .1, .1),
    c(3, 1.4, 1.4, 1.4, 1.4, 1.4),
    rep(1, 6),
    c(2, 1, .2, 0, 0, 0),
    rep(0, 6),
    c(.01, .02, .03, .01, .01, .01),
    c(0, 0, 10, 1, 0, 0),
    c(0, 1, 0, 8, 1, 0),
    c(0, 1, 0, 0, 6, 1),
    c(0, 1, 0, 0, 1, 9),
    c(2, 3, 1, 5, 4, 2),
    c(5, 2, 3, 1, 4, 2)
  )
  colnames(mat) <- LETTERS[1:6]
  starts <- seq(100, by = 200, length.out = nrow(mat))
  rownames(mat) <- paste0("chr1:", starts, "-", starts + 99)
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(starts, starts + 99))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat), rowRanges = gr,
    colData = S4Vectors::DataFrame(nCells = seq(10, 60, 10), row.names = colnames(mat)))
  snps <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(120, 320, 330, 520, 5000), width = 1))
  S4Vectors::mcols(snps)$names <- paste0("rs", 1:5)
  S4Vectors::mcols(snps)$signal <- c("signal-1", "signal-1", "signal-10",
                                    "signal-10", "signal-no-peak")
  list(se = se, snps = snps, mat = mat)
}

clear_report_fixture <- function(out) {
  x <- clear_fixture()
  se <- addSEMetrics(prepareSE(x$se))
  se <- flagAccessiblePeaks(se)
  mat <- SummarizedExperiment::assay(se, "raw_l2")
  lineage <- addLineageMap(mat, 3L)
  overlap <- addGWASOverlap(se, x$snps)
  peaks <- addSpecPeakTable(overlap$gwas_mat, x$snps, overlap$se_gwas, out,
                            cell_lineage = lineage$cell_lineage)
  list(se_name = "toy tissue", trait_name = "toy trait", mat_spec = mat,
       gwas_mat = overlap$gwas_mat, se_gwas = overlap$se_gwas, snps = x$snps,
       peak_cell = peaks, lineage_obj = lineage,
       input_summary = data.frame(per_cell_mode = "no", depth = "none",
                                   access_filter_mode = "none", access_filter_cut = NA_real_))
}
