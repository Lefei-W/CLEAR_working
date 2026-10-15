test_that("report joins shared peaks to exact signals and retains uncovered signals", {
  out <- withr::local_tempdir()
  args <- clear_report_fixture(out)
  data <- do.call(prepareSpecHtmlReportData, args)
  expect_equal(nrow(data$map), 4L)
  counts <- setNames(data$signal_summary$n_peaks, data$signal_summary$signal)
  expect_equal(unname(counts[c("signal-1", "signal-10", "signal-no-peak")]), c(2L, 2L, 0L))
  expect_equal(data$signal_summary$coverage[data$signal_summary$signal == "signal-no-peak"],
               "no_retained_peak")
  expect_equal(data$signal_summary$cell_category[data$signal_summary$signal == "signal-10"],
               "no_high")
  expect_equal(rowSums(data$cell_share), rep(1, 3), ignore_attr = TRUE)
  expect_equal(rowSums(data$lineage_share), rep(1, 3), ignore_attr = TRUE)

  args$trait_name <- "<script>alert('test')</script>"
  file <- file.path(out, "report", "index.html")
  do.call(writeSpecHtmlReport, c(args, list(
    output_file = file, mapping_file = file.path(out, "map.txt"),
    signal_file = file.path(out, "signals.txt"))))
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_gt(file.info(file)$size, 10000)
  expect_match(html, "Correlation and global enrichment", fixed = TRUE)
  expect_match(html, "data-signal-keys=", fixed = TRUE)
  expect_match(html, "High / dominant agreement", fixed = TRUE)
  expect_match(html, "&lt;script&gt;", fixed = TRUE)
  expect_false(grepl("<script>alert", html, fixed = TRUE))
  expect_false(grepl('src="https?://', html))
  expect_equal(nrow(read.delim(file.path(out, "map.txt"))), 4L)
  expect_equal(nrow(read.delim(file.path(out, "signals.txt"))), 3L)
})

test_that("report remains useful without overlapping peaks or enrichment", {
  out <- withr::local_tempdir()
  args <- clear_report_fixture(out)
  args$gwas_mat <- args$gwas_mat[0, , drop = FALSE]
  args$se_gwas <- args$se_gwas[0, , drop = FALSE]
  args$peak_cell <- addSpecPeakTable(args$gwas_mat, args$snps, args$se_gwas, out,
                                    cell_lineage = args$lineage_obj$cell_lineage)
  file <- file.path(out, "empty.html")
  do.call(writeSpecHtmlReport, c(args, list(output_file = file)))
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_match(html, "No retained variant peaks", fixed = TRUE)
  expect_match(html, "no enrichment result", fixed = TRUE)
  expect_match(html, "no_retained_peak", fixed = TRUE)
})

test_that("global enrichment is available without plotting or attached dplyr", {
  out <- withr::local_tempdir()
  args <- clear_report_fixture(out)
  ge <- addGlobalEnrichment(
    args$se_gwas, args$snps, args$gwas_mat, out, NULL, "tissue", "trait")
  expect_equal(nrow(ge), ncol(args$gwas_mat))
  expect_true(all(c("z", "fdr", "celltype") %in% names(ge)))
  expect_equal(ge$z, rep(0, nrow(ge)), tolerance = 1e-12)
  expect_error(addGlobalEnrichment(args$se_gwas, args$snps,
                                   args$gwas_mat[3:1, ], out, NULL, "tissue", "trait"),
                "rows must match")
  none <- addGlobalEnrichment(
    args$se_gwas, GenomicRanges::shift(args$snps, 10000L), args$gwas_mat,
    out, NULL, "tissue", "trait")
  expect_equal(nrow(none), 0L)
  expect_true(all(c("z", "fdr", "celltype") %in% names(none)))
  expect_equal(nrow(read.delim(file.path(out, "global_enrichment.txt"))), 0L)
})
