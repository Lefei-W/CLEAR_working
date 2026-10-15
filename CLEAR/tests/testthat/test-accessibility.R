test_that("default processing preserves original L2 scores and retains every peak", {
  x <- clear_fixture()
  se <- addSEMetrics(x$se)
  expected <- t(apply(x$mat, 1, function(v) if (sum(v^2) == 0) v else v / sqrt(sum(v^2))))
  expect_equal(SummarizedExperiment::assay(se, "raw"), x$mat)
  expect_equal(SummarizedExperiment::assay(se, "norm"), x$mat)
  expect_equal(SummarizedExperiment::assay(se, "raw_l2"), expected, tolerance = 1e-14)
  filtered <- filterAccessiblePeaks(flagAccessiblePeaks(se))
  expect_equal(nrow(filtered), nrow(se))
  expect_equal(SummarizedExperiment::assay(filtered, "raw_l2"), expected, tolerance = 1e-14)
  expect_identical(S4Vectors::metadata(filtered)$CLEAR$access_filter$mode, "none")
})

test_that("normalisation is opt-in, recorded, and validates cell counts", {
  x <- clear_fixture()
  counts <- seq(10, 60, 10)
  per_cell <- addSEMetrics(x$se, per_cell = "yes")
  expect_equal(SummarizedExperiment::assay(per_cell, "norm"), sweep(x$mat, 2, counts, "/"))
  scaled <- addSEMetrics(x$se, per_cell = "yes", depth = "colsum")
  expect_equal(unname(colSums(SummarizedExperiment::assay(scaled, "norm"))),
               rep(stats::median(colSums(sweep(x$mat, 2, counts, "/"))), 6))
  expect_identical(S4Vectors::metadata(scaled)$CLEAR$normalisation$depth, "colsum")
  expect_error(addSEMetrics(x$se, per_cell = "yes", n_cells = rep(0, 6)), "positive")
  expect_error(addSEMetrics(x$se, per_cell = "yes", n_cells = setNames(1:6, letters[1:6])), "cover")
  SummarizedExperiment::colData(x$se) <- S4Vectors::DataFrame(row.names = colnames(x$mat))
  expect_error(addSEMetrics(x$se, per_cell = "yes"), "Cell counts are required")
  expect_no_error(addSEMetrics(x$se))
})

test_that("accessibility filtering uses maximum signal and recomputes ranks", {
  x <- clear_fixture()
  se <- addSEMetrics(x$se)
  flagged <- flagAccessiblePeaks(se, min_max_quantile = .2)
  filtered <- filterAccessiblePeaks(flagged)
  expect_lt(nrow(filtered), nrow(se))
  expect_false(rownames(se)[5] %in% rownames(filtered))
  expected <- apply(SummarizedExperiment::assay(filtered, "raw_l2"), 2,
                     function(v) rank(-v, ties.method = "average"))
  expect_equal(SummarizedExperiment::assay(filtered, "raw_l2_rank"), expected)
  expect_error(flagAccessiblePeaks(se, min_max = 1, min_max_quantile = .1), "either")
  expect_error(flagAccessiblePeaks(se, min_max_quantile = 1), "\\[0, 1\\)")
  expect_error(filterAccessiblePeaks(flagAccessiblePeaks(se, min_max = 1000)), "every peak")
  zero_cut <- flagAccessiblePeaks(se, min_max = 0)
  expect_false(SummarizedExperiment::rowData(zero_cut)$pass_access[5])
})

test_that("filter losses are recorded in the input summary", {
  out <- withr::local_tempdir()
  x <- clear_fixture()
  all <- flagAccessiblePeaks(addSEMetrics(x$se), min_max = 2)
  se <- filterAccessiblePeaks(all)
  overlap <- addGWASOverlap(se, x$snps, se_all = all)
  expect_equal(overlap$n_dropped_by_filter, 1L)
  summary <- getInputSummary("tissue", "trait", "se.rds", "snps.rds",
                             se, x$snps, overlap, out)
  expect_equal(summary$n_gwas_peaks_dropped_by_filter, 1L)
  expect_equal(summary$access_filter_mode, "absolute")
  expect_equal(summary$access_filter_cut, 2)
})

test_that("quantile normalisation averages reference values over tied blocks", {
  mat <- matrix(c(1, 1, 3, 4, 8, 8), nrow = 3,
                 dimnames = list(paste0("p", 1:3), c("A", "B")))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  processed <- addSEMetrics(se, depth = "quantile")
  expect_equal(unname(SummarizedExperiment::assay(processed, "norm")),
               matrix(c(3.5, 3.5, 5.5, 2.5, 5, 5), nrow = 3))
  expect_equal(rowSums(SummarizedExperiment::assay(processed, "raw_l2")^2),
               rep(1, 3), ignore_attr = TRUE)
  expect_error(addSEMetrics(se, per_cell = "auto"), "Cell counts are required")
})
