test_that("high and dominance are independent, full-precision annotations", {
  x <- clear_fixture()
  se <- addSEMetrics(x$se)
  mat <- SummarizedExperiment::assay(se, "raw_l2")
  calls <- specPeakCore(mat)
  expect_true(calls$is_high[1])
  expect_false(calls$is_dominant_r20[1])
  expect_false(calls$is_high[2])
  expect_true(calls$is_dominant_r20[2])
  expect_equal(calls$ratio[4], 2)
  expect_false(calls$is_dominant_r20[4])
  expect_identical(calls$is_dominant_r20, unname(rowSums(getDominanceMask(mat)) > 0))
  expect_false(calls$is_dominant_r20[5])
  expect_false(calls$is_high[5])
  expect_true(is.na(calls$ratio[5]))
  expect_true(all(getSpecTiers(mat)$moderate[3, ]))
  expect_equal(rowSums(getSpecTiers(mat)$high), as.numeric(calls$is_high), ignore_attr = TRUE)
  expect_error(getSpecTiers(x$mat), "unit L2")
})

test_that("combined table contains both resolutions and exact signal memberships", {
  out <- withr::local_tempdir()
  args <- clear_report_fixture(out)
  peaks <- args$peak_cell
  expect_equal(nrow(peaks), 3L)
  expect_true(all(c("top_cell", "top_lineage", "is_high", "is_high_lineage",
                    "is_dominant_r20", "is_dominant_r20_lineage",
                    "call_agreement", "call_agreement_lineage",
                    "access_max", "nearest_gene") %in% names(peaks)))
  expect_identical(peaks$call_agreement, c("high_only", "dominant_only", "neither"))
  expect_true(all(peaks$is_high_lineage[peaks$is_high]))
  expect_equal(peaks$locus[2], "signal-1, signal-10")
  expect_equal(peaks$gwas_snps[2], "rs2, rs3")
  expect_equal(peaks$access_max, c(3, 3, 1))
  expect_true(file.exists(file.path(out, "spec_peak_summary.txt")))
  expect_false(file.exists(file.path(out, "specificity_summary_l2.txt")))
  expect_error(addSpecPeakTable(args$gwas_mat, args$snps, args$se_gwas, out,
                                gtf_path = file.path(out, "missing.gtf")), "GTF file")
})

test_that("zero and one overlapping peak retain table and matrix dimensions", {
  out <- withr::local_tempdir()
  args <- clear_report_fixture(out)
  for (n in c(0L, 1L)) {
    m <- args$gwas_mat[seq_len(n), , drop = FALSE]
    se <- args$se_gwas[seq_len(n), , drop = FALSE]
    lin <- getLineageMatrix(m, args$lineage_obj$cell_lineage)
    expect_equal(dim(lin), c(n, 3L))
    peaks <- addSpecPeakTable(m, args$snps, se, out,
                              cell_lineage = args$lineage_obj$cell_lineage)
    expect_equal(nrow(peaks), n)
    expect_true("call_agreement_lineage" %in% names(peaks))
    expect_true(file.exists(file.path(out, "spec_peak_summary.txt")))
  }
})

test_that("legacy summary functions remain usable from the package namespace", {
  out <- withr::local_tempdir()
  args <- clear_report_fixture(out)
  cell <- addSpecificity(args$gwas_mat, args$snps, 1 / sqrt(2), 1 / sqrt(6), out,
                         cell_lineage = args$lineage_obj$cell_lineage, se_gwas = args$se_gwas)
  lin <- addLineageSpecificity(args$gwas_mat, args$snps, out,
                               args$lineage_obj$cell_lineage, se_gwas = args$se_gwas)
  cmp <- getSpecificityComparison(cell, lin, out)
  expect_equal(nrow(cmp), nrow(args$gwas_mat))
  expect_identical(cell$is_dominant, args$peak_cell$is_dominant_r20)
  expect_equal(addDominanceLineage(args$gwas_mat, args$lineage_obj$cell_lineage)$ratio,
               args$peak_cell$ratio_lineage)
})

test_that("optional GTF annotation uses coding and lncRNA TSSs", {
  out <- withr::local_tempdir()
  args <- clear_report_fixture(out)
  gtf <- file.path(out, "genes.gtf")
  writeLines(c(
    'chr1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id "g1"; gene_name "CODING"; gene_type "protein_coding";',
    'chr1\ttest\tgene\t300\t400\t.\t+\t.\tgene_id "g2"; gene_name "LNCRNA"; gene_type "lncRNA";'
  ), gtf)
  peaks <- addSpecPeakTable(args$gwas_mat, args$snps, args$se_gwas, out,
                            cell_lineage = args$lineage_obj$cell_lineage, gtf_path = gtf)
  expect_equal(peaks$nearest_gene[1:2], c("CODING", "LNCRNA"))
  expect_equal(peaks$dist_to_tss[1:2], c(0, 0))
  expect_equal(peaks$tss_bin[1:2], rep("core_prom", 2))
})

test_that("pooling does not turn exact half-mass boundaries into high calls", {
  uniform <- matrix(rep(1 / sqrt(6), 6), 1,
                     dimnames = list("peak", LETTERS[1:6]))
  lin <- getLineageMatrix(uniform, setNames(rep(c("L1", "L2"), each = 3), LETTERS[1:6]))
  expect_false(any(getSpecTiers(lin)$high))
  expect_true(all(getSpecTiers(lin)$moderate))
  expect_false(specPeakCore(lin)$is_high)
  s <- 1 / sqrt(2) + 1e-12
  above <- matrix(c(s, sqrt(1 - s^2)), 1, dimnames = list("peak", c("A", "B")))
  expect_true(specPeakCore(above)$is_high)
  expect_true(any(getDominanceMask(matrix(c(2 + 1e-12, 1), 1))))
})
