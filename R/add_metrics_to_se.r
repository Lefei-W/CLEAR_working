add_metrics_to_se <- function(se, input_assay = NULL, overwrite = FALSE, verbose = TRUE) {
  stopifnot(inherits(se, "SummarizedExperiment"))
  
  if (is.null(input_assay)) input_assay <- choose_input_assay(se)
  if (!input_assay %in% assayNames(se)) stop("Assay '", input_assay, "' not found in SE.")
  
  if (!overwrite && has_any_metrics(se)) {
    if (verbose) message("   … metrics already present; skipping add_metrics_to_se()")
    return(se)
  }
  
  mat <- SummarizedExperiment::assay(se, input_assay)
  if (ncol(mat) < 2) {
    warning("Need ≥2 columns to compute metrics; skipping.")
    return(se)
  }
  
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  stopifnot(!anyNA(mat), all(is.finite(mat)))
  
  rn <- rownames(se); if (is.null(rn)) rn <- seq_len(nrow(se))
  cn <- colnames(se); if (is.null(cn)) cn <- seq_len(ncol(se))
  
  ms <- compute_metrics(mat, rn, cn)
  for (nm in names(ms)) {
    a <- ms[[nm]]
    dimnames(a) <- list(rn, cn)
    SummarizedExperiment::assay(se, nm) <- a
  }
  
  if (verbose) message("   ✓ metrics added from assay: ", input_assay)
  se
}