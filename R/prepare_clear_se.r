# assays; rowRanges; rownames; peakid in rowData; genome info Ready for compute metrics; normalise the RSE structure for CLEAR
# Should cover most reasonable cases of downloaded SE objects
prepare_clear_se <- function(se, genome = "hg38") {
  if (!length(assayNames(se))) stop("No assay found in SE.")
  # Make sure to add the row ranges as ArchR
  rr <- rowRanges(se)
  if (is.null(rr) || length(rr) == 0L) { # Downloaded SE may not have rowRanges
    rd <- as.data.frame(rowData(se))
    if (all(c("seqnames", "start", "end") %in% names(rd))) {
      rr <- GRanges(seqnames = rd$seqnames, ranges = IRanges(start = rd$start, end = rd$end)) # Create GRanges from rowData GENERATE RSE
      rowRanges(se) <- rr
    } else {
      stop("No rowRanges and no usable seqnames/start/end in rowData.")
    }
  }
  
  if (is.null(rownames(se)) || anyNA(rownames(se)) || any(rownames(se) == "")) { # Check for rownames, if not create from rowRanges
    rr <- rowRanges(se)
    rownames(se) <- paste0(as.character(seqnames(rr)), ":", start(rr), "-", end(rr))
  }
  
  if (!"peakid" %in% names(rowData(se))) rowData(se)$peakid <- rownames(se)
  
  try(seqlevelsStyle(rowRanges(se)) <- "UCSC", silent = TRUE)
  try(GenomeInfoDb::genome(rowRanges(se)) <- genome, silent = TRUE)
  
  se
}

has_any_metrics <- function(se) any(DERIVED_ASSAYS %in% assayNames(se)) # Check if processed by compute metrics 

choose_input_assay <- function(se) { # In-house breast data processed into different assays, select some 
  prefs <- c("CPM","PeakMatrix_perbp","PeakMatrix")
  hit <- prefs[prefs %in% assayNames(se)]
  if (length(hit)) return(hit[1])
  assayNames(se)[1]
}