#### 2. RSE processing ####

# Assign the peak IDs 
ensure_peakids <- function(se) {
  if (is.null(rownames(se)) || anyNA(rownames(se))) {
    rr <- rowRanges(se) # should all have rowRanges as a RSE object from ArchR normalised 
    rownames(se) <- paste0(as.character(seqnames(rr)), ":", start(rr), "-", end(rr))
  }
  if (!"peakid" %in% names(rowData(se))) rowData(se)$peakid <- rownames(se) # assign peakid as rowData
  se
}

