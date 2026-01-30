index_se_library <- function(root = SE_LIBRARY_ROOT, genome = c("hg19","hg38")) {
  genome <- match.arg(genome)
  pat <- paste0("\\.", genome, "\\.(se\\.)?rds$")
  paths <- list.files(root, pattern = pat, recursive = TRUE, full.names = TRUE)
  if (!length(paths)) return(tibble::tibble())
  
  root_norm <- normalizePath(root, mustWork = FALSE)
  mk_rel <- function(p) sub(paste0("^", root_norm, "(/|\\\\)?"), "", normalizePath(p, mustWork = FALSE))
  
  meta <- lapply(paths, function(p) {
    rel_dir <- mk_rel(dirname(p))
    elems <- unlist(strsplit(rel_dir, "[/\\\\]", perl = TRUE))
    n <- length(elems)
    last_is_bin <- n >= 3 && grepl("^bin-", elems[n]) # nomencalture: bin-<size>
    if (last_is_bin) {
      bin   <- elems[n]
      assay <- elems[n-1]
      input <- elems[n-2]
    } else {
      bin   <- NA_character_
      assay <- elems[n]
      input <- elems[n-1]
    }
    data.frame(input = input, assay = assay, bin = bin,
               genome = genome, path = p, stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(meta)
}