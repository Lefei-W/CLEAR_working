# This is the full script from HPC
# Would need to modularize it for use in package
#### Packages and Setup #######

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(SummarizedExperiment); library(DelayedArray)
  library(ggsci); library(GenomicRanges); library(S4Vectors); library(readr); library(Matrix)
  library(preprocessCore); library(rtracklayer); library(GenomeInfoDb); library(tools)
})

source('~/plottheme.R')   # must define theme_Publication()

wd <- '/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/'
setwd(wd)

SE_LIBRARY_ROOT <- "output/snATAC_SE_library" # nolint

input_files <- list(
  union_peak = list(
    path   = "data/PeakCalls_Xiong/union_GRSE.rds",
    assays = c("PeakMatrix", "union_bw_signal_perbp")
  ),
  breast_final_union_331832 = list(
    path   = '/working/lab_jonathb/lefeiW/projects/CLEAR_breast_demo/results/union_peak331832_CLEAR_se.rds',
    assays = c("PeakMatrix")
  ),
  breast_final_union_531279 = list(
    path   = '/working/lab_jonathb/lefeiW/projects/CLEAR_breast_demo/results/union_peak531279_CLEAR_se.rds',
    assays = c("PeakMatrix")
  ),
  merged = list(
    path   = "data/merged_peakMatrix_perBP_normalised.rds",
    assays = c("PeakMatrix_perbp", "merged_bw_signal")
  ),
  pbmc10k = list(
    path   = "/working/joint_projects/bc_risk_locus_multiomics/published_scATAC/processed/pbmc10k_CLEAR_ready.se.rds",
    assays = c("PeakMatrix_perbp")
  ),
  hten_dcis = list(
    path   = "/working/joint_projects/bc_risk_locus_multiomics/published_scATAC/processed/hten_dcis.rds",
    assays = c("PeakMatrix_perbp")
  ),
  zhang2021_mammary = list(
    path   = "output/Zhang2021_mammary.hg38.rds",
    assays = c("CPM")
  ),
  zhang2021_mammary_ArchR = list(
    path   = "/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/output/snATAC_SE_library/zhang2021_mammary/zhang2021_mammary_ArchR_se.rds",
    assays = c("PeakMatrix")
  ),
  zhang2021_mammary_CPM_peak_filtered = list(
    path   = "/working/joint_projects/bc_risk_locus_multiomics/published_scATAC/downloaded/scATAC_corpus/Zhang2021/Zhang2021-mammary_tissue/output/Zhang2021_mammary_CPM_peak_filtered.RDS",
    assays = c("CPM")
  ),
  zhang2021_pancreras_ArchR = list(
    path   = "/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/output/snATAC_SE_library/zhang2021_pancreas/zhang2021_pancreas_ArchR_se.rds",
    assays = c("PeakMatrix")
  ),
  regner2025_mammary_ArchR = list(
    path   = "/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/output/snATAC_SE_library/regner2025_mammary_ArchR/regner_2025_mammary_ArchR.rds",
    assays = c("PeakMatrix")
  ),
  kanemaru2023_heart_ArchR = list(
    path   = "/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/output/snATAC_SE_library/kanemaru2023_heart_ArchR/Kanemaru_2023_heart_ArchR_se.rds",
    assays = c("PeakMatrix")
  )
)

DROP_CTS_DEFAULT <- c("Basal2", "Intermediate") # noisy in-house clusters

CHAIN_FILE <- "/working/lab_jonathb/lefeiW/projects/sc_eQTL/results/Cleaned_BC_subtype_GWAS/hg38ToHg19.over.chain"



#### 1. CLEAR core specificity ####

safe_l2 <- function(x) { # on the rows (peaks)
  s <- sum(x^2)
  if (s == 0) x else x / sqrt(s)  # l2 euclidean normalization avoiding row sum of 0 (happen with the Human-scATAC-Corpus)
}

add_rank <- function(M, rn, cn) { # ranking assign to cell types (columns)
  R <- apply(M, 2, function(x) rank(x, ties.method = "min")) # rank 1 is the smallest, and tied values as the minimal rank
  dimnames(R) <- list(rn, cn)
  R # ranking matrix 
}

## compute all single-matrix transformations + row-level summaries
## Returns a named list of matrices (same dims as input)
.transform_single_mat <- function(mat, mat_raw, rn, cn) {
  k <- ncol(mat)

  # --- per-row transformations (each returns peak × celltype matrix) ---
  zscaled   <- t(apply(mat, 1, scale));       dimnames(zscaled)  <- list(rn, cn)
  l2_norm   <- t(apply(mat, 1, safe_l2));     dimnames(l2_norm)  <- list(rn, cn)

  # Proportion: signal(peak,ct) / row_sum  (simplex)
  row_sums  <- rowSums(mat)
  row_sums[row_sums == 0] <- 1  # avoid 0/0
  prop_mat  <- mat / row_sums;               dimnames(prop_mat) <- list(rn, cn)

  # Composite scores (combine specificity + magnitude)
  log_raw   <- log2(mat_raw + 1);             dimnames(log_raw)  <- list(rn, cn)
  comp1     <- l2_norm + log_raw;             dimnames(comp1)    <- list(rn, cn)
  comp2     <- l2_norm * log_raw;             dimnames(comp2)    <- list(rn, cn)

  # Tau (tissue specificity index)
  calc_tau  <- function(x) if (max(x) == 0) 0 else sum(1 - (x / max(x))) / (length(x) - 1)
  tau_vec   <- apply(mat, 1, calc_tau)
  tau_mat   <- matrix(rep(tau_vec, k), nrow = nrow(mat), ncol = k, dimnames = list(rn, cn))
  comp1_tau <- tau_mat + log_raw;             dimnames(comp1_tau) <- list(rn, cn)
  comp2_tau <- tau_mat * log_raw;             dimnames(comp2_tau) <- list(rn, cn)

  # Gini entropy on L2-normalised squared scores  H = -sum(x^2 * log(x^2))
  l2sq <- l2_norm^2
  row_entropy <- apply(l2sq, 1, function(x) {
    x <- x[x > 0]
    -sum(x * log(x))
  })
  gini_mat <- matrix(rep(row_entropy, k), nrow = nrow(mat), ncol = k, dimnames = list(rn, cn))

  # Effective number of cell types (inverse Simpson on L2-normalised scores)
  #   N_eff = 1 / sum(x_i^4)   where x_i are already L2-normed (sum(x^2)=1)
  row_neff <- 1 / rowSums(l2_norm^4)
  row_neff[!is.finite(row_neff)] <- 0
  neff_mat <- matrix(rep(row_neff, k), nrow = nrow(mat), ncol = k, dimnames = list(rn, cn))

  list(
    zscaled   = zscaled,
    logged    = log_raw,
    l2_norm   = l2_norm,
    prop      = prop_mat,
    comp1     = comp1,
    comp2     = comp2,
    tau       = tau_mat,
    comp1_tau = comp1_tau,
    comp2_tau = comp2_tau,
    gini      = gini_mat,
    neff      = neff_mat
  )
}

## Names of the metrics produced by .transform_single_mat (order preserved)
.METRIC_NAMES <- c(
  "zscaled","logged","l2_norm","prop","comp1","comp2",
  "tau","comp1_tau","comp2_tau","gini","neff"
)

compute_metrics <- function(mat, rn, cn) {
  mat_raw <- mat
  dimnames(mat_raw) <- list(rn, cn)

  # Quantile normalisation (CHEERS-inspired, equalises empirical distributions across cell types)
  mat_qn <- preprocessCore::normalize.quantiles(mat)
  dimnames(mat_qn) <- list(rn, cn)

  # --- Transformations on raw PeakMatrix ---
  raw_out <- .transform_single_mat(mat_raw, mat_raw, rn, cn)

  # --- Transformations on QN matrix (composites still use log2(raw+1) for magnitude) ---
  qn_out  <- .transform_single_mat(mat_qn, mat_raw, rn, cn)

  # --- Build output: raw_*, qn_*, and all ranks ---
  out <- list()

  # The raw PeakMatrix and QN matrix themselves
  out[["quantile_normalized"]] <- mat_qn

  for (nm in .METRIC_NAMES) {
    out[[paste0("raw_", nm)]]      <- raw_out[[nm]]
    out[[paste0("raw_", nm, "_rank")]] <- add_rank(raw_out[[nm]], rn, cn)
    out[[paste0("qn_",  nm)]]      <- qn_out[[nm]]
    out[[paste0("qn_",  nm, "_rank")]] <- add_rank(qn_out[[nm]], rn, cn)
  }

  out
}

## All derived assay names (auto-generated from .METRIC_NAMES)
DERIVED_ASSAYS <- c(
  "quantile_normalized",
  as.vector(outer(
    c("raw_", "qn_"),
    c(.METRIC_NAMES, paste0(.METRIC_NAMES, "_rank")),
    paste0
  ))
)



#### 2. RSE processing ####

# Assign the peak IDs
ensure_peakids <- function(se) {
  if (is.null(rownames(se)) || anyNA(rownames(se))) {
    rr <- rowRanges(se) # should all have rowRanges as a RSE object from ArchR normalised 
    rownames(se) <- paste0(as.character(seqnames(rr)), ":", start(rr), "-", end(rr))
  }
  if (!"peakid" %in% names(rowData(se))) rowData(se)$peakid <- rownames(se)
  se
}

# ArchR getGroupSE extension check for rowRanges peakID and genome Info
prepare_clear_se <- function(se, genome = "hg38") {
  if (!length(assayNames(se))) stop("No assay found in SE.")
  # Make sure to add the row ranges 
  rr <- rowRanges(se)
  if (is.null(rr) || length(rr) == 0L) {
    rd <- as.data.frame(rowData(se))
    if (all(c("seqnames", "start", "end") %in% names(rd))) {
      rr <- GRanges(seqnames = rd$seqnames, ranges = IRanges(start = rd$start, end = rd$end))
      rowRanges(se) <- rr
    } else {
      stop("No rowRanges and no usable seqnames/start/end in rowData.")
    }
  }
  
  if (is.null(rownames(se)) || anyNA(rownames(se)) || any(rownames(se) == "")) {
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
    last_is_bin <- n >= 3 && grepl("^bin-", elems[n])
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

build_se_library_if_needed <- function(
    liftover = TRUE, # Can liftover the annotations to hg19, to adapt the GWAS SNPs list in hg19   NOTE: all OpenTargets finemapped traits are in hg38
    drop_cts = DROP_CTS_DEFAULT, # User specified to drop certain cell types 
    se_library_root = SE_LIBRARY_ROOT,
    chain_file = CHAIN_FILE, # Only loaded the hg38tohg19 as scATAC libraries should be in hg38
    genome_from = "hg38",
    genome_to = "hg19",
    add_metrics = TRUE,
    metrics_overwrite = FALSE,
    force = FALSE
) {
  dir.create(se_library_root, recursive = TRUE, showWarnings = FALSE)
  if (liftover) chain <- import.chain(chain_file)
  
  for (input_name in names(input_files)) {
    se_path <- input_files[[input_name]]$path
    assays  <- input_files[[input_name]]$assays
    
    # Give different suffix for genome build but generally set to FALSE
    for (assay_name in assays) {
      message(">> Checking input: ", input_name, " | assay: ", assay_name)
      out_dir <- file.path(se_library_root, input_name, assay_name)
      out38   <- file.path(out_dir, paste0("SE_", input_name, "_", assay_name, ".hg38.rds"))
      out19   <- file.path(out_dir, paste0("SE_", input_name, "_", assay_name, ".hg19.rds"))
      
      if (!force) {
        if (liftover && file.exists(out38) && file.exists(out19)) {
          message("   ✓ Both hg38 & hg19 exist  skipping."); next
        }
        if (!liftover && file.exists(out38)) {
          message("   ✓ hg38 exists  skipping."); next
        }
      }
      
      message("   … building new SE(s)")
      se <- readRDS(se_path)
      se <- prepare_clear_se(se, genome = genome_from)
      
      if (length(drop_cts)) {
        keep <- setdiff(colnames(se), drop_cts)
        if (length(keep) < ncol(se)) {
          dropped <- setdiff(colnames(se), keep)
          se <- se[, keep, drop = FALSE]
          message("   … dropped: ", paste(dropped, collapse = ", "))
        }
      }
      
      if (add_metrics) {
        se <- add_metrics_to_se(se, overwrite = metrics_overwrite, verbose = TRUE)
      }
      
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(se, out38); message("   ✓ Saved hg38: ", out38)
      
      if (liftover) {
        rr38 <- rowRanges(se)
        grl  <- liftOver(rr38, chain)
        one  <- which(S4Vectors::elementNROWS(grl) == 1L)
        if (length(one) > 0L) {
          rr19 <- unlist(grl[one], use.names = FALSE)
          seqlevelsStyle(rr19) <- "UCSC"
          genome(rr19) <- genome_to
          se19 <- se[one, ]
          rowRanges(se19) <- rr19
          saveRDS(se19, out19); message("   ✓ Saved hg19: ", out19)
        } else {
          warning("No 1:1 mappings found for ", input_name, "/", assay_name)
        }
      }
    }
  }
  message("Library build complete.")
}


#### 3. Credible-set utilities ####

# GRanges  
detect_credset_cols <- function(credset) {
  cs_cols <- names(mcols(credset))
  pick_first <- function(cands) {
    hit <- cands[cands %in% cs_cols]
    if (length(hit)) hit[1] else NA_character_
  }
  list(
    ccv    = pick_first(c("names","SNP","rsid","rsID","variant","marker","ID")),
    signal = pick_first(c("signal","locus","cs","credible_set","signal_id","cluster","region"))
  )
}


#### 4. Overlap RSE with the credible set ####
# For a given trait, take the credible set (GRanges) and for all the RSE inthe SE library, find the overlaps --> subset of the peak matrix with CCV metadata

build_overlap_rse_for_library <- function(
    TRAIT_CODE, # Specified for trait name 
    credset_gr,
    se_index = index_se_library(genome = "hg38"),
    out_root = file.path("output/CLEAR_OVERLAPS", paste0(TRAIT_CODE, "__hg38")), # Only the hg38 now for OpenTargets
    assays_keep = NULL,
    force = FALSE
) {
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  ## detect metadata fields
  cols <- detect_credset_cols(credset_gr)
  cs_chr <- function(col) {
    if (is.na(col)) return(rep(NA_character_, length(credset_gr)))
    as.character(mcols(credset_gr)[[col]])
  }
  
  ## main loop for all the RSE indexed # but need to be specified in the beginning
  for (i in seq_len(nrow(se_index))) {  
    se_path    <- se_index$path[i]
    settings_id <- tools::file_path_sans_ext(basename(se_path))
    genome_tag <- if (grepl("\\.hg19(\\.se)?\\.rds$", se_path)) "hg19" else "hg38" # Regular expression to get useful namebases 
    
    out_dir <- file.path(out_root, settings_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    out_rds <- file.path(
      out_dir,
      paste0("OVERLAPS_", TRAIT_CODE, "__", settings_id, ".", genome_tag, ".rds")
    )
    
    if (!force && file.exists(out_rds)) {
      message(">> Skipping overlap RSE (exists): ", basename(out_rds))
      next
    }
    
    message(">> Building overlap RSE for: ", settings_id)
    
    ## Load SE and normalize identifiers
    se <- readRDS(se_path)
    se <- ensure_peakids(se)
    try(seqlevelsStyle(rowRanges(se)) <- "UCSC", silent = TRUE)
    
    N_total <- nrow(se) # need to keep a record of how many peaks are there for mean rank calculation, because we only keeping the overlapped subset of the original matrix
    
    ## Credible set overlap
    ov <- findOverlaps(se, credset_gr, ignore.strand = TRUE)
    
    if (length(ov) == 0L) {
      warning("  No overlaps for ", settings_id, "  saving empty RSE skeleton.")
      se_sub <- se[0, ]
      metadata(se_sub)$CLEAR <- list( # Keep some metadata, useful for naming and calculation, store everything required for CLEAR step two enrichment 
        trait          = TRAIT_CODE,
        settings_id    = settings_id,
        genome         = genome_tag,
        overlap_only   = TRUE,
        N_total        = N_total,
        N_overlap      = 0L,
        assays_present = assayNames(se_sub),
        created        = Sys.time()
      )
      saveRDS(se_sub, out_rds)
      next
    }
    
    ## Subset to overlapping peaks
    hit_rows <- sort(unique(queryHits(ov)))
    se_sub   <- se[hit_rows, , drop = FALSE]
    
    ## Map original SE indices to subset indices
    sub_map <- setNames(seq_along(hit_rows), hit_rows)
    
    ## Get those SE indices and respective matched CCV and signals 
    rsid_vec   <- cs_chr(cols$ccv)
    signal_vec <- cs_chr(cols$signal)
    
    kk <- data.frame(
      orig_row = queryHits(ov),
      i        = sub_map[as.character(queryHits(ov))],
      rsid     = rsid_vec[subjectHits(ov)],
      signal   = signal_vec[subjectHits(ov)],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::group_by(i) %>%
      dplyr::summarise(
        CCVs      = paste(unique(na.omit(rsid)), collapse = ","),
        signal    = paste(unique(na.omit(signal)), collapse = ","),
        n_CCVs    = dplyr::n_distinct(na.omit(rsid)),
        n_signals = dplyr::n_distinct(na.omit(signal)),
        .groups   = "drop"
      )
    
    ## Merge overlap metadata into rowData
    rd <- as.data.frame(rowData(se_sub))
    rd$peakid <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se_sub)
    rd$.__i__ <- seq_len(nrow(rd))   # sequential index for join
    
    rd <- dplyr::left_join(rd, kk, by = c(".__i__" = "i")) %>%
      dplyr::mutate(
        CCVs      = ifelse(is.na(CCVs), "", CCVs),
        signal    = ifelse(is.na(signal), "", signal),
        n_CCVs    = ifelse(is.na(n_CCVs), 0L, n_CCVs),
        n_signals = ifelse(is.na(n_signals), 0L, n_signals)
      ) %>%
      dplyr::select(-.__i__)  ## remove join key, no -i_orig needed
    
    rowData(se_sub) <- S4Vectors::DataFrame(rd, row.names = rownames(se_sub))
    
    ## Drop unwanted assays
    base_assays <- c("CPM", "PeakMatrix_perbp", "PeakMatrix")
    keep <- if (is.null(assays_keep)) {
      intersect(c(base_assays, DERIVED_ASSAYS), assayNames(se_sub))
    } else {
      intersect(assays_keep, assayNames(se_sub))
    }
    drop <- setdiff(assayNames(se_sub), keep)
    if (length(drop)) for (nm in drop) assays(se_sub)[[nm]] <- NULL
    
    metadata(se_sub)$CLEAR <- list(
      trait          = TRAIT_CODE,
      settings_id    = settings_id,
      genome         = genome_tag,
      overlap_only   = TRUE,
      N_total        = N_total,
      N_overlap      = nrow(se_sub),
      assays_present = assayNames(se_sub),
      created        = Sys.time()
    )
    
    saveRDS(se_sub, out_rds)
    message("  ✓ Saved overlap RSE: ", out_rds,
            " (rows=", nrow(se_sub), ", N_total=", N_total, ")")
  }
  
  message("All overlap RSEs written to: ", out_root)
}


#### 5. Mean-rank & plot DF from overlaps ####

meanrank_from_overlap <- function(se_overlap, rank_type) {
  if (!rank_type %in% assayNames(se_overlap)) {
    warning("Assay '", rank_type, "' not found in overlap SE.")
    return(NULL)
  }
  md <- metadata(se_overlap)$CLEAR
  if (is.null(md$N_total)) stop("Overlap SE missing CLEAR$N_total metadata.")
  
  N <- as.numeric(md$N_total)
  M <- SummarizedExperiment::assay(se_overlap, rank_type)
  n <- nrow(M)
  if (n == 0L) return(NULL)
  
  ## CHEERS adapted mean ranking test for global enrichment 
  mu    <- (N + 1) / 2
  sigma <- sqrt((N^2 - 1) / (12 * n))
  mr_vals <- if (inherits(M, "Matrix")) Matrix::colMeans(M) else colMeans(as.matrix(M))
  
  out <- data.frame(
    celltype           = colnames(M),
    observed_mean_rank = as.numeric(mr_vals),
    n_peaks            = n,
    exp_rank           = mu,
    stringsAsFactors   = FALSE, check.names = FALSE
  )
  out$obs.exp  <- out$observed_mean_rank / mu
  out$p        <- pnorm(out$observed_mean_rank, mu, sigma, lower.tail = FALSE)
  out$fdr      <- p.adjust(out$p, method = "BH")
  out$positive <- out$obs.exp > 1 & out$fdr < 0.05
  rownames(out) <- out$celltype
  out
}

make_plot_df_from_overlap <- function(se_overlap) {
  rd <- as.data.frame(rowData(se_overlap))
  if (!"peakid" %in% names(rd)) {
    rr <- rowRanges(se_overlap)
    rd$peakid <- paste0(as.character(seqnames(rr)), ":", start(rr), "-", end(rr))
  }
  rd$CCVs[is.na(rd$CCVs)]     <- ""
  rd$signal[is.na(rd$signal)] <- ""
  
  assays <- intersect(
    DERIVED_ASSAYS,
    assayNames(se_overlap)
  )
  if (!length(assays)) return(tibble::tibble())
  
  Reduce(dplyr::full_join, lapply(assays, function(a) {
    M <- SummarizedExperiment::assay(se_overlap, a)
    M <- suppressWarnings(as.matrix(M)); storage.mode(M) <- "double"
    df <- cbind(rd, as.data.frame(M, check.names = FALSE))
    long <- tidyr::pivot_longer(
      df,
      cols      = tidyselect::all_of(colnames(se_overlap)),
      names_to  = "celltype",
      values_to = a
    )
    long[[a]] <- suppressWarnings(as.numeric(long[[a]]))
    long[, c("peakid","CCVs","signal","celltype", a)]
  }))
}


# Stacked barplots of per-variant cell-type contributions for highly specific peaks
plot_celltype_stackbars <- function(
    se_overlap,
    value_metric    = "qn_l2_norm",
    value_threshold = 0.5,
    outfile,
    max_peaks_per_ct = 100
) {
  stopifnot(!missing(outfile))
  
  if (!value_metric %in% assayNames(se_overlap)) {
    warning("Assay '", value_metric, "' not found; skipping stack bar plots.")
    return(invisible(FALSE))
  }
  
  val_mat <- SummarizedExperiment::assay(se_overlap, value_metric)
  if (inherits(val_mat, "Matrix")) val_mat <- as.matrix(val_mat)
  storage.mode(val_mat) <- "double"
  
  # plot squared contributions so stacks sum to ~1 per peak (if L2-normalized per peak)
  val_plot <- val_mat^2
  
  # Labels
  rd <- as.data.frame(rowData(se_overlap))
  peakid     <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se_overlap)
  ccv_label  <- if ("CCVs"  %in% names(rd)) rd$CCVs  else ""
  sig_label  <- if ("signal" %in% names(rd)) rd$signal else ""
  base_label <- ifelse(nzchar(ccv_label), ccv_label, peakid)
  labels     <- ifelse(nzchar(sig_label), paste0(base_label, " (", sig_label, ")"), base_label)
  
  n_pages <- 0L
  grDevices::pdf(outfile, width = 10, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  for (ct in colnames(val_mat)) {
    
    # --- selection: peaks where THIS cell type's specificity score >= threshold ---
    sel <- which(val_mat[, ct] >= value_threshold)
    if (!length(sel)) next
    
    # cap plotted peaks (keep highest-scoring)
    if (length(sel) > max_peaks_per_ct) {
      ord <- order(val_mat[sel, ct], decreasing = TRUE)[seq_len(max_peaks_per_ct)]
      sel <- sel[ord]
    }
    
    submat <- val_plot[sel, , drop = FALSE]
    
    df <- as.data.frame(submat, check.names = FALSE)
    df$peak_label <- labels[sel]
    
    # order peaks by current ct's squared contribution (descending)
    target_vals <- df[[ct]]
    peak_levels <- df$peak_label[order(target_vals, decreasing = TRUE, na.last = NA)]
    if (!length(peak_levels)) peak_levels <- df$peak_label
    
    long <- tidyr::pivot_longer(
      df,
      cols      = tidyselect::all_of(colnames(val_plot)),
      names_to  = "celltype",
      values_to = "value"
    )
    
    long$peak_label <- factor(long$peak_label, levels = rev(unique(peak_levels)))
    long$highlight  <- long$celltype == ct
    
    # stack order: ct at bottom, others by total squared contribution
    ct_sums <- colSums(submat, na.rm = TRUE)
    other_levels <- names(sort(ct_sums, decreasing = TRUE))
    other_levels <- setdiff(other_levels, ct)
    long$celltype <- factor(long$celltype, levels = c(ct, other_levels))
    
    p <- ggplot(long) +
      geom_bar(
        aes(x = peak_label, y = value, fill = celltype, alpha = highlight),
        stat = "identity"
      ) +
      scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1), guide = "none") +
      labs(
        title = paste0(ct, " | ", value_metric, "\u00B2 | ", value_metric, " >= ", value_threshold),
        x     = "Peak / CCV",
        y     = paste0(value_metric, "\u00B2")
      ) +
      theme_Publication() +
      theme(
        axis.text.x  = element_text(angle = 60, hjust = 1, vjust = 1, size = 4),
        axis.text.y  = element_text(size = 6),
        legend.text  = element_text(size = 6),
        legend.title = element_text(size = 8),
        plot.title   = element_text(size = 11)
      )
    
    print(p)
    n_pages <- n_pages + 1L
  }
  
  if (n_pages == 0L) {
    warning("No cell types had peaks with ", value_metric, " >= ", value_threshold, "; PDF is empty.")
  }
  message("  ✓ Saved cell-type stackbars (", n_pages, " pages): ", outfile)
  invisible(n_pages > 0L)
}


# Peak specificity distribution: how many cell types is each peak "active" in at various thresholds
plot_peak_specificity_distribution <- function(
    se,
    metric       = "qn_l2_norm",
    thresholds   = seq(0, 0.9, 0.1),
    outfile      = NULL,
    title_prefix = ""
) {
  if (!metric %in% assayNames(se)) {
    warning("Assay '", metric, "' not found; skipping specificity distribution.")
    return(invisible(NULL))
  }
  
  mat <- SummarizedExperiment::assay(se, metric)
  if (inherits(mat, "Matrix")) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  n_ct <- ncol(mat)
  
  # For each threshold, count how many cell types each peak exceeds
  df_list <- lapply(thresholds, function(t) {
    k <- rowSums(mat > t)
    tab <- table(factor(k, levels = 0:n_ct))
    data.frame(
      threshold    = t,
      n_active_cts = as.integer(names(tab)),
      n_peaks      = as.integer(tab),
      stringsAsFactors = FALSE
    )
  })
  df <- dplyr::bind_rows(df_list)
  df$threshold_label <- factor(
    paste0("T=", df$threshold),
    levels = paste0("T=", sort(unique(df$threshold)))
  )
  
  p <- ggplot(df, aes(x = n_active_cts, y = n_peaks, fill = threshold_label)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ threshold_label, scales = "free_y") +
    labs(
      title = paste0(title_prefix, "Peak specificity distribution (", metric, ")"),
      x     = paste0("# cell types with ", metric, " > threshold"),
      y     = "# peaks"
    ) +
    theme_Publication() +
    theme(axis.text.x = element_text(size = 7))
  
  if (!is.null(outfile)) {
    grDevices::pdf(outfile, width = 12, height = 8)
    print(p)
    grDevices::dev.off()
    message("  ✓ Saved specificity distribution: ", outfile)
  }
  
  invisible(list(data = df, plot = p))
}


# Correlation plot triplets:
#   Panel 1: raw PeakMatrix (reference, col means per cell type)
#   Panel 2: transformed metric on PeakMatrix  (raw_*)
#   Panel 3: transformed metric on PeakMatrix_QN (qn_*)
# One PDF page per metric × cell-type-pair
plot_metric_correlations <- function(
    se,
    input_assay     = NULL,
    metrics         = c("l2_norm","prop","comp1","comp2","tau","gini","neff"),
    outfile         = NULL,
    sample_n        = 5000,    # down-sample peaks for speed
    title_prefix    = ""
) {
  # Determine the raw input assay used for compute_metrics
  if (is.null(input_assay)) input_assay <- choose_input_assay(se)
  if (!input_assay %in% assayNames(se)) {
    warning("Input assay '", input_assay, "' not found in SE; skipping correlation plots.")
    return(invisible(NULL))
  }

  mat_ref <- SummarizedExperiment::assay(se, input_assay)
  if (inherits(mat_ref, "Matrix")) mat_ref <- as.matrix(mat_ref)
  storage.mode(mat_ref) <- "double"

  # Subsample rows for plotting
  idx <- if (nrow(mat_ref) > sample_n) sort(sample.int(nrow(mat_ref), sample_n)) else seq_len(nrow(mat_ref))

  if (!is.null(outfile)) {
    grDevices::pdf(outfile, width = 15, height = 5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  n_pages <- 0L

  for (m in metrics) {
    raw_name <- paste0("raw_", m)
    qn_name  <- paste0("qn_",  m)

    has_raw <- raw_name %in% assayNames(se)
    has_qn  <- qn_name  %in% assayNames(se)
    if (!has_raw && !has_qn) next

    raw_mat <- if (has_raw) {
      mm <- SummarizedExperiment::assay(se, raw_name)
      if (inherits(mm, "Matrix")) mm <- as.matrix(mm); storage.mode(mm) <- "double"; mm
    } else NULL

    qn_mat <- if (has_qn) {
      mm <- SummarizedExperiment::assay(se, qn_name)
      if (inherits(mm, "Matrix")) mm <- as.matrix(mm); storage.mode(mm) <- "double"; mm
    } else NULL

    # Build tidy data for all cell types at once (vectorised)
    ref_sub <- mat_ref[idx, , drop = FALSE]
    cts     <- colnames(mat_ref)
    n_sub   <- length(idx)

    df_list <- lapply(cts, function(ct) {
      d <- data.frame(
        celltype = ct,
        ref      = ref_sub[, ct],
        stringsAsFactors = FALSE
      )
      d$raw_val <- if (!is.null(raw_mat)) raw_mat[idx, ct] else NA_real_
      d$qn_val  <- if (!is.null(qn_mat))  qn_mat[idx, ct]  else NA_real_
      d
    })
    df <- dplyr::bind_rows(df_list)

    # Panel 1: ref vs raw
    cor_raw <- if (!is.null(raw_mat)) {
      round(cor(as.vector(ref_sub), as.vector(raw_mat[idx, ]), use = "complete.obs"), 3)
    } else NA
    # Panel 2: ref vs qn
    cor_qn <- if (!is.null(qn_mat)) {
      round(cor(as.vector(ref_sub), as.vector(qn_mat[idx, ]), use = "complete.obs"), 3)
    } else NA
    # Panel 3: raw vs qn
    cor_rq <- if (!is.null(raw_mat) && !is.null(qn_mat)) {
      round(cor(as.vector(raw_mat[idx, ]), as.vector(qn_mat[idx, ]), use = "complete.obs"), 3)
    } else NA

    p1 <- ggplot(df, aes(x = ref, y = raw_val)) +
      geom_point(alpha = 0.15, size = 0.3) +
      facet_wrap(~ celltype, scales = "free") +
      labs(
        title = paste0(title_prefix, m, " | ", input_assay, " vs raw_", m, " (r=", cor_raw, ")"),
        x = input_assay, y = paste0("raw_", m)
      ) + theme_Publication()

    p2 <- ggplot(df, aes(x = ref, y = qn_val)) +
      geom_point(alpha = 0.15, size = 0.3) +
      facet_wrap(~ celltype, scales = "free") +
      labs(
        title = paste0(title_prefix, m, " | ", input_assay, " vs qn_", m, " (r=", cor_qn, ")"),
        x = input_assay, y = paste0("qn_", m)
      ) + theme_Publication()

    p3 <- ggplot(df, aes(x = raw_val, y = qn_val)) +
      geom_point(alpha = 0.15, size = 0.3) +
      facet_wrap(~ celltype, scales = "free") +
      labs(
        title = paste0(title_prefix, m, " | raw_", m, " vs qn_", m, " (r=", cor_rq, ")"),
        x = paste0("raw_", m), y = paste0("qn_", m)
      ) + theme_Publication()

    print(p1); print(p2); print(p3)
    n_pages <- n_pages + 3L
  }

  if (n_pages == 0L) warning("No matching assay pairs found for correlation plots.")
  message("  ✓ Saved metric correlation plots (", n_pages, " pages): ", outfile)
  invisible(n_pages > 0L)
}


###############################
## 6. Main CLEAR-from-overlaps runner
###############################

run_trait_from_overlap_library <- function(
    TRAIT_CODE,
    overlap_root   = file.path("output/CLEAR_OVERLAPS", paste0(TRAIT_CODE, "__hg38")),
    out_root       = file.path("output/CLEAR_FROM_OVERLAPS", paste0(TRAIT_CODE, "__hg38")),
    rank_types     = c("qn_l2_norm_rank","qn_comp1_rank","qn_comp2_rank","qn_tau_rank","qn_comp1_tau_rank","qn_comp2_tau_rank"),
    spec_metric    = "qn_l2_norm",
    conc_threshold = 0.8,   # sum of L2^2 across the top-k cell types must exceed this
    max_k          = 4,     # check k = 1, 2, ..., max_k dominant cell types
    force          = FALSE
) {
  stopifnot(dir.exists(overlap_root))
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  paths <- list.files(
    overlap_root,
    pattern = "^OVERLAPS_.*\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(paths)) {
    warning("No overlap RSEs found.")
    return(invisible())
  }
  
  for (p in paths) {
    se <- readRDS(p)
    md <- metadata(se)$CLEAR
    settings_id <- md$settings_id
    trait_dir   <- file.path(out_root, settings_id)
    dir.create(trait_dir, recursive = TRUE, showWarnings = FALSE)
    input_id <- paste0(TRAIT_CODE, "__", settings_id)
    
    final_file <- file.path(trait_dir, paste0(input_id, "_summary.tsv"))
    if (!force && file.exists(final_file)) {
      message(">> Skipping (already finished): ", settings_id)
      next
    }
    
    message(">> Running from OVERLAPS | SE=", settings_id)
    
    ## 1. Mean Rank Test & PDF generation
    pdf_path <- file.path(trait_dir, paste0(input_id, "_meanrank.pdf"))
    grDevices::pdf(pdf_path, width = 10, height = 5)
    
    core_enriched_cts <- character(0)
    all_meanrank <- list()
    
    for (rt in rank_types) {
      mr <- meanrank_from_overlap(se, rt)
      if (is.null(mr)) next
      
      thresh  <- -log10(0.05 / max(1, nrow(mr)))
      sig_idx <- which(-log10(mr$p) >= thresh & !is.na(mr$p))
      if (length(sig_idx)) {
        core_enriched_cts <- union(core_enriched_cts, mr$celltype[sig_idx])
      }
      
      p_ <- ggplot(mr) +
        geom_hline(yintercept = thresh, linetype = "dashed", color = "gray40") +
        geom_bar(
          aes(
            x    = reorder(celltype, -log10(p)),
            y    = -log10(p),
            fill = (-log10(p) >= thresh)
          ),
          stat = "identity"
        ) +
        scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "firebrick")) +
        coord_flip() +
        ggtitle(paste(rt, "|", input_id)) +
        theme_Publication() +
        theme(legend.position = "none")
      
      print(p_)
      mr$rank_type <- rt
      all_meanrank[[rt]] <- mr
    }
    
    grDevices::dev.off()  # <- clean close of the PDF
    if (length(all_meanrank)) {
      readr::write_tsv(
        dplyr::bind_rows(all_meanrank),
        file.path(trait_dir, paste0(input_id, "_meanrank.tsv"))
      )
    }
    
    ## 2. Specificity-based peak–cell-type assignments
    ##    For each peak, find if k = 1..max_k cell types concentrate >= conc_threshold
    ##    of the L2^2 signal (i.e. sum of top-k squared scores >= conc_threshold),
    ##    with each contributing cell type scoring >= sqrt(conc_threshold / k).
    
    if (!spec_metric %in% assayNames(se)) {
      warning("Assay '", spec_metric, "' not found; skipping specificity summaries for ", settings_id)
      readr::write_tsv(tibble::tibble(), final_file)
      next
    }
    
    spec_mat <- SummarizedExperiment::assay(se, spec_metric)
    if (inherits(spec_mat, "Matrix")) spec_mat <- as.matrix(spec_mat)
    storage.mode(spec_mat) <- "double"
    spec_sq <- spec_mat^2  # squared L2 norms; rows sum to ~1
    
    rd <- as.data.frame(rowData(se))
    peakids  <- if ("peakid" %in% names(rd)) rd$peakid else rownames(se)
    ccv_col  <- if ("CCVs"   %in% names(rd)) rd$CCVs   else rep("", nrow(se))
    sig_col  <- if ("signal"  %in% names(rd)) rd$signal  else rep("", nrow(se))
    
    specific_peaks <- list()
    
    for (k in seq_len(max_k)) {
      min_per_ct <- sqrt(conc_threshold / k)
      
      for (i in seq_len(nrow(spec_sq))) {
        row_sorted <- sort(spec_sq[i, ], decreasing = TRUE)
        top_k_sum  <- sum(row_sorted[seq_len(min(k, length(row_sorted)))])
        
        if (top_k_sum < conc_threshold) next
        
        # which cell types pass the per-ct threshold?
        passing_cts <- names(which(spec_mat[i, ] >= min_per_ct))
        if (length(passing_cts) < 1L || length(passing_cts) > k) next
        
        for (ct in passing_cts) {
          specific_peaks[[length(specific_peaks) + 1L]] <- data.frame(
            peakid        = peakids[i],
            celltype      = ct,
            signal        = sig_col[i],
            CCVs          = ccv_col[i],
            l2_score      = spec_mat[i, ct],
            l2_sq         = spec_sq[i, ct],
            top_k_sq_sum  = top_k_sum,
            k             = k,
            min_per_ct    = min_per_ct,
            core_enriched = ct %in% core_enriched_cts,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    if (length(specific_peaks)) {
      spec_df <- dplyr::bind_rows(specific_peaks) |>
        dplyr::distinct(peakid, celltype, .keep_all = TRUE)  # keep smallest k
      
      readr::write_tsv(
        spec_df,
        file.path(trait_dir, paste0(input_id, "_specific_peaks.tsv"))
      )
      
      ## Per-k summary
      spec_summary <- spec_df |>
        dplyr::group_by(k) |>
        dplyr::summarise(
          n_peaks        = dplyr::n_distinct(peakid),
          n_celltypes    = dplyr::n_distinct(celltype),
          n_signals      = dplyr::n_distinct(signal[nzchar(signal)]),
          n_CCVs         = dplyr::n_distinct(CCVs[nzchar(CCVs)]),
          n_core_ct      = sum(core_enriched),
          n_residue_ct   = sum(!core_enriched),
          .groups        = "drop"
        )
      readr::write_tsv(spec_summary, final_file)
    } else {
      readr::write_tsv(
        tibble::tibble(
          k = seq_len(max_k), n_peaks = 0L, n_celltypes = 0L,
          n_signals = 0L, n_CCVs = 0L, n_core_ct = 0L, n_residue_ct = 0L
        ),
        final_file
      )
    }
    
    ## 3. Cell-type stacked bar plots for peaks with high specificity scores
    stack_metrics <- c("qn_l2_norm", "qn_comp1", "qn_comp2")
    for (vm in stack_metrics) {
      stack_pdf <- file.path(trait_dir, paste0(input_id, "_", vm, "_stackbars.pdf"))
      try(
        plot_celltype_stackbars(
          se_overlap      = se,
          value_metric    = vm,
          value_threshold = 0.5,
          outfile         = stack_pdf
        ),
        silent = FALSE
      )
    }
    
    ## 4. Peak specificity distribution
    spec_pdf <- file.path(trait_dir, paste0(input_id, "_specificity_distribution.pdf"))
    try(
      plot_peak_specificity_distribution(
        se           = se,
        metric       = "qn_l2_norm",
        outfile      = spec_pdf,
        title_prefix = paste0(input_id, " | ")
      ),
      silent = FALSE
    )

    ## 5. Metric correlation triplets (raw PeakMatrix vs raw_* vs qn_*)
    corr_pdf <- file.path(trait_dir, paste0(input_id, "_metric_correlations.pdf"))
    try(
      plot_metric_correlations(
        se           = se,
        outfile      = corr_pdf,
        title_prefix = paste0(input_id, " | ")
      ),
      silent = FALSE
    )
    
    message(">> Finished (from overlaps) | SE=", settings_id)
  }
}




###############################
## 7. Config & run block
###############################

## 1) Build SE library (if not already done)
build_se_library_if_needed(
  liftover          = FALSE,   # change to TRUE if you want hg19 as well
  add_metrics       = TRUE,
  metrics_overwrite = FALSE,
  force             = FALSE
)

GENOME <- "hg38"

## Choose one trait at a time by setting TRAIT_CODE and CREDSET_FILEMAP
TRAIT_CODE <- "bcac_overall_opentarget"
CREDSET_FILEMAP <- c(
  bcac_overall_opentarget = "data/opentarget_trait/breast_carcinoma_2017_michailidou_cs210__GCST004988/breast_carcinoma_2017_michailidou_cs210__GCST004988_CLEAR_credible_set_variants.gr.rds"
)

TRAIT_CODE <- "bcac_overall_FM_Fachal_2020"
CREDSET_FILEMAP <- c(
  bcac_overall_FM_Fachal_2020 = "data/BCAC_FM_GR.rds"
)

TRAIT_CODE <- "t1d_santiago"
CREDSET_FILEMAP <- c(
  t1d_santiago = "data/opentarget_trait/working/t1d_santiago.rds"
)

TRAIT_CODE <- "JT_interval_GCST90179157"
CREDSET_FILEMAP <- c(
  JT_interval_GCST90179157 = "data/opentarget_trait/_gwas_by_bin/100_200/GCST90179157/GCST90179157_CLEAR_credible_set_variants.gr.rds"
)

TRAIT_CODE <- "FINNGEN_R12_I9_AF"
CREDSET_FILEMAP <- c(
  FINNGEN_R12_I9_AF= "data/opentarget_trait/_gwas_by_bin/100_200/FINNGEN_R12_I9_AF/FINNGEN_R12_I9_AF_CLEAR_credible_set_variants.gr.rds"
)

TRAIT_CODE <- "AD_GCST90301303"
CREDSET_FILEMAP <- c(
  AD_GCST90301303= "data/opentarget_trait/_gwas_by_bin/20_50/GCST90301303/GCST90301303_CLEAR_credible_set_variants.gr.rds"
)

se_index <- index_se_library(genome = GENOME)

credset_path <- CREDSET_FILEMAP[[TRAIT_CODE]]
stopifnot(!is.null(credset_path), file.exists(credset_path))

credset_gr <- readRDS(credset_path)
GenomeInfoDb::seqlevelsStyle(credset_gr) <- "UCSC"

## 2) Build overlap RSEs
overlap_root <- file.path("output/CLEAR_OVERLAPS", paste0(TRAIT_CODE, "__", GENOME))
build_overlap_rse_for_library(
  TRAIT_CODE = TRAIT_CODE,
  credset_gr = credset_gr,
  se_index   = se_index,
  out_root   = overlap_root,
  assays_keep= NULL,   # keep base + derived
  force      = FALSE
)

## 3) Run CLEAR from overlaps
out_root <- file.path("output/CLEAR_FROM_OVERLAPS", paste0(TRAIT_CODE, "__", GENOME))
run_trait_from_overlap_library(
  TRAIT_CODE     = TRAIT_CODE,
  overlap_root   = overlap_root,
  out_root       = out_root,
  rank_types     = c("qn_l2_norm_rank","qn_comp1_rank","qn_comp2_rank"),
  spec_metric    = "qn_l2_norm",
  conc_threshold = 0.8,
  max_k          = 4,
  force          = TRUE
)

message("CLEAR run from overlaps completed for TRAIT=", TRAIT_CODE, " | genome=", GENOME)
