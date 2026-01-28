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

SE_LIBRARY_ROOT <- "output/snATAC_SE_library"

input_files <- list(
  union_peak = list(
    path   = "data/PeakCalls_Xiong/union_GRSE.rds",
    assays = c("PeakMatrix", "union_bw_signal_perbp")
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

compute_metrics <- function(mat, rn, cn) {
  mat_qn <- preprocessCore::normalize.quantiles(mat)  # quantile normalisation based on CHEERS (keep empirical distri across cell types)
  dimnames(mat_qn) <- list(rn, cn)
  mat <- mat_qn # use this normalised matrix for all other specificity calculations 
  
  scaled_mat  <- t(apply(mat, 1, scale));       dimnames(scaled_mat)  <- list(rn, cn) # Z scaling
  l2norm_mat  <- t(apply(mat, 1, safe_l2));     dimnames(l2norm_mat)  <- list(rn, cn) # l2 !!!!!!
  compsc1_mat <- l2norm_mat + log2(mat + 1);    dimnames(compsc1_mat) <- list(rn, cn) # Composite calculations 
  compsc2_mat <- l2norm_mat * log2(mat + 1);    dimnames(compsc2_mat) <- list(rn, cn)
  
  # Tau calclulations 
  calc_tau <- function(x) if (max(x) == 0) 0 else sum(1 - (x / max(x))) / (length(x) - 1)
  tau_vec  <- apply(mat, 1, calc_tau)
  tau_mat  <- matrix(
    rep(tau_vec, ncol(mat)),
    nrow = nrow(mat), ncol = ncol(mat),
    dimnames = list(rn, cn)
  )
  
  compsc1_tau <- tau_mat + log2(mat + 1); dimnames(compsc1_tau) <- list(rn, cn)
  compsc2_tau <- tau_mat * log2(mat + 1); dimnames(compsc2_tau) <- list(rn, cn)
  
  list(
    zscaled             = scaled_mat,
    l2_norm             = l2norm_mat,
    comp1               = compsc1_mat,
    comp2               = compsc2_mat,
    spec_rank           = add_rank(l2norm_mat,  rn, cn),
    comp1_rank          = add_rank(compsc1_mat, rn, cn),
    comp2_rank          = add_rank(compsc2_mat, rn, cn),
    tau                 = tau_mat,
    tau_rank            = add_rank(tau_mat, rn, cn),
    compsc1_tau         = compsc1_tau,
    compsc1_tau_rank    = add_rank(compsc1_tau, rn, cn),
    compsc2_tau         = compsc2_tau,
    compsc2_tau_rank    = add_rank(compsc2_tau, rn, cn)
  )
}

DERIVED_ASSAYS <- c(
  "zscaled","l2_norm","comp1","comp2","spec_rank","comp1_rank","comp2_rank",
  "tau","tau_rank","compsc1_tau","compsc1_tau_rank","compsc2_tau","compsc2_tau_rank"
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
    c("l2_norm","spec_rank","comp1","comp1_rank","comp2","comp2_rank",
      "tau","tau_rank","compsc1_tau","compsc1_tau_rank","compsc2_tau","compsc2_tau_rank"),
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


###############################
## 6. Main CLEAR-from-overlaps runner
###############################

run_trait_from_overlap_library <- function(
    TRAIT_CODE,
    overlap_root = file.path("output/CLEAR_OVERLAPS", paste0(TRAIT_CODE, "__hg38")),
    out_root     = file.path("output/CLEAR_FROM_OVERLAPS", paste0(TRAIT_CODE, "__hg38")),
    rank_types   = c("spec_rank","comp1_rank","comp2_rank","tau_rank","compsc1_tau_rank","compsc2_tau_rank"),
    score_types  = c("spec_rank","comp1","comp2","tau_rank","compsc1_tau_rank","compsc2_tau_rank"),
    top_cutoff   = 0.10,
    force        = FALSE
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
    
    final_file <- file.path(trait_dir, paste0(input_id, "_highrank_summary.tsv"))
    if (!force && file.exists(final_file)) {
      message(">> Skipping (already finished): ", settings_id)
      next
    }
    
    message(">> Running from OVERLAPS | SE=", settings_id)
    
    ## 1. Mean Rank Test & PDF generation
    pdf_path <- file.path(trait_dir, paste0(input_id, "_meanrank.pdf"))
    grDevices::pdf(pdf_path, width = 10, height = 5)
    
    core_enriched_cts <- character(0)
    
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
      readr::write_tsv(
        mr,
        file.path(trait_dir, paste0(input_id, "_", rt, "_meanrank.tsv"))
      )
    }
    
    grDevices::dev.off()  # <- clean close of the PDF
    
    ## 2. High-Rank Summaries
    plot_df <- make_plot_df_from_overlap(se) |>
      dplyr::filter(!is.na(signal) & nzchar(signal))
    
    # --- CRITICAL CHECK 1: any rows with signal? ---
    if (nrow(plot_df) == 0) {
      warning("Filtered plot_df is empty (no valid signal ID). Skipping high-rank summaries for ", settings_id)
      
      score_summary <- dplyr::bind_rows(
        lapply(
          score_types,
          function(st)
            tibble::tibble(
              score_type     = st,
              n_CCVs         = 0L,
              n_signals      = 0L,
              n_peaks        = 0L,
              n_celltypes    = 0L,
              residue_signal = 0L,
              core_ct_signal = 0L
            )
        )
      )
      readr::write_tsv(score_summary, final_file)
      next
    }
    
    summary_list <- list()
    
    for (score_type in score_types) {
      if (!score_type %in% names(plot_df)) next
      
      df <- plot_df |>
        dplyr::select(
          peakid, celltype, signal, CCVs,
          rank = !!rlang::sym(score_type)
        ) |>
        dplyr::mutate(rank = suppressWarnings(as.numeric(rank))) |>
        dplyr::group_by(celltype) |>
        dplyr::mutate(
          thresh   = if (all(is.na(rank))) NA_real_
          else stats::quantile(rank, probs = top_cutoff, na.rm = TRUE),
          highrank = ifelse(is.na(thresh), FALSE, rank <= thresh)
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(core_enriched = celltype %in% core_enriched_cts)
      
      readr::write_tsv(
        df,
        file.path(trait_dir, paste0(input_id, "_", score_type, "_highrank_all.tsv"))
      )
      
      df_top <- dplyr::filter(df, highrank)
      
      # --- CRITICAL CHECK 2: any high-rank peaks? ---
      if (nrow(df_top) == 0L) {
        summary_list[[score_type]] <- tibble::tibble(
          score_type     = score_type,
          n_CCVs         = 0L,
          n_signals      = 0L,
          n_peaks        = 0L,
          n_celltypes    = 0L,
          residue_signal = 0L,
          core_ct_signal = 0L
        )
        next
      }
      
      per_signal <- df_top |>
        dplyr::group_by(signal) |>
        dplyr::summarise(
          n_peaks     = dplyr::n_distinct(peakid),
          any_core_ct = any(core_enriched),
          .groups     = "drop"
        )
      
      summary_list[[score_type]] <- tibble::tibble(
        score_type     = score_type,
        n_CCVs         = dplyr::n_distinct(df_top$CCVs),
        n_signals      = dplyr::n_distinct(df_top$signal),
        n_peaks        = dplyr::n_distinct(df_top$peakid),
        n_celltypes    = dplyr::n_distinct(df_top$celltype),
        residue_signal = sum(!per_signal$any_core_ct),
        core_ct_signal = sum(per_signal$any_core_ct)
      )
      
      readr::write_tsv(
        df_top,
        file.path(trait_dir, paste0(input_id, "_", score_type, "_highrank_peaks.tsv"))
      )
    }
    
    score_summary <- dplyr::bind_rows(summary_list)
    readr::write_tsv(score_summary, final_file)
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
  TRAIT_CODE   = TRAIT_CODE,
  overlap_root = overlap_root,
  out_root     = out_root,
  rank_types   = c("spec_rank","comp1_rank","comp2_rank","tau_rank","compsc1_tau_rank","compsc2_tau_rank"),
  score_types  = c("spec_rank","comp1","comp2","tau_rank","compsc1_tau_rank","compsc2_tau_rank"),
  top_cutoff   = 0.10,
  force        = TRUE
)

message("CLEAR run from overlaps completed for TRAIT=", TRAIT_CODE, " | genome=", GENOME)