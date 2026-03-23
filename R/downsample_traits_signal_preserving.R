#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(S4Vectors)
  library(ggplot2)
})

# in the 19_03 run of CLEAR within signal coherent analysis for breast_x_t2d and breast_x_heart_failure, there are a lot of significantly
# coherent loci identified. Consider testing for:
# 1. What are the heart failure loci that are coherent in the heart compare to those in breast? Similarly compare for the BRCA status in breast!
# 2. Global enrichment was wildly significant in many cell types in the t2d across many tissues, breast/heart/kidney, suspect this is the case
#    for traits with loads of SNPs, global enrichment not finding biology, but the polygenic nature of the trait. Examine the GWAS vector
#     a. SNPs per signal 
#     b. enrichment is driven by those loci with many SNPs. Rank is column-wise, suspect too many promoter GWAS peaks?

# Prefer traits with enough loci (for example >20 <500 (but new breast would have more than 500) signals),
# more importantly NEED TO AVOID traits with many huge signals (for example lots of loci with >500-700 SNPs),
# then use CLEAR to prioritize coherent loci/cell-type programs.

# Justify to control for domination by very large loci, LD-driven inflation (capping, stratification, weighting)
# But seems a bit too artificial and may break some per-locus SNP architecture (LD block we potentially want to fix?)

# The key is to do sensitivity analysis across multiple schedules and make sure the downsampling is stable. 

# ---------------------------
# Runtime config (RStudio)
# ---------------------------
# Input for the height 
interactive_cfg <- list(
  input_dir = "/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/height_finngen_cs1903__FINNGEN_R12_HEIGHT_IRN",
  output_dir = "/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/downsampling_test_height",
  file_pattern = "_CLEAR_credible_set_variants\\.gr\\.rds$",
  target_fraction = 0.5, # this is the downsampling fraction. keep half of the SNPs in each signal
  min_per_signal = 20L, # minimum number of SNPs per signal. But this way only downsampling signal with more than 20 SNPs
  max_per_signal = 200, # constrain highly heterogenous loci and LD block
  seed = 1L,
  schedule = "medium", # one of: global, mild, medium, strong
  n_repeats = 3L
)

usage <- paste(
  "Usage:",
  "Rscript downsample_traits_signal_preserving.R <input_dir> <output_dir>",
  "[file_pattern='\\\\.rds$'] [target_fraction=0.5] [min_per_signal=1] [max_per_signal=Inf] [seed=1] [schedule='medium'] [n_repeats=1]",
  sep = "\n"
)

parse_runtime_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # CLI mode (bash / terminal)
  if (length(args) > 0) {
    if (length(args) < 2) {
      stop(usage)
    }
    
    max_per_signal_raw <- if (length(args) >= 6) args[[6]] else "Inf"
    max_per_signal <- suppressWarnings(as.numeric(max_per_signal_raw))
    if (is.na(max_per_signal)) max_per_signal <- Inf
    
    cfg <- list(
      input_dir = args[[1]],
      output_dir = args[[2]],
      file_pattern = if (length(args) >= 3) args[[3]] else "\\.rds$",
      target_fraction = if (length(args) >= 4) as.numeric(args[[4]]) else 0.5,
      min_per_signal = if (length(args) >= 5) as.integer(args[[5]]) else 1L,
      max_per_signal = max_per_signal,
      seed = if (length(args) >= 7) as.integer(args[[7]]) else 1L,
      schedule = if (length(args) >= 8) as.character(args[[8]]) else "medium",
      n_repeats = if (length(args) >= 9) as.integer(args[[9]]) else 1L
    )
    return(cfg)
  }
  
  # Interactive mode (RStudio)
  if (interactive()) {
    if (!nzchar(interactive_cfg$input_dir) || !nzchar(interactive_cfg$output_dir)) {
      stop(
        "Interactive mode detected. Please set interactive_cfg$input_dir and interactive_cfg$output_dir at the top of the script."
      )
    }
    return(interactive_cfg)
  }
  
  stop(usage)
}

cfg <- parse_runtime_args()

input_dir <- cfg$input_dir
output_dir <- cfg$output_dir
file_pattern <- cfg$file_pattern
target_fraction <- as.numeric(cfg$target_fraction)
min_per_signal <- as.integer(cfg$min_per_signal)
max_per_signal <- as.numeric(cfg$max_per_signal)
seed <- as.integer(cfg$seed)
schedule <- as.character(cfg$schedule)
n_repeats <- as.integer(cfg$n_repeats)

if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (is.na(target_fraction) || target_fraction <= 0 || target_fraction > 1) {
  stop("target_fraction must be in (0, 1].")
}
if (is.na(min_per_signal) || min_per_signal < 1L) {
  stop("min_per_signal must be >= 1.")
}
if (is.na(seed)) {
  stop("seed must be an integer.")
}
if (is.na(n_repeats) || n_repeats < 1L) {
  stop("n_repeats must be >= 1.")
}
if (is.na(max_per_signal)) {
  max_per_signal <- Inf
}
if (!is.infinite(max_per_signal) && max_per_signal < min_per_signal) {
  stop("max_per_signal must be >= min_per_signal or Inf.")
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

SNP_BIN_LEVELS <- c("1-5", "5-10", "10-20", "20-100", "100-200", "200-400", "400-600", ">600")

SCHEDULE_TABLE <- data.frame(
  schedule = rep(c("mild", "medium", "strong"), each = length(SNP_BIN_LEVELS)),
  snp_bin = rep(SNP_BIN_LEVELS, times = 3),
  keep_fraction = c(
    # mild
    1.00, 0.90, 0.75, 0.60, 0.45, 0.35, 0.25, 0.20,
    # medium
    1.00, 0.80, 0.65, 0.50, 0.35, 0.25, 0.18, 0.12,
    # strong
    1.00, 0.70, 0.50, 0.35, 0.25, 0.18, 0.12, 0.08
  ),
  stringsAsFactors = FALSE
)

if (!schedule %in% c("global", unique(SCHEDULE_TABLE$schedule))) {
  stop("schedule must be one of: global, mild, medium, strong")
}

extract_signal <- function(gr) {
  mc <- S4Vectors::mcols(gr)
  if (!"signal" %in% colnames(mc)) {
    stop("Object is missing mcols(signal).")
  }
  sig <- as.character(mc$signal)
  sig[is.na(sig) | !nzchar(sig)] <- "UNASSIGNED_SIGNAL"
  sig
}

summarize_locus_counts <- function(signal_vec) {
  locus_counts <- as.integer(table(signal_vec))
  data.frame(
    n_signals = length(locus_counts),
    n_snps = length(signal_vec),
    snps_per_signal_mean = mean(locus_counts),
    snps_per_signal_median = as.numeric(stats::median(locus_counts)),
    snps_per_signal_q90 = as.numeric(stats::quantile(locus_counts, probs = 0.9, names = FALSE)),
    snps_per_signal_q95 = as.numeric(stats::quantile(locus_counts, probs = 0.95, names = FALSE)),
    snps_per_signal_max = max(locus_counts)
  )
}

compute_keep_count <- function(n_i, f, min_keep, max_keep) {
  raw <- n_i * f
  base <- floor(raw)
  p_up <- raw - base
  bump <- stats::rbinom(length(n_i), size = 1L, prob = p_up)
  keep <- base + bump
  
  keep <- pmax(keep, pmin(min_keep, n_i))
  if (!is.infinite(max_keep)) {
    keep <- pmin(keep, pmin(max_keep, n_i))
  }
  keep
}

downsample_signal_preserving <- function(gr, target_frac, min_keep, max_keep, seed_val, schedule_name) {
  set.seed(seed_val)
  
  signal_vec <- extract_signal(gr)
  idx_by_signal <- split(seq_along(signal_vec), signal_vec)
  n_by_signal <- vapply(idx_by_signal, length, integer(1))
  frac_by_signal <- fraction_by_signal(n_by_signal, schedule_name, target_frac)
  names(frac_by_signal) <- names(n_by_signal)
  
  keep_by_signal <- compute_keep_count(
    n_i = n_by_signal,
    f = frac_by_signal,
    min_keep = min_keep,
    max_keep = max_keep
  )
  
  keep_idx <- integer(0)
  for (sig in names(idx_by_signal)) {
    idx <- idx_by_signal[[sig]]
    k <- keep_by_signal[[sig]]
    if (k >= length(idx)) {
      keep_idx <- c(keep_idx, idx)
    } else {
      keep_idx <- c(keep_idx, sample(idx, size = k, replace = FALSE))
    }
  }
  
  keep_idx <- sort(unique(keep_idx))
  gr_down <- gr[keep_idx]
  
  list(
    gr = gr_down,
    keep_idx = keep_idx,
    keep_by_signal = keep_by_signal,
    n_by_signal = n_by_signal,
    frac_by_signal = frac_by_signal
  )
}

snp_bin_requested <- function(x) {
  out <- rep(NA_character_, length(x))
  out[x >= 1 & x <= 5] <- "1-5"
  out[x > 5 & x <= 10] <- "5-10"
  out[x > 10 & x <= 20] <- "10-20"
  out[x > 20 & x <= 100] <- "20-100"
  out[x > 100 & x <= 200] <- "100-200"
  out[x > 200 & x <= 400] <- "200-400"
  out[x > 400 & x <= 600] <- "400-600"
  out[x > 600] <- ">600"
  factor(out, levels = SNP_BIN_LEVELS)
}

fraction_by_signal <- function(n_by_signal, schedule_name, global_fraction) {
  if (schedule_name == "global") {
    return(rep(global_fraction, length(n_by_signal)))
  }
  
  bins <- as.character(snp_bin_requested(n_by_signal))
  sched <- SCHEDULE_TABLE[SCHEDULE_TABLE$schedule == schedule_name, c("snp_bin", "keep_fraction")]
  f_map <- setNames(sched$keep_fraction, sched$snp_bin)
  f <- unname(f_map[bins])
  f[is.na(f)] <- global_fraction
  f
}

build_dist_requested <- function(signal_orig, signal_down, trait_tag, rep_tag, schedule) {
  b_orig <- table(snp_bin_requested(as.integer(table(signal_orig))))
  b_down <- table(snp_bin_requested(as.integer(table(signal_down))))
  
  dist_requested <- rbind(
    data.frame(
      trait = trait_tag,
      set = "original",
      snp_bin = SNP_BIN_LEVELS,
      n_signals = as.integer(b_orig[SNP_BIN_LEVELS]),
      replicate = rep_tag,
      schedule = schedule,
      stringsAsFactors = FALSE
    ),
    data.frame(
      trait = trait_tag,
      set = "downsampled",
      snp_bin = SNP_BIN_LEVELS,
      n_signals = as.integer(b_down[SNP_BIN_LEVELS]),
      replicate = rep_tag,
      schedule = schedule,
      stringsAsFactors = FALSE
    )
  )
  dist_requested$n_signals[is.na(dist_requested$n_signals)] <- 0L
  dist_requested
}

plot_signal_size_bins_dodged <- function(dist_requested, out_pdf, trait_tag) {
  plot_df <- dist_requested
  plot_df$snp_bin <- factor(
    plot_df$snp_bin,
    levels = SNP_BIN_LEVELS
  )
  plot_df$set <- factor(plot_df$set, levels = c("original", "downsampled"))
  
  p <- ggplot(plot_df, aes(x = snp_bin, y = n_signals, fill = set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_text(
      aes(label = n_signals),
      position = position_dodge(width = 0.8),
      vjust = -0.25,
      size = 3
    ) +
    theme_bw() +
    scale_fill_manual(values = c(original = "#4C78A8", downsampled = "#F58518")) +
    labs(
      title = paste0("Signal size distribution (dodged): ", trait_tag),
      x = "SNPs per signal bin",
      y = "Number of signals",
      fill = "Set"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(out_pdf, p, width = 9, height = 5)
}

trait_files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE, recursive = TRUE)
trait_files <- sort(trait_files)
if (!length(trait_files)) {
  stop("No files matched pattern in input_dir: ", file_pattern)
}

summary_rows <- list()
distribution_rows <- list()

for (i in seq_along(trait_files)) {
  f <- trait_files[[i]]
  trait_dir_name <- basename(dirname(f))
  if (nzchar(trait_dir_name)) {
    trait_tag <- sanitize_name(trait_dir_name)
  } else {
    trait_tag <- sanitize_name(tools::file_path_sans_ext(basename(f)))
  }
  trait_out_dir <- file.path(output_dir, trait_tag)
  dir.create(trait_out_dir, recursive = TRUE, showWarnings = FALSE)
  message("Processing [", i, "/", length(trait_files), "]: ", basename(f))
  
  obj <- readRDS(f)
  if (!inherits(obj, "GenomicRanges")) {
    warning("Skipping non-GRanges file: ", f)
    next
  }
  
  signal_orig <- extract_signal(obj)
  summary_orig <- summarize_locus_counts(signal_orig)
  
  for (rep_idx in seq_len(n_repeats)) {
    rep_tag <- sprintf("rep%02d", rep_idx)
    rep_seed <- seed + (i * 1000L) + rep_idx
    
    down <- downsample_signal_preserving(
      gr = obj,
      target_frac = target_fraction,
      min_keep = min_per_signal,
      max_keep = max_per_signal,
      seed_val = rep_seed,
      schedule_name = schedule
    )
    
    signal_down <- extract_signal(down$gr)
    summary_down <- summarize_locus_counts(signal_down)
    
    # Sanity check: enforce that downsampled signal sizes obey max_per_signal cap.
    max_down_per_signal <- max(as.integer(table(signal_down)))
    if (!is.infinite(max_per_signal) && max_down_per_signal > max_per_signal) {
      warning(
        "Downsampled signal size exceeded max_per_signal for ", trait_tag,
        " (", rep_tag, "): observed ", max_down_per_signal, " > ", max_per_signal
      )
    }
    
    n_signals_retained <- length(unique(signal_down))
    n_signals_orig <- length(unique(signal_orig))
    signal_retention <- n_signals_retained / n_signals_orig
    snp_retention <- length(signal_down) / length(signal_orig)
    
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      trait = trait_tag,
      input_file = f,
      replicate = rep_tag,
      seed_used = rep_seed,
      schedule = schedule,
      target_fraction = target_fraction,
      min_per_signal = min_per_signal,
      max_per_signal = ifelse(is.infinite(max_per_signal), NA_real_, max_per_signal),
      original_n_snps = summary_orig$n_snps,
      original_n_signals = summary_orig$n_signals,
      downsampled_n_snps = summary_down$n_snps,
      downsampled_n_signals = summary_down$n_signals,
      snp_retention = snp_retention,
      signal_retention = signal_retention,
      original_snps_per_signal_mean = summary_orig$snps_per_signal_mean,
      downsampled_snps_per_signal_mean = summary_down$snps_per_signal_mean,
      original_snps_per_signal_q95 = summary_orig$snps_per_signal_q95,
      downsampled_snps_per_signal_q95 = summary_down$snps_per_signal_q95,
      stringsAsFactors = FALSE
    )
    
    dist_requested <- build_dist_requested(signal_orig, signal_down, trait_tag, rep_tag, schedule)
    distribution_rows[[length(distribution_rows) + 1L]] <- dist_requested
    
    locus_count_orig <- as.data.frame(table(signal_orig), stringsAsFactors = FALSE)
    colnames(locus_count_orig) <- c("signal", "n_snps")
    locus_count_orig$trait <- trait_tag
    locus_count_orig$set <- "original"
    locus_count_orig$replicate <- rep_tag
    locus_count_orig$schedule <- schedule
    
    locus_count_down <- as.data.frame(table(signal_down), stringsAsFactors = FALSE)
    colnames(locus_count_down) <- c("signal", "n_snps")
    locus_count_down$trait <- trait_tag
    locus_count_down$set <- "downsampled"
    locus_count_down$replicate <- rep_tag
    locus_count_down$schedule <- schedule
    
    write.table(
      rbind(locus_count_orig, locus_count_down),
      file = file.path(trait_out_dir, paste0(trait_tag, "_", schedule, "_", rep_tag, "_locus_signal_counts.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    write.table(
      dist_requested,
      file = file.path(trait_out_dir, paste0(trait_tag, "_", schedule, "_", rep_tag, "_signal_size_bins_requested.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    plot_signal_size_bins_dodged(
      dist_requested = dist_requested,
      out_pdf = file.path(trait_out_dir, paste0(trait_tag, "_", schedule, "_", rep_tag, "_signal_size_bins_requested.pdf")),
      trait_tag = paste0(trait_tag, " | ", schedule, " | ", rep_tag)
    )
    
    saveRDS(
      down$gr,
      file = file.path(trait_out_dir, paste0(trait_tag, "_", schedule, "_", rep_tag, "_downsampled_signal_preserving.gr.rds"))
    )
    
    frac_tbl <- data.frame(
      trait = trait_tag,
      signal = names(down$n_by_signal),
      n_snps_original = as.integer(down$n_by_signal),
      snp_bin = as.character(snp_bin_requested(as.integer(down$n_by_signal))),
      keep_fraction_used = as.numeric(down$frac_by_signal),
      keep_n = as.integer(down$keep_by_signal),
      replicate = rep_tag,
      schedule = schedule,
      stringsAsFactors = FALSE
    )
    write.table(
      frac_tbl,
      file = file.path(trait_out_dir, paste0(trait_tag, "_", schedule, "_", rep_tag, "_fraction_by_signal.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

if (!length(summary_rows)) {
  stop("No valid GRanges trait files were processed.")
}

summary_df <- do.call(rbind, summary_rows)
dist_df <- do.call(rbind, distribution_rows)

trait_tags <- unique(summary_df$trait)
run_tag <- if (length(trait_tags) == 1) trait_tags else paste0("all_traits_n", length(trait_tags))

write.table(
  summary_df,
  file = file.path(output_dir, paste0(run_tag, "_", schedule, "_trait_downsampling_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  dist_df,
  file = file.path(output_dir, paste0(run_tag, "_", schedule, "_trait_locus_size_distribution.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

schedule_export <- if (schedule == "global") {
  data.frame(schedule = "global", snp_bin = SNP_BIN_LEVELS, keep_fraction = target_fraction, stringsAsFactors = FALSE)
} else {
  SCHEDULE_TABLE[SCHEDULE_TABLE$schedule == schedule, ]
}
write.table(
  schedule_export,
  file = file.path(output_dir, paste0(run_tag, "_", schedule, "_schedule_table.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Done.")
message("  Summary: ", file.path(output_dir, paste0(run_tag, "_", schedule, "_trait_downsampling_summary.tsv")))
message("  Distribution: ", file.path(output_dir, paste0(run_tag, "_", schedule, "_trait_locus_size_distribution.tsv")))
message("  Schedule table: ", file.path(output_dir, paste0(run_tag, "_", schedule, "_schedule_table.tsv")))
message("  Per-trait subdirectories + outputs written to: ", output_dir)
