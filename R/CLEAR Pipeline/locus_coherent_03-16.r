

# 2026-03-15, sunday Jonno's run on GWAS signal coherence across peaks within a locus (signals with more than 1 GWAS peak)

# Peak-level dominance - geometric specificity within a regulatory element
# Locus-level mean dominance - consistent cell-type contribution across candidate elements
# Locus-level rank concordance - whether peaks agree on relative ordering
# Locus-level lineage structure

library(dplyr)
library(tidyr)
library(ggplot2)
library(corpcor)
library(ComplexHeatmap)
library(gridExtra)
library(SummarizedExperiment)


# function for L2 normalisation of initial matrix (input matrix e.g. ArchR 'PeakMatrix')
L2_mat <- function(mat){
  
  row_norms <- sqrt(rowSums(mat^2))
  row_norms[row_norms == 0] <- 1
  mat_l2 <- mat / row_norms
  return(mat_l2)
}

# define funciotn to classify "dominant" cells [ie. is the top most specific cell type twice the signal as the second]
compute_dominance <- function(mat, ratio_thresh = 2) {
  
  message('identifying dominant cell types...')
  cell_names <- colnames(mat)
  
  result <- t(apply(mat, 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    max_val <- x[ord[1]]
    second_val <- x[ord[2]]
    ratio <- max_val / second_val
    dominant_ct <- cell_names[ord[1]]
    c(dominant_ct = dominant_ct,
      ratio = ratio,
      is_dominant = ratio > ratio_thresh)
  }))
  
  result <- as.data.frame(result)
  result$ratio <- as.numeric(result$V2)
  result$is_dominant <- result$V3 == "TRUE"
  result$peak <- rownames(mat)
  return(result)
  
}

# locus level, with permutation testing to control for configuration wrt promoters
compute_locus_rank_concordance_vec <- function(
    gwas_mat,
    full_mat,        # L2 mat
    peak_to_locus,   # named vector: peak -> locus ID
    joint_bin_full,  # named vector for full peaks (TSS + density combined)
    n_perm = 100
) {
  
  # -----------------------
  # Precompute rank matrices
  # -----------------------
  rank_gwas <- t(apply(gwas_mat, 1, rank))
  rank_full <- t(apply(full_mat, 1, rank))
  
  # Split GWAS peaks by locus
  locus_split <- split(rownames(gwas_mat),
                       peak_to_locus[rownames(gwas_mat)])
  
  # Precompute full peak bins (joint bins)
  full_bins <- split(names(joint_bin_full), joint_bin_full)
  
  results <- lapply(names(locus_split), function(locus) {
    
    peaks <- locus_split[[locus]]
    k <- length(peaks)
    
    if (k < 2) {
      return(data.frame(
        locus = locus,
        n_peaks = k,
        concordance = NA,
        mean = NA,
        sd = NA,
        z_score = NA,
        p_value = NA,
        dominant_cell_locus = NA,
        dominant_score_locus = NA,
        dominant_lineage = NA,
        lineage_fraction = NA))
    }
    
    # -------------------
    # Observed concordance
    # -------------------
    obs_mat <- rank_gwas[peaks, , drop = FALSE]
    obs_cor <- cor(t(obs_mat), method = "spearman")
    obs_conc <- mean(obs_cor[upper.tri(obs_cor)])
    
    # -------------------
    # Dominant cell
    # -------------------
    mean_profile <- colMeans(gwas_mat[peaks, , drop = FALSE])
    dominant_cell_locus <- names(which.max(mean_profile))
    dominant_score_locus <- max(mean_profile)
    
    # -------------------
    # Lineage aggregation
    # -------------------
    lineage_profile <- tapply(mean_profile, cell_lineage, sum)
    dominant_lineage <- names(which.max(lineage_profile))
    lineage_fraction <- max(lineage_profile) / sum(lineage_profile)
    
    # -------------------
    # Joint bin composition (subset from full)
    # -------------------
    locus_joint_bins <- joint_bin_full[peaks]
    locus_bins <- table(locus_joint_bins)
    
    # # -------------------
    # # OLD !!! Permutation null (matched by joint-bin composition)
    # # -------------------
    # perm_vals <- replicate(n_perm, {
    #
    #   sampled_peaks <- unlist(lapply(names(locus_bins), function(bin) {
    #     sample(full_bins[[bin]], locus_bins[[bin]])
    #   }))
    #
    #   perm_mat <- rank_full[sampled_peaks, , drop = FALSE]
    #   perm_cor <- cor(t(perm_mat), method = "spearman")
    #   mean(perm_cor[upper.tri(perm_cor)])
    # })
    
    # -------------------
    # Permutation null (matched by joint-bin composition)
    # -------------------
    
    # Preallocate matrix of sampled peaks
    sampled_matrix <- matrix(NA_character_, nrow = k, ncol = n_perm)
    
    # Fill it bin-by-bin (vectorised across permutations)
    row_index <- 1
    for (bin in names(locus_bins)) {
      
      n_bin <- locus_bins[[bin]]
      
      draws <- replicate(
        n_perm,
        sample(full_bins[[bin]], n_bin),
        simplify = "matrix")
      
      sampled_matrix[row_index:(row_index + n_bin - 1), ] <- draws
      row_index <- row_index + n_bin
    }
    
    # Compute concordance for all permutations
    perm_vals <- apply(sampled_matrix, 2, function(sampled_peaks) {
      perm_mat <- rank_full[sampled_peaks, , drop = FALSE]
      perm_cor <- cor(t(perm_mat), method = "spearman")
      mean(perm_cor[upper.tri(perm_cor)])
    })
    
    perm_mean <- mean(perm_vals)
    perm_sd   <- sd(perm_vals)
    
    data.frame(
      locus = locus,
      n_peaks = k,
      concordance = obs_conc,
      mean = perm_mean,
      sd = perm_sd,
      z_score = (obs_conc - perm_mean) / perm_sd,
      p_value = mean(perm_vals >= obs_conc),
      dominant_cell_locus = dominant_cell_locus,
      dominant_score_locus = dominant_score_locus,
      dominant_lineage = dominant_lineage,
      lineage_fraction = lineage_fraction
    )
  })
  do.call(rbind, results)
}

# function to get TSS positions from GTF file [gtf is GENCODE.vXX]
get_tss <- function(gtf){
  
  genes <- gtf[gtf$type == "gene"]
  genes <- genes[genes$gene_type %in% c("lncRNA","protein_coding")]
  tss_gr <- GenomicRanges::promoters(
    genes,
    upstream = 0,
    downstream = 1
  )
  S4Vectors::mcols(tss_gr) <- S4Vectors::mcols(tss_gr)[,c("gene_type","gene_name")]
  
  return(tss_gr)
}

# peak density (peak distribution is a potential confounder)
compute_peak_density <- function(peaks_gr, window = 50000) {
  
  counts <- countOverlaps(
    peaks_gr,
    resize(peaks_gr, width = window, fix = "center")
  )
  
  # subtract self
  counts - 1
}

# bin peaks by density
make_density_bins <- function(density_vec, nbins = 4) {
  
  bins <- cut(
    density_vec,
    breaks = unique(quantile(density_vec,
                             probs = seq(0, 1, length.out = nbins + 1),
                             na.rm = TRUE)),
    include.lowest = TRUE
  )
  
  bins
}

prepare_clear_se <- function(se, genome = "hg38") {
  rd <- as.data.frame(rowData(se))
  
  # Build GRanges from rowData (ArchR groupSE always has seqnames/start/end)
  if (!all(c("seqnames", "start", "end") %in% names(rd))) {
    stop("rowData must contain seqnames, start, end columns (ArchR PeakMatrix format expected)") # Check what Signac does
  }
  
  gr <- GRanges(
    seqnames = rd$seqnames,
    ranges   = IRanges::IRanges(start = as.integer(rd$start), end = as.integer(rd$end))
  )
  
  # Keep rowData columns 
  extra_cols <- setdiff(names(rd), c("seqnames", "start", "end", "width", "strand"))
  if (length(extra_cols)) mcols(gr) <- rd[, extra_cols, drop = FALSE]
  
  # Rebuild as RSE
  se <- SummarizedExperiment(
    assays    = as.list(assays(se)),
    rowRanges = gr,
    colData   = colData(se)
  )
  
  # Rownames as genomic coords; peakid in rowData
  rownames(se) <- paste0(as.character(seqnames(gr)), ":", start(gr), "-", end(gr))
  rowData(se)$peakid <- rownames(se)
  
  try(seqlevelsStyle(rowRanges(se)) <- "UCSC", silent = TRUE)
  try(GenomeInfoDb::genome(rowRanges(se)) <- genome, silent = TRUE)
  
  se
}

# GTF file
gtf_file <- '~/R_projects/nearest_gene/data/gencode.v46.basic.annotation.gtf.gz'
gtf_file <- '../CLEAR_data/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz'
gtf <- rtracklayer::import(gtf_file)

# CCVs
snps <- readRDS("/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/BCAC_FM_GR.rds")

# ATAC SE
se <- readRDS("/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/breast_330017peak.rds")
se <- prepare_clear_se(se)
# L2 normalise ArchR peak matrix
mat_l2 <- L2_mat(assay(se,1))
assay(se, 'l2_norm') <- mat_l2

# intersect SNPs
gwas_mat <- assay(unique(se[S4Vectors::subjectHits(GenomicRanges::findOverlaps(snps, se))]), 'l2_norm')

## plot correlation
assay(se,1) %>% cor() %>% heatmap()
assay(se,2) %>% cor() %>% heatmap()
gwas_mat %>% cor() %>% heatmap()

# summarise patterns
high_thresh <- 1/sqrt(2)   #1 (or 2-cell) theoretical boundary (1 or 2 cells must be > 1/sqrt(3); 3 can be = 1/sqrt(3))
low_thresh  <- 1/sqrt(ncol(mat_l2))

# full L2norm Mat
specificity_summary_FULL <- data.frame(
  peak = rownames(mat_l2),
  n_high = rowSums(mat_l2 > high_thresh),
  n_mid  = rowSums(mat_l2 > low_thresh & mat_l2 <= high_thresh),
  n_low  = rowSums(mat_l2 <= low_thresh))

# GWAS L2norm Mat
specificity_summary_GWAS <- data.frame(
  peak = rownames(gwas_mat),
  n_high = rowSums(gwas_mat > high_thresh),
  n_mid  = rowSums(gwas_mat > low_thresh & gwas_mat <= high_thresh),
  n_low  = rowSums(gwas_mat <= low_thresh))

# this is a sanity check that the GWAS peaks include high spec
prop_df <- as.data.frame(prop.table(table(specificity_summary_GWAS$n_high)))
colnames(prop_df) <- c("n_high", "proportion")

# single stacked bar
ggplot(prop_df, aes(x = "GWAS peaks", y = proportion, fill = n_high)) +
  geom_bar(stat = "identity", width = 0.6) + theme_bw() + ggsci::scale_fill_bmj() +
  labs(title = "Distribution of High-Specific Cell Counts per Peak",
       x = "", y = "Proportion",
       fill = "# Highly specific\ncell types") + coord_flip()

# ID dominant cell types
dominance_all <- compute_dominance(mat_l2, ratio_thresh = 2)
dominance_gwas <- compute_dominance(gwas_mat, ratio_thresh = 2)

all_counts <- dominance_all %>% filter(is_dominant) %>% dplyr::count(V1)
gwas_counts <- dominance_gwas %>% filter(is_dominant) %>% dplyr::count(V1)

counts <- left_join(all_counts, gwas_counts, by = "V1")
colnames(counts) <- c('cell','genome-wide','GWAS')
counts <- tidyr::gather(counts, set, n, -cell)
counts <- counts %>% group_by(set) %>% mutate(prop = n/sum(n, na.rm = T))

# plot
ggplot(counts, aes(x = reorder(cell, ifelse(set == "GWAS", prop, NA), na.rm = TRUE), y = prop, fill = set)) +
  geom_bar(stat = "identity", position = 'dodge') +
  theme_bw() + scale_fill_manual(values = c("GWAS"="red","genome-wide"="gray")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell type", y = "Proportion dominant peaks (ratio > 2)",
       title = "Cell-Type Dominant peaks (L2mat)") + coord_flip()

# Annotate summary with SNP info
# Make GRanges directly from peak names
peak_gr <- GRanges(rownames(specificity_summary_GWAS))

# Find overlaps with GWAS SNPs
hits <- findOverlaps(peak_gr, snps)

# Annotate peaks with overlapping SNPs
snp_list <- tapply(names(snps)[subjectHits(hits)], queryHits(hits), paste, collapse=",")
locus_vec <- tapply(mcols(snps)$signal[subjectHits(hits)],
                    queryHits(hits),
                    function(x) paste(unique(x), collapse = ", "))

specificity_summary_GWAS$gwas_snps[as.numeric(names(snp_list))] <- snp_list
specificity_summary_GWAS$gwas_locus[as.numeric(names(locus_vec))] <- locus_vec


######################
# Locus level analysis
######################

# get TSS positions
tss_gr <- get_tss(gtf)
peaks_gr <- rowRanges(se)
nearest_hits <- distanceToNearest(peaks_gr, tss_gr)
dist2tss <- mcols(nearest_hits)$distance

# attach to peaks
mcols(peaks_gr)$dist2tss <- dist2tss

# categorise
tss_bin <- cut(
  dist2tss,
  breaks = c(0, 1000, 5000, 50000, 200000, Inf),
  labels = c("core_prom", "prox_prom", "near_reg", "distal", "long_range"),
  include.lowest = TRUE)

names(tss_bin) <- rownames(mat_l2)

# named vectors for input - TSS
tss_bin_full <- tss_bin
tss_bin_gwas <- tss_bin[rownames(gwas_mat)]
peak_to_locus <- specificity_summary_GWAS$gwas_locus
names(peak_to_locus) <- specificity_summary_GWAS$peak

# named vectors for input - peak density
x <- compute_peak_density(rowRanges(se))
density_bin_full <- make_density_bins(x)
names(density_bin_full) <- rownames(mat_l2)

# combine tss and peak density for input into lucus permutation
joint_bin_full <- interaction(tss_bin_full, density_bin_full, drop = TRUE)
joint_bin_full <- as.character(joint_bin_full)
names(joint_bin_full) <- names(tss_bin_full)

# Lineage aggregation
# cluster cell types
cell_cor <- cor(mat_l2, method = "pearson")   # or spearman
hc <- hclust(as.dist(1 - cell_cor), method = "average")
plot(hc)
# choose k lineages
k <- 4
cell_lineage <- cutree(hc, k = k)

# Analyse locus
locus_results <- compute_locus_rank_concordance_vec(
  gwas_mat = gwas_mat,
  full_mat = mat_l2,
  peak_to_locus = peak_to_locus,
  joint_bin_full = joint_bin_full,
  n_perm = 1000  # real analysis should use 1000
)

# summarise
locus_results$significant <- locus_results$p_value < 0.05

# plot
ggplot(locus_results %>% filter(n_peaks>1), aes(x = concordance, y = -log10(p_value))) +
  geom_point(aes(size = n_peaks, color = significant), alpha = 0.8) +
  theme_bw() + ggsci::scale_color_igv() +
  labs(x = "Mean pairwise rank concordance",
       y = "-log10(p-value)",
       title = "Locus-level regulatory concordance")

# categorise
locus_results <- locus_results %>%
  mutate(category = case_when(
    n_peaks == 1 ~ "Single_peak",
    z_score > 2 & p_value < 0.05 & dominant_score_locus > 0.5 ~ "Highly_coherent_specific",
    z_score > 2 & p_value < 0.05 ~ "Coherent_moderate",
    z_score < 1 ~ "Incoherent",
    TRUE ~ "Intermediate"))

# different highlights
ggplot(locus_results %>% filter(n_peaks>1), aes(x = concordance, y = -log10(p_value))) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_point(aes(size = n_peaks, color = category), alpha = 0.8) +
  theme_bw() + ggsci::scale_color_igv() +
  labs(x = "Mean pairwise rank concordance",
       y = "-log10(p-value)",
       title = "Locus-level regulatory concordance")

# counts
locus_results %>% filter(n_peaks>1) %>%
  select(locus, n_peaks, dominant_cell_locus, dominant_lineage, category) %>%
  arrange(category)


#  | Column                | Interpretation                        |
#  | --------------------- | ------------------------------------- |
#  | observed              | coherence among GWAS peaks within a locus`
#  | dominant_cell_locus ` | Which cell type most peaks prioritise |
#  | dominant_score_locus  | Max score                 |

