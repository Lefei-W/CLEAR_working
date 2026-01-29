# All processed inputs to be listed here (define name, path, assays)
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

DROP_CTS_DEFAULT <- c("Basal2", "Intermediate") # noisy in-house clusters to drop and test effects on enrichment analyses

CHAIN_FILE <- "/working/lab_jonathb/lefeiW/projects/sc_eQTL/results/Cleaned_BC_subtype_GWAS/hg38ToHg19.over.chain"
