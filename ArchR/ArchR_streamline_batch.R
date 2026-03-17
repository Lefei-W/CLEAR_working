# Self-contained ArchR batch pipeline for multiple tissues and samples.

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(GenomicRanges)
  library(ggplot2)
  library(qpdf)
  library(cowplot)
  library(patchwork)
  library(ComplexHeatmap)
  library(dplyr)
})

addArchRGenome("hg38")
addArchRThreads(24) 
source('~/plottheme.R')

DEFAULT_OUTPUT_BASE <- "/working/lab_jonathb/lefeiW/projects/CLEAR_data/ArchR_batch_processed"

print_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript ArchR_streamline_batch.R <input_tsv> [output_base] [genome] [threads] [minTSS] [minFrags] [tissue_filter] [macs2_path]\n\n",
    "Arguments:\n",
    "  input_tsv   Delimited file with columns: tissue, sample, frag_path, [workdir]\n",
    "  output_base Base output directory (one subfolder per tissue).\n",
    "              Default: /working/lab_jonathb/lefeiW/projects/CLEAR_data/ArchR_batch_processed\n",
    "  genome      ArchR genome, default: hg38\n",
    "  threads     Number of threads, default: 24\n",
    "  minTSS      Minimum TSS enrichment, default: 4\n",
    "  minFrags    Minimum fragments, default: 3200\n",
    "  tissue_filter Optional tissue name to run only one tissue\n",
    "  macs2_path  Path to MACS2 executable for peak calling\n\n",
    "Example:\n",
    "  Rscript ArchR_streamline_batch.R inputs.tsv /path/to/output hg38 24 4 3200 breast_regner_2025_healthy /software/MACS2/MACS2-20241021/Python-3.7.7-venv/bin/macs2\n"
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
  print_usage()
  quit(status = 0)
}

if (length(args) < 1) {
  print_usage()
  stop("Missing required argument: input_tsv")
}

input_tsv <- args[1]

get_arg_or_default <- function(i, default) {
  if (length(args) >= i && nzchar(args[i])) {
    return(args[i])
  }
  default
}

output_base <- get_arg_or_default(2, DEFAULT_OUTPUT_BASE)
genome <- get_arg_or_default(3, "hg38")
threads <- as.integer(get_arg_or_default(4, "24"))
minTSS <- as.numeric(get_arg_or_default(5, "4"))
minFrags <- as.integer(get_arg_or_default(6, "3200"))
tissue_filter <- get_arg_or_default(7, "")
path_to_macs2 <- get_arg_or_default(8, "/software/MACS2/MACS2-20241021/Python-3.7.7-venv/bin/macs2")

if (!file.exists(input_tsv)) {
  stop("Input TSV not found: ", input_tsv)
}

dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

ArchR::addArchRGenome(genome)
ArchR::addArchRThreads(threads)

inputs <- read.table(input_tsv, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("tissue", "sample", "frag_path")
missing_cols <- setdiff(required_cols, colnames(inputs))
if (length(missing_cols) > 0) {
  stop("Input TSV is missing columns: ", paste(missing_cols, collapse = ", "))
}

if (!"workdir" %in% colnames(inputs)) {
  inputs$workdir <- ""
}

if (nzchar(tissue_filter)) {
  inputs <- inputs[inputs$tissue == tissue_filter, , drop = FALSE]
  if (nrow(inputs) == 0) {
    stop("No rows found for tissue_filter: ", tissue_filter)
  }
}

resolve_path <- function(path, base_dir) {
  if (!nzchar(base_dir)) {
    return(normalizePath(path, mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, path), mustWork = FALSE)
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

run_tissue <- function(tissue_name, tissue_df) {
  tissue_dir <- file.path(output_base, tissue_name)
  dir.create(tissue_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tissue_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

  # Keep all outputs together for the tissue
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tissue_dir)

  tissue_df$frag_path <- mapply(
    FUN = resolve_path,
    path = tissue_df$frag_path,
    base_dir = ifelse(is.na(tissue_df$workdir), "", as.character(tissue_df$workdir)),
    USE.NAMES = FALSE
  )
  frag_paths <- setNames(tissue_df$frag_path, tissue_df$sample)
  sample_names <- tissue_df$sample

  arrow_files <- createArrowFiles(
    inputFiles = frag_paths,
    sampleNames = sample_names,
    minTSS = minTSS,
    minFrags = minFrags,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
  )

  project_dir <- paste0("ArchR_", sanitize_name(tissue_name))
  proj <- ArchRProject(
    ArrowFiles = arrow_files,
    outputDirectory = project_dir,
    copyArrows = TRUE
  )

  proj <- addDoubletScores(
    input = proj,
    k = 10,
    knnMethod = "UMAP",
    LSIMethod = 1
  )

  proj <- filterDoublets(proj)

  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI_ArchR",
    iterations = 2,
    clusterParams = list(
      resolution = c(0.2),
      sampleCells = 10000,
      n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30
  )

  proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI_ArchR",
    name = "Harmony_ArchR",
    groupBy = "Sample",
    force = TRUE
  )

  proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI_ArchR",
    method = "Seurat",
    name = "Clusters_ArchR",
    resolution = 0.8
  )

  proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "IterativeLSI_ArchR",
    name = "UMAP_ArchR",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
  )

  proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "Harmony_ArchR",
    name = "UMAPHarmony_ArchR",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
  )

  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony_ArchR")
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters_ArchR", embedding = "UMAPHarmony_ArchR")
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP_ArchR")
  p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters_ArchR", embedding = "UMAP_ArchR")

  pdf(file.path("plots", "UMAP_default.pdf"), width = 20, height = 8)
  print(p1 | p2 | p3 | p4)
  dev.off()

  markers_gs <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters_ArchR",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )

  heatmap_gs <- plotMarkerHeatmap(
    seMarker = markers_gs,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    transpose = TRUE
  )

  p_heat <- ComplexHeatmap::draw(
    heatmap_gs,
    heatmap_legend_side = "bot",
    annotation_legend_side = "bot"
  )

  pdf(file.path("plots", "geneMarkers_heatmap.pdf"), width = 12, height = 6)
  print(p_heat)
  dev.off()

  proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "Clusters_ArchR"
  )

  proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "Clusters_ArchR",
    force = TRUE,
    pathToMacs2 = path_to_macs2,
    maxPeaks = 300000,
    cutOff = 0.01
  )

  proj <- addPeakMatrix(proj)

  groupse <- getGroupSE(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters_ArchR"
  )
  saveRDS(groupse, file.path(tissue_dir, "peakcalls_Clusters_ArchR_groupSE_PeakMatrix.rds"))

  saveRDS(proj, file.path(tissue_dir, paste0(tissue_name, "_ArchRProject.rds")))
}

tissues <- unique(inputs$tissue)
for (tissue_name in tissues) {
  tissue_df <- inputs[inputs$tissue == tissue_name, , drop = FALSE]
  run_tissue(tissue_name, tissue_df)
}


