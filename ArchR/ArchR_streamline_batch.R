# Self-contained ArchR batch pipeline for multiple tissues and samples.

suppressPackageStartupMessages({
  library(ArchR), library(BSgenome.Hsapiens.UCSC.hg38), library(Biostrings), library(GenomicRanges),
  library(ggplot2), library(qpdf), library(cowplot), library(patchwork), library(ComplexHeatmap),library(dplyr)
})

addArchRGenome("hg38")
addArchRThreads(24) 
source('~/plottheme.R')

print_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript ArchR_streamline_batch.R <input_tsv> <output_base> [genome] [threads] [minTSS] [minFrags]\n\n",
    "Arguments:\n",
    "  input_tsv   Tab-delimited file with columns: tissue, sample, frag_path, workdir\n",
    "  output_base Base output directory (one subfolder per tissue)\n",
    "  genome      ArchR genome, default: hg38\n",
    "  threads     Number of threads, default: 24\n",
    "  minTSS      Minimum TSS enrichment, default: 4\n",
    "  minFrags    Minimum fragments, default: 3200\n\n",
    "Example:\n",
    "  Rscript ArchR_streamline_batch.R inputs.tsv /path/to/output hg38 24 4 3200\n"
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
  print_usage()
  quit(status = 0)
}

if (length(args) < 2) {
  print_usage()
  stop("Missing required arguments: input_tsv and output_base")
}

input_tsv <- args[1]
output_base <- args[2]

get_arg_or_default <- function(i, default) {
  if (length(args) >= i && nzchar(args[i])) {
    return(args[i])
  }
  default
}

genome <- get_arg_or_default(3, "hg38")
threads <- as.integer(get_arg_or_default(4, "24"))
minTSS <- as.numeric(get_arg_or_default(5, "4"))
minFrags <- as.integer(get_arg_or_default(6, "3200"))

if (!file.exists(input_tsv)) {
  stop("Input TSV not found: ", input_tsv)
}

ArchR::addArchRGenome(genome)
ArchR::addArchRThreads(threads)

inputs <- read.table(input_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
required_cols <- c("tissue", "sample", "frag_path", "workdir")
missing_cols <- setdiff(required_cols, colnames(inputs))
if (length(missing_cols) > 0) {
  stop("Input TSV is missing columns: ", paste(missing_cols, collapse = ", "))
}

resolve_path <- function(path, base_dir) {
  is_abs_unix <- grepl("^/", path)
  is_abs_win <- grepl("^[A-Za-z]:\\\\", path)
  is_abs_unc <- grepl("^\\\\\\\\", path)
  if (is_abs_unix || is_abs_win || is_abs_unc) {
    return(normalizePath(path, mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, path), mustWork = FALSE)
}

run_tissue <- function(tissue_name, tissue_df) {
  workdirs <- unique(tissue_df$workdir)
  if (length(workdirs) != 1) {
    stop("Multiple workdir values found for tissue: ", tissue_name)
  }
  workdir <- workdirs[1]

  tissue_dir <- file.path(output_base, tissue_name)
  dir.create(tissue_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tissue_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

  # Keep all outputs together for the tissue
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tissue_dir)

  tissue_df$frag_path <- vapply(
    tissue_df$frag_path,
    resolve_path,
    character(1),
    base_dir = workdir
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

  proj <- ArchRProject(
    ArrowFiles = arrow_files,
    outputDirectory = file.path(tissue_dir, "ArchRProject"),
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

  saveRDS(proj, file.path(tissue_dir, paste0(tissue_name, "_ArchRProject.rds")))
}

tissues <- unique(inputs$tissue)
for (tissue_name in tissues) {
  tissue_df <- inputs[inputs$tissue == tissue_name, , drop = FALSE]
  run_tissue(tissue_name, tissue_df)
}
