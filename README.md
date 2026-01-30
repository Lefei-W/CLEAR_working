# CLEAR – Cell-type-informed Locus-specific Enrichment Analysis by Ranking

CLEAR is a framework for quantifying and ranking cell-type-specific regulatory contributions at individual GWAS risk loci by integrating fine-mapped variants with single-cell chromatin accessibility data. It provides locus-level, quantitative measures of cell-type involvement beyond global enrichment or single “risk-driving” cell type identification.

---

## Table of Contents

- [Overview](#overview)
- [Method Abstract](#method-abstract)
- [Input Data](#input-data)
- [ArchR Peak Signal Matrix Construction](#archr-peak-signal-matrix-construction)
- [Specificity Metrics](#specificity-metrics)
- [Enrichment Framework](#enrichment-framework)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Outputs](#outputs)

---

## Overview

Genome-wide association studies (GWAS) have been highly successful in identifying germline risk loci for complex traits, yet the majority of associated variants reside in non-coding regions and are believed to act through cell-type-specific regulatory mechanisms. Single-cell ATAC-seq (scATAC-seq) enables the construction of cell-type-specific (CTS) regulatory landscapes, providing an opportunity to map risk variants to regulatory elements in a cell-type-resolved manner.

Existing integrative methods combining GWAS with single-cell genomics, epigenomics, or spatial transcriptomics have primarily focused on identifying key risk-driving cell types at a genome-wide level, often recovering well-established cells of origin. However, these approaches generally lack quantitative measurements of cell-type contribution at individual loci and are not designed to highlight unexpected or secondary cell types of action. Such non-canonical cell-type effects, including immune or stromal contributions in solid tumors, may represent important but underexplored biological mechanisms and potential therapeutic targets.

CLEAR addresses this gap by providing a locus-level, quantitative framework that ranks cell-type specificity of regulatory elements overlapping fine-mapped GWAS variants. By comparing specificity ranks of variant-overlapping peaks against appropriate null models, CLEAR estimates how strongly each cell type contributes to regulatory activity at a given locus. The framework is modular, scalable, and applicable to diverse complex traits and single-cell regulatory datasets.

---

## Method Abstract

![CLEAR graphical abstract](images/CLEAR_AGTA_lightning-1.png)

---

## Input Data

<!-- TODO: GWAS fine-mapped variants (~1500 traits from OpenTargets), scATAC-seq peak × cell-type matrices (ArchR processed) -->

- GWAS fine-mapped variants from [OpenTargets Genetics](https://platform.opentargets.org/)
- scATAC-seq peak × cell-type matrices generated using [ArchR](https://www.archrproject.com/)

### Data availability

<!-- TODO: subsection of the input data -->

---

## ArchR Peak Signal Matrix Construction

<!-- TODO: ArchR processing pipeline -->

- Summarised ArchR essential processing steps [ArchR processing](https://docs.google.com/document/d/12TJ8RZt97AsLWTHmqDecG8jwKzTnrnOjEpVfefXgW1M/edit?tab=t.0)

---

## Specificity Metrics

<!-- TODO: L2 norm, z-score, Tau, composite metrics, etc. -->

---

## Enrichment Framework

![CLEAR graphical workflow](images/CLEAR_flowchart-1.png)

---

## Usage

---

## Project Structure

<!-- TODO: to be adapted from new CLEAR folder -->
<!-- TODO: scripts to be modularised and stored under /R -->

---

## Outputs

---
