# GWAS input for CLEAR

## 1. What CLEAR expects

A single **`GRanges`** object (saved as `.rds`) on the **same genome
build as your snATAC SE** (usually hg38).

Minimum:

| Field | Required | Notes |
|---|---|---|
| `seqnames`, `start`, `end` | yes | SNPs: `start == end`. |
| `names` (or `names(snps)`) | recommended | rsID or `chr:pos:ref:alt`. |
| `mcols(snps)$signal` | optional | Locus / credible-set ID. Only needed for locus-grouped plots. |

CLEAR only reads `names` and `signal`; everything else in `mcols`
is ignored but preserved.

```r
library(GenomicRanges)
snps <- GRanges(seqnames = "chr10", ranges = IRanges(123456, 123456),
                names = "rs123", signal = "locus1")
saveRDS(snps, "my_trait_GR.rds")
```

---

## 2. If you only have a list of variants

No locus / fine-mapping info needed — drop the `signal` column. Locus
grouped plots will collapse to one peak per locus, everything else
runs unchanged.

```r
library(readr); library(GenomicRanges)

bed <- read_tsv("variants.tsv")        # columns: chr, pos, rsid
snps <- GRanges(
  seqnames = bed$chr,
  ranges   = IRanges(start = bed$pos, end = bed$pos),
  names    = bed$rsid
)
saveRDS(snps, "my_variants_GR.rds")
```

---

## 3. From a trait of interest -> Open Targets Genetics API

If you just have a disease/trait name, pull credible sets directly
from Open Targets Genetics (GraphQL):

```r
# OpenTargets API 
```

Open Targets coordinates are GRCh38 (hg38). Replace `studyId` with the
GWAS Catalog accession for your trait.

---

## 4. From a list of index/lead SNPs -> PLINK LD expansion

If you have just lead SNPs and want all variants in LD r² >= 0.8 with
locus labels:

```bash
# leads.tsv: one rsID per line; ref panel = 1000G EUR on hg38
plink \
  --bfile 1000G_EUR_hg38 \
  --r2 \
  --ld-snp-list leads.tsv \
  --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0.8 \
  --out leads_ld
# -> leads_ld.ld with columns: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
```

```r
library(readr); library(GenomicRanges)

ld <- read_table("leads_ld.ld")      # one row per (lead, tag) pair

snps <- GRanges(
  seqnames = paste0("chr", ld$CHR_B),
  ranges   = IRanges(start = ld$BP_B, end = ld$BP_B),
  names    = ld$SNP_B,
  signal   = ld$SNP_A                # lead SNP id = locus label
)
saveRDS(snps, "leads_ld08_GR.rds")
```

Confirm with `seqlevelsStyle(snps)` and `genome(snps) <- "hg38"` before
saving if your reference uses different chromosome naming.
