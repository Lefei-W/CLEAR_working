# CLEAR Locus Coherence Methods

## Overview

The CLEAR pipeline measures whether GWAS-associated regulatory elements (peaks) within a genetic locus show consistent cell-type specificity patterns. Two complementary locus-level metrics are computed, each capturing a different aspect of regulatory coherence.

---

## 1. Rank Concordance

### Question
> Do the peaks within a locus agree on the **full ordering** of cell types by specificity?

### Calculation

Given a locus with $k$ GWAS-overlapping peaks and $C$ cell types:

**Step 1 — Row-wise ranking.** For each peak $i$, replace its L2-normalised specificity vector with within-peak ranks across the $C$ cell types:

$$r_{ic} = \text{rank}(s_{ic}) \quad \text{among } \{s_{i1}, s_{i2}, \ldots, s_{iC}\}$$

where $s_{ic}$ is the L2-normalised specificity of peak $i$ in cell type $c$. This converts continuous specificity scores to ordinal positions (1 = lowest, $C$ = highest) **within each peak**.

**Step 2 — Pairwise Spearman correlation.** For every pair of peaks $(i, j)$ in the locus, compute the Spearman rank correlation across the $C$ cell types:

$$\rho_{ij} = \text{cor}_{\text{Spearman}}(\mathbf{r}_i, \mathbf{r}_j)$$

This measures whether peak $i$ and peak $j$ rank cell types in the same order. Note: since the data is already ranked, this is equivalent to a Pearson correlation on the ranks.

**Step 3 — Locus concordance score.** Average all pairwise correlations:

$$\text{Concordance} = \frac{2}{k(k-1)} \sum_{i < j} \rho_{ij}$$

### Score interpretation

| Value | Meaning |
|-------|---------|
| +1.0 | All peaks rank cell types identically |
| 0.0 | No systematic agreement (random) |
| −1.0 | Peaks rank cell types in opposite orders |

A locus can have high concordance even if no single cell type dominates — what matters is that peaks **agree on the relative ordering** (e.g., all peaks rank cell A > cell C > cell B).

### Dominant cell type annotation

The dominant cell type is identified from the **mean L2 specificity profile** across the locus peaks:

$$\bar{s}_c = \frac{1}{k} \sum_{i=1}^{k} s_{ic}$$

$$\text{dominant cell} = \arg\max_c \; \bar{s}_c$$

Lineage-level dominance sums $\bar{s}_c$ within each lineage group $g$:

$$L_g = \sum_{c \in g} \bar{s}_c, \qquad \text{lineage fraction} = \frac{\max_g L_g}{\sum_g L_g}$$

---

## 2. Consensus Concordance

### Question
> Do the peaks within a locus agree on **which single cell type (or lineage)** is dominant?

### Calculation

Uses `dominance_consensus()` — a weighted plurality vote across peaks.

**Step 1 — Group specificity per peak.** For each peak $i$, sum L2 specificity scores within each group $g$ (group = individual cell type at fine level, or lineage cluster at coarse level):

$$G_{ig} = \sum_{c \in g} s_{ic}$$

**Step 2 — Peak assignment and weighting.** Assign each peak to its top group:

$$g_i^* = \arg\max_g \; G_{ig}$$

Weight each peak's vote by both its specificity strength and its overall accessibility:

$$w_i = G_{ig_i^*} \times \bar{a}_i$$

where $\bar{a}_i$ = mean raw accessibility of peak $i$ across cell types (the L2 weight). This ensures that peaks which are both **specific** and **accessible** have the strongest influence.

**Step 3 — Plurality tally.** Sum weights by assigned group:

$$T_g = \sum_{i : g_i^* = g} w_i$$

**Step 4 — Consensus score.** The fraction of total weight captured by the winning group:

$$\text{Consensus} = \frac{\max_g \; T_g}{\sum_g T_g}$$

### Score interpretation

| Value | Meaning |
|-------|---------|
| 1.0 | All peaks vote for the same group |
| $1/n_{\text{groups}}$ | Votes spread uniformly (no agreement) |

The consensus score is computed at two resolutions:
- **Fine (cell-type level):** each cell type is its own group — harder to achieve high scores.
- **Lineage level:** cell types grouped by data-driven hierarchical clustering — coarser, easier to achieve agreement.

---

## 3. Permutation Testing

Both methods use permutation testing to assess whether observed coherence exceeds what is expected by chance. The key challenge is that peaks near promoters and in dense regulatory regions may show correlated specificity patterns for structural reasons unrelated to biology. The permutation null controls for this by matching on **genomic context**.

### Joint bin construction

Each genome-wide peak is classified into a **joint bin** defined by two properties:

**TSS distance bin** — distance to the nearest gene transcription start site:

| Bin label | Distance range |
|-----------|---------------|
| `core_prom` | 0 – 1 kb |
| `prox_prom` | 1 – 5 kb |
| `near_reg` | 5 – 50 kb |
| `distal` | 50 – 200 kb |
| `long_range` | > 200 kb |

**Peak density bin** — local peak density (number of other peaks within a 50 kb window), divided into quantile-based bins.

The joint bin is the interaction (cross-product) of these two:

$$\text{joint\_bin}(p) = \text{tss\_bin}(p) \times \text{density\_bin}(p)$$

### Rank concordance permutation (bin-matched peak sampling)

For a locus with $k$ peaks distributed across bins $\{b_1, b_2, \ldots\}$ with counts $\{n_1, n_2, \ldots\}$:

1. For each bin $b_j$, sample $n_j$ peaks uniformly from all genome-wide peaks in bin $b_j$.
2. Combine to form a permuted set of $k$ peaks.
3. Compute concordance on the permuted set using the genome-wide ranked matrix.
4. Repeat $N_{\text{perm}}$ times.

This preserves the **exact bin composition** (e.g., "2 promoter peaks + 3 distal peaks in a medium-density region") of the observed locus.

### Consensus permutation (locus-matched whole-locus sampling)

The consensus method samples **whole loci** rather than individual peaks, to preserve within-locus spatial structure:

1. Define genome-wide loci by clustering nearby peaks (within 100 kb on the same chromosome).
2. Exclude all GWAS-overlapping loci from the pool.
3. For each permutation:
   - Compute the observed locus's bin proportion vector: $\mathbf{p}_{\text{obs}} = (n_1/k, \; n_2/k, \; \ldots)$
   - Compute the L1 (Manhattan) distance to every non-GWAS locus's bin proportions:
   
   $$d(\text{obs}, \ell) = \sum_j \left| p_{\text{obs},j} - p_{\ell,j} \right|$$
   
   - Select the top-20 closest-matching loci (lowest $d$).
   - Randomly pick one; sample $k$ peaks from it.
   - Compute consensus on the sampled peaks using the genome-wide L2 matrix.
4. Repeat $N_{\text{perm}}$ times.

### Statistical summary

For both methods, the null distribution $\{\theta_1^*, \theta_2^*, \ldots, \theta_{N_\text{perm}}^*\}$ is summarised as:

**Z-score:**

$$Z = \frac{\theta_{\text{obs}} - \bar{\theta}^*}{\text{SD}(\theta^*)}$$

**Empirical p-value (rank concordance, with pseudocount):**

$$p = \frac{\sum_{i=1}^{N_\text{perm}} \mathbf{1}[\theta_i^* \geq \theta_{\text{obs}}] + 1}{N_\text{perm} + 1}$$

**Empirical p-value (consensus, without pseudocount):**

$$p = \frac{\sum_{i=1}^{N_\text{perm}} \mathbf{1}[\theta_i^* \geq \theta_{\text{obs}}]}{N_\text{perm, eff}}$$

where $N_\text{perm, eff}$ is the number of non-NA permutations (some permutations may fail if no matching loci are available).

---

## 4. Locus Categorisation

After permutation testing, loci are classified into interpretive categories:

| Category | Criteria |
|----------|----------|
| **Single_peak** | Only 1 peak in the locus (no pairwise comparison possible) |
| **Highly_coherent_specific** | $Z > 2$, $p < 0.05$, and strong specificity score (dominant_score > 0.5 for concordance, or consensus_lineage > 0.5 for consensus) |
| **Coherent_moderate** | $Z > 2$ and $p < 0.05$, but specificity not exceptionally strong |
| **Incoherent** | $Z < 1$ — peaks show no more agreement than expected by chance |
| **Intermediate** | Between coherent and incoherent |

---

## 5. Key Differences Between the Two Methods

| Aspect | Rank Concordance | Consensus |
|--------|-----------------|-----------|
| **What it measures** | Agreement on full cell-type ranking | Agreement on which single group dominates |
| **Metric** | Mean pairwise Spearman $\rho$ | Weighted plurality fraction |
| **Range** | [−1, +1] | [1/$n_\text{groups}$, 1] |
| **Sensitivity** | Sensitive to ordering across ALL cell types | Only sensitive to the top group per peak |
| **Weights** | Unweighted (all peaks equal) | Peaks weighted by specificity × accessibility |
| **Resolution** | Cell-type level only | Cell-type level AND lineage level |
| **Permutation unit** | Individual peaks (by bin) | Whole loci (by bin composition) |
| **Permutation matching** | Exact bin counts | L1 distance on bin proportions (top-20 pool) |

### When they diverge

- **High concordance, moderate consensus:** Peaks agree on the full ordering but no single cell type captures most of the weight (e.g., two cell types are consistently co-ranked at the top).
- **High consensus, moderate concordance:** Peaks agree that one cell type is dominant, but disagree on the ordering of the remaining cell types.
- **Both high:** Strong, consistent cell-type specificity across all peaks in the locus.
- **Both low:** No regulatory coherence — peaks point in different directions.
