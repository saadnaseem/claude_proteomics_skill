# Stage 1 — QC summary

## Goal
Sanity-check the data before any biological interpretation. Surface anomalies up front.

## What to compute

In a notebook cell, compute and display:

1. **Coverage**: `n_proteins`, `n_with_protein_name` (non-null `Protein.Names`), `n_unique_protein_groups`.
2. **Missingness**: count rows where any of `log2_mean_A`, `log2_mean_B`, `log2_std_A`, `log2_std_B` is NaN. If >5% missing, flag it.
3. **Std dev sanity**: rows where `log2_std_A == 0 OR log2_std_B == 0` (single-replicate inferred → unreliable t-test). Count and flag.
4. **p-value distribution**: histogram of `p-value` (50 bins). A flat distribution with a spike near 0 = healthy. A U-shape or flat-only = something's wrong with the test or normalization.
5. **log2FC distribution**: histogram and percentiles (1%, 5%, 25%, 50%, 75%, 95%, 99%) of `log2_Fold_change_A/B`. Note any extreme outliers (|log2FC| > 8).
6. **Volcano preview**: scatter of `log2_Fold_change_A/B` vs `-log10(p_adjusted(BH))`, no thresholds applied yet. Save to `figures/01_volcano_preview.png`.
7. **Significance counts** at standard cutoffs:
   - p_adj < 0.05: total, up (log2FC>0), down (log2FC<0)
   - p_adj < 0.05 AND |log2FC| > 1: total, up, down
   - p_adj < 0.01 AND |log2FC| > 1: total, up, down

## Output to user

One line per metric. Numbers-first. Example:
```
Stage 1 (QC):
  Proteins: 3,142   |   With names: 2,891 (92%)   |   Missingness: 47 rows (1.5%)
  Single-rep proteins (std=0): 12 (flagged but kept)
  p-value distribution: healthy (spike near 0, flat tail)
  log2FC range: [-7.2, +6.8]   |   |log2FC|>8 outliers: 0
  Significance @ p_adj<0.05, |log2FC|>1: 247 total (147 up, 100 down)
  Volcano preview saved: figures/01_volcano_preview.png
```

If anything is anomalous (high missingness, U-shaped p-dist, etc.), surface it as a separate "⚠️ Anomalies:" block and recommend whether to proceed or stop. **Do not pause unless there's a real problem** — this is a no-pause stage by default.
