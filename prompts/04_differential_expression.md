# Stage 4 — Differential expression call

## Goal
Apply DE cutoffs to produce the canonical up/down sets that all downstream stages consume.

## Default cutoffs
- `p_adjusted(BH) < 0.05`
- `|log2_Fold_change_A/B| > 1` (i.e., 2× change)

These are conservative defaults appropriate for microbial proteomics with replicates. Override only if user requests.

## What to do

1. Compute the boolean masks:
   ```python
   sig_up   = (df['p_adjusted(BH)'] < 0.05) & (df['log2_Fold_change_A/B'] >  1)
   sig_down = (df['p_adjusted(BH)'] < 0.05) & (df['log2_Fold_change_A/B'] < -1)
   sig_any  = sig_up | sig_down
   ```
2. Save the DE table:
   - `tables/04_de_proteins.csv` with columns: Protein, Protein.Group, Protein.Names, Protein.Description, log2_Fold_change_A/B, p-value, p_adjusted(BH), direction (up/down/ns)
3. Annotated volcano plot:
   - Color: significant up (red), down (blue), ns (grey)
   - Label top 15 by `|log2FC| × -log10(p_adj)` ranking (impact score)
   - Save `figures/04_volcano.png` (300 DPI) and embed in notebook
4. MA plot (mean log2 intensity vs log2FC):
   - x = `(log2_mean_A + log2_mean_B) / 2`, y = `log2_Fold_change_A/B`
   - Same color scheme
   - Save `figures/04_ma_plot.png`
5. Summary table cell:
   ```
   Direction  | Count | % of detected
   Up (sig)   |  147  | 4.7%
   Down (sig) |  100  | 3.2%
   NS         | 2895  | 92.1%
   ```

## ✋ Gate

```
Stage 4 — DE call complete:
  Cutoffs: p_adj < 0.05, |log2FC| > 1
  Up:   147 proteins
  Down: 100 proteins
  Top 5 up:   <names with log2FC>
  Top 5 down: <names with log2FC>
  Volcano: figures/04_volcano.png

Proceed with these cutoffs, or adjust? (yes / set p_adj=X, log2FC=Y)
```

If user adjusts, recompute and re-gate.

## Notes
- Record applied cutoffs in `manifest.json` under `de_cutoffs`.
- Pass the up/down protein lists (UniProt accessions, deduplicated, taking the first accession of each Protein.Group) to all downstream stages.
