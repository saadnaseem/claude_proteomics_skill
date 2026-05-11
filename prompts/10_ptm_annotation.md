# Stage 10 — UniProt PTM site annotation

## Goal
For DE proteins, list known post-translational modification sites from UniProt and check whether any are at functionally important residues (active sites, binding sites). Surface DE proteins with known regulatory PTMs as candidates for follow-up.

## Scope clarification
This input format does NOT contain site-level PTM quantitation (no `Modified.Sequence` column). So this stage is **annotation only**: we attach UniProt's curated PTM knowledge to the DE set. If the user later supplies FragPipe or MaxQuant `Phospho (STY)Sites.txt` outputs, that's a separate v0.2 extension.

## What to do

1. From Stage 5 annotations, extract `ptm_features` for every DE protein. Categories of interest:
   - Phosphorylation (`MOD_RES` with phospho descriptions)
   - Acetylation (often N-terminal or lysine)
   - Methylation
   - Glycosylation (`CARBOHYD`, less common in microbes but check)
   - Lipidation (`LIPID`)
   - Disulfide bonds
   - Cleavage sites (`PROPEP`, `SIGNAL`)
2. **For each PTM**, record: protein, gene_name, position, type, description, and whether the position falls within ±5 residues of an active or binding site.
3. **Highlight regulatory PTMs**: any phospho/acetyl on a kinase/phosphatase, sigma factor, transcription factor, or response regulator in the DE set is high-priority. Flag in the report.
4. Save `tables/10_ptm_annotations.csv`.
5. **Bar chart**: count of PTMs per type in up vs down sets. Save `figures/10_ptm_summary.png`.

## Microbial-specific notes

- Bacterial PTMs are less catalogued than eukaryotic. Don't be surprised by sparse coverage.
- Phospho-Ser/Thr/Tyr in bacteria often signal two-component-system response regulators.
- Acetylation of central metabolism enzymes (e.g., AceB, GltA) is a known regulatory mechanism — flag if present.
- N-terminal Met excision and acetylation are very common; not biologically interesting unless specific to DE pattern.

## No ✋ gate — proceed to Stage 11.

Output to user:
```
Stage 10 — PTMs:
  DE proteins with annotated PTMs: 89/247 (36%)
  Most common: phosphorylation (52 sites, 28 proteins)
  Regulatory candidates (PTM near active site):
    - <gene>: Ser142 phospho, 4Å from active site Asp140
    - <gene>: Lys87 acetyl on TCA enzyme — known regulatory site
  Saved: tables/10_ptm_annotations.csv, figures/10_ptm_summary.png
```
