# Stage 9 — AlphaFold structural inspection

## Goal
For the top biological hits, inspect predicted structures and confidence scores. Look for: high-confidence functional domains, low-confidence regions (often regulatory / disordered), known binding sites.

## ✋ Pre-stage gate

```
Stage 9 — AlphaFold structures:
  Will fetch top N hits ranked by impact (|log2FC| × -log10(p_adj))
  Default N = 10 structures (~2min download, ~30MB)
  Will compute pLDDT distributions and render with py3Dmol
Proceed with N=10? (yes / set N=X / skip)
```

## API
- AlphaFold DB: `https://alphafold.ebi.ac.uk/files/AF-<UniProt_acc>-F1-model_v4.pdb` (PDB) or `.cif`
- Metadata: `https://alphafold.ebi.ac.uk/api/prediction/<UniProt_acc>` (JSON)

## Selection criteria

Rank DE proteins by **impact score** = `|log2FC| × -log10(p_adj)`. Take top N (default 10), excluding:
- Hits with no AlphaFold entry (note in report)
- Trivial cases (e.g., ribosomal proteins where structures are well-known) — only if hypothesis-relevant

Also include any protein the user named in Stage 2 motivation, even if not in the top N.

## What to do per protein

1. Fetch the PDB and metadata. Cache to `annotations/alphafold/<acc>.pdb`.
2. Compute pLDDT statistics: mean, median, % high-confidence (pLDDT > 70), % very-low (pLDDT < 50).
3. Render with `py3Dmol` colored by pLDDT (default AF coloring), inline in notebook. Save HTML snapshot.
4. Annotate with UniProt features (Stage 5 data):
   - Active sites
   - Binding sites (substrate, cofactor, metal)
   - Modified residues (PTMs)
   - Domains (Pfam-style)
5. Render a static PNG via PyMOL/ChimeraX is too heavy; instead use `Bio.PDB` to compute simple metrics (radius of gyration, n_residues, n_chains) and a thumbnail via py3Dmol+selenium IF available, else just save the interactive HTML.

## Output

`tables/09_structural.csv`:
```
acc, gene_name, log2FC, p_adj, plddt_mean, plddt_median, pct_high_conf, n_active_sites, n_binding_sites, n_ptm_sites, alphafold_url
```

Per-protein markdown in notebook: structure thumbnail/HTML, key metrics, biological interpretation note.

## Output to user

```
Stage 9 — AlphaFold:
  Fetched: 10 structures
  Mean pLDDT: 87.3 (range 62.1 - 94.8)
  Notable:
    - <gene>: pLDDT 91.2, 2 active sites, predicted phosphorylation at Ser142 — see figures/09_<acc>.html
    - <gene>: pLDDT 68.4 (low), large disordered N-term — likely regulatory
  Saved: figures/09_*.html
```

## No ✋ gate after — proceed to Stage 10.
