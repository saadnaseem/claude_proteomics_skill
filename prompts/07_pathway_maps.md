# Stage 7 — KEGG pathway maps with up/down overlay

## Goal
For each significantly enriched KEGG pathway from Stage 6, overlay the DE proteins onto the canonical pathway diagram, colored by direction and magnitude.

## Tools
- `bioservices.KEGG` for pathway data
- `gseapy.scripts` has `pathways_overlay` helpers, OR
- Direct KEGG REST: `https://rest.kegg.jp/get/<pathway_id>/image` for the static map, plus per-gene coloring via the KEGG Pathview-style URL: `https://www.kegg.jp/pathway/<pathway_id>+<color_codes>`

## What to do

1. **Select pathways**: take the union of (top 5 enriched up + top 5 enriched down + any pathway specifically named in the Stage 2 hypothesis). Cap at 10 to keep runtime reasonable.
2. **For each pathway**:
   - Fetch pathway gene mapping: `KEGG().link('<pathway_id>', '<kegg_org_code>')`
   - Build color codes: red (#FF6666) for up-regulated genes in the DE set, blue (#6666FF) for down, white for present-but-not-DE
   - Either:
     - (a) Build KEGG color URL and fetch the rendered image
     - (b) Use `gseapy.kegg_pathway` if available
     - (c) As fallback, render with `networkx` showing only the DE-overlapping subgraph
   - Save `figures/07_kegg_<pathway_id>.png`
3. **Build a summary table** `tables/07_pathway_overlay.csv`:
   - pathway_id, pathway_name, n_genes_in_pathway, n_de_up, n_de_down, n_de_total, % coverage, summary
4. **Per-pathway markdown** in the notebook: pathway name, image, list of up/down genes with log2FC

## Special cases
- **Organism not in KEGG**: skip and log to `manifest.json:degraded_stages`. Continue to Stage 8.
- **Pathway image fetch fails**: render the gene-level subgraph with networkx as a fallback.

## No ✋ gate — proceed to Stage 8.

Output to user:
```
Stage 7 — Pathway maps:
  Rendered: 8 pathways
  Notable:
    - ppu00020 (TCA cycle): 8/12 enzymes upregulated (gltA, sucA, ...)
    - ppu03010 (Ribosome): 12/54 ribosomal proteins downregulated
  Saved: figures/07_kegg_*.png
```
