# Stage 6 — GO / KEGG over-representation analysis

## Goal
Identify functional categories statistically over-represented in up-regulated and down-regulated sets vs the background of detected proteins.

## Background set
**Critical**: the background is **all detected proteins** (the full xlsx), not the entire genome. Using the genome inflates significance.

## Tools
- `gseapy.enrich` for GO/KEGG with custom background
- Or `goatools.go_enrichment.GOEnrichmentStudyNS` for GO with full DAG awareness

Prefer `gseapy.enrich` for simplicity. Use `goatools` only if user wants term-DAG visualization.

## What to do

1. **Build gene sets**:
   - Up: list of UniProt accessions (or gene names — `gseapy` accepts either, gene names preferred)
   - Down: list
   - Background: all detected proteins
2. **Map UniProt → gene symbol** using the annotation table from Stage 5. For organisms `gseapy` doesn't natively support, use custom GMT files built from the Stage 5 GO/KEGG columns.
3. **GO enrichment** for up and down separately, three ontologies (BP, MF, CC):
   - `gseapy.enrich(gene_list=up_genes, gene_sets='GO_Biological_Process_2023', background=bg_genes, organism='<organism>', outdir='tables/06_go_up_bp/')`
   - Filter results: `Adjusted P-value < 0.05`, term size 5–500 (avoid trivial / overly broad terms)
   - Save top 30 to CSV
4. **KEGG enrichment** for up and down:
   - Use the organism's KEGG code (e.g., `ppu` for *P. putida* KT2440)
   - Build pathway gene sets from KEGG REST: `https://rest.kegg.jp/link/<kegg_org_code>/pathway`
   - Same filtering
5. **Visualization**:
   - Bar plot of top 15 enriched terms, up and down (separate panels), -log10(p_adj) on x-axis
   - Save `figures/06_enrichment_up.png`, `figures/06_enrichment_down.png`
   - Optional: dotplot via `gseapy.dotplot`
6. **Self-critique**: read `prompts/critic.md` and ask: are the enriched terms consistent with the Stage 2 hypothesis? If not, note the discrepancy in the notebook and flag it for the synthesis report.

## Output table format

`tables/06_enrichment.csv`:
```
direction, source, term_id, term_name, n_overlap, n_term, n_background, p_value, p_adj, fold_enrichment, genes
up, GO_BP, GO:0006099, tricarboxylic acid cycle, 8, 12, 3142, 1.2e-7, 4.3e-5, 18.4, "gltA;sucA;..."
```

## No ✋ gate — proceed to Stage 7.

Output to user:
```
Stage 6 — Enrichment:
  GO BP enriched (up):   23 terms (top: TCA cycle, p_adj=4.3e-5)
  GO BP enriched (down): 18 terms (top: ribosome biogenesis, p_adj=2.1e-4)
  KEGG enriched (up):    7 pathways
  KEGG enriched (down):  4 pathways
  Hypothesis consistency: ✓ (TCA cycle up, ribosome down — matches H)
```
