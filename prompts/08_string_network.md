# Stage 8 — STRING interaction network + hub analysis

## Goal
Build the protein-protein interaction subgraph for the DE set, identify hub proteins, and detect dense functional modules.

## ✋ Pre-stage gate

Before starting (this stage hits the STRING API and can take 30-60s for large sets):
```
Stage 8 — STRING network:
  Will query STRING for 247 DE proteins (organism taxid <id>)
  Estimated: ~30s, ~3 API calls
  Confidence threshold: 0.7 (high)
Proceed? (yes / lower threshold / skip stage)
```

## API
STRING REST: `https://string-db.org/api/[output_format]/[method]?[parameters]`

Key endpoints:
- Get IDs: `https://string-db.org/api/json/get_string_ids?identifiers=<UniProt accs separated by %0d>&species=<taxid>`
- Network: `https://string-db.org/api/json/network?identifiers=<STRING IDs>&species=<taxid>&required_score=700&network_type=physical`
- Enrichment: `https://string-db.org/api/json/enrichment?identifiers=<>&species=<taxid>` (we mostly skip — Stage 6 covered this)

Cite STRING version in the manifest and report.

## What to do

1. **Map UniProt → STRING IDs** for the DE set (single batch call).
2. **Fetch interactions** at `required_score=700`. Default to `network_type='physical'` (PPIs) but optionally also do `'functional'` for context.
3. **Build NetworkX graph**:
   - Nodes = proteins, attrs: gene_name, log2FC, p_adj, direction (up/down)
   - Edges = STRING interactions, attrs: combined_score, evidence channels
4. **Hub analysis**:
   - Compute degree centrality and betweenness for each node
   - Identify top 10 hubs (high degree among DE set — likely regulators / scaffolds)
   - Save `tables/08_hubs.csv`
5. **Module detection**:
   - Run `networkx.algorithms.community.louvain_communities()` for community detection
   - For each community ≥5 nodes, list dominant GO terms (intersect with Stage 5 annotations)
   - Save `tables/08_communities.csv`
6. **Visualization**:
   - Static: `nx.draw_spring_layout`, node color by direction (red up, blue down), node size by `|log2FC|`, edge width by combined_score. Save `figures/08_network.png`.
   - Interactive: build a `pyvis.Network()` HTML at `figures/08_network.html` — far more useful for exploration. Add tooltips with gene name + log2FC.

## No ✋ gate after — proceed to Stage 9.

Output to user:
```
Stage 8 — STRING network:
  Queried: 247 proteins → 218 mapped to STRING
  Edges: 1,432 at score≥700
  Hubs (top 5): <name (degree, log2FC)>, ...
  Communities: 12 (3 with ≥10 nodes)
  Largest community: 28 proteins, dominant GO: amino acid biosynthesis
  Static: figures/08_network.png   Interactive: figures/08_network.html
```
