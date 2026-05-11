# Stage 5 — UniProt annotation

## Goal
Enrich the DE protein set with functional metadata from UniProt: GO terms, KEGG pathways, EC numbers, COG categories, gene names, descriptions, subcellular location.

## API
UniProt REST: `https://rest.uniprot.org/uniprotkb/<accession>?format=json` (one accession at a time) OR
batch: `https://rest.uniprot.org/uniprotkb/search?query=accession:(P12345 OR Q67890)&format=json&size=500&fields=accession,gene_names,protein_name,go,xref_kegg,ec,xref_string,xref_alphafolddb,cc_subcellular_location,ft_mod_res,ft_lipid,ft_carbohyd,organism_id`

**Always use batch mode** — 100 accessions per query, paginated. Rate limit yourself to ~5 req/sec.

## What to do

1. Initialize `requests-cache`:
   ```python
   import requests_cache
   requests_cache.install_cache('annotations/uniprot_cache', backend='sqlite', expire_after=86400*30)
   ```
2. Build the accession list: union of up + down DE proteins. Take the FIRST accession from each `Protein.Group` (semicolon-split).
3. Query UniProt in batches of 100. Use `tenacity` for retries on 429/500.
4. Parse responses into a DataFrame with columns:
   - accession, gene_name (primary), protein_name, ec_number, kegg_ids (list), go_bp (list), go_mf (list), go_cc (list), string_id, alphafold_id, subcellular_location, ptm_features (list of dicts: position, type, description)
5. Save to `tables/05_annotations.csv` (lists JSON-serialized).
6. Compute coverage: % of DE proteins with at least one GO term, KEGG ID, AlphaFold structure. Report in notebook.

## Special cases
- **Obsolete accessions**: UniProt may redirect. Follow `secondaryAccession` and update the table.
- **Trembl vs SwissProt**: prefer SwissProt entries. Note in coverage report how many are Trembl-only (less curated).
- **Protein groups**: for each group, the first accession is the representative; record all members in a `group_members` column for reference.

## No ✋ gate — proceed to Stage 6.

Output to user:
```
Stage 5 — UniProt annotation:
  Queried: 247 accessions
  Cache hits: 0 (first run) / N (re-run)
  GO coverage: 232/247 (94%)
  KEGG coverage: 198/247 (80%)
  AlphaFold coverage: 241/247 (98%)
  SwissProt: 215   Trembl-only: 32
```
