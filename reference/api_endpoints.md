# External API reference

All endpoints are unauthenticated unless noted. Always cache responses, always rate-limit, always retry with exponential backoff.

## UniProt
**Base**: `https://rest.uniprot.org/`

| Endpoint | Purpose | Example |
|---|---|---|
| `/uniprotkb/search` | Batch query | `?query=accession:(P12345 OR Q67890)&format=json&size=500&fields=accession,gene_names,protein_name,go,xref_kegg,ec,xref_string,xref_alphafolddb,cc_subcellular_location,ft_mod_res,ft_lipid,ft_carbohyd,organism_id` |
| `/uniprotkb/<acc>` | Single record | `/uniprotkb/P12345?format=json` |
| `/proteomes/search` | Find proteome ID | `?query=organism_id:160488&format=json` |
| `/taxonomy/search` | Resolve organism name → taxid | `?query=Pseudomonas+putida+KT2440&format=json&size=5` |

Rate limit: ~5 req/sec. Batch endpoints accept up to 500 accessions per query.

## STRING
**Base**: `https://string-db.org/api/`

| Endpoint | Purpose |
|---|---|
| `/json/get_string_ids` | Map identifiers (UniProt → STRING) |
| `/json/network` | Get interaction edges |
| `/json/enrichment` | (we usually skip — Stage 6 covers it) |
| `/json/version` | Get DB version (record in report) |

Required parameter: `species=<NCBI taxid>`. Confidence threshold: `required_score=700` (high) by default.

Rate limit: ~1 req/sec, ~2,000 proteins per call.

## KEGG
**Base**: `https://rest.kegg.jp/`

| Endpoint | Purpose | Example |
|---|---|---|
| `/find/genome/<query>` | Find KEGG organism code | `/find/genome/Pseudomonas+putida` |
| `/list/pathway/<org>` | All pathways for organism | `/list/pathway/ppu` |
| `/link/<org>/pathway` | Pathway → gene mapping | `/link/ppu/pathway` |
| `/get/<pathway_id>/image` | Pathway diagram PNG | `/get/ppu00020/image` |
| `/conv/uniprot/<org>` | UniProt ↔ KEGG gene ID | `/conv/uniprot/ppu` |

Rate limit: ~3 req/sec. They block abusive IPs.

## AlphaFold DB
**Base**: `https://alphafold.ebi.ac.uk/`

| Endpoint | Purpose |
|---|---|
| `/api/prediction/<UniProt_acc>` | JSON metadata |
| `/files/AF-<UniProt_acc>-F1-model_v4.pdb` | PDB structure file |
| `/files/AF-<UniProt_acc>-F1-model_v4.cif` | mmCIF (preferred for >2000 residues) |
| `/files/AF-<UniProt_acc>-F1-predicted_aligned_error_v4.json` | PAE matrix |

No published rate limit. Be polite: ~2 req/sec.

## kgMicrobe / Monarch
**Base**: `https://api.monarchinitiative.org/api/`

| Endpoint | Purpose |
|---|---|
| `/bioentity/taxon/<id>/genes` | Genes annotated to taxon |
| `/bioentity/gene/<id>/phenotypes` | Phenotypes for gene |
| `/search/entity/<query>` | Free-text search |

Direct kg-microbe data dumps: `https://kg-hub.berkeleybop.io/kg-microbe/` (TSV/JSON). Often more reliable than the live API.

If unavailable, fall back to: BacDive (`https://bacdive.dsmz.de/`), KEGG organism page, IMG/M.

## QuickGO (GO annotations)
**Base**: `https://www.ebi.ac.uk/QuickGO/services/`

Useful when UniProt's GO data is sparse:
- `/annotation/search?taxonId=<taxid>&geneProductId=<acc>&format=json`

## PubMed (literature search via Entrez)
**Base**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`

For Deep Research stage when WebSearch isn't enough:
- `/esearch.fcgi?db=pubmed&term=<query>&retmode=json&retmax=20&sort=relevance`
- `/esummary.fcgi?db=pubmed&id=<pmid>&retmode=json`

## STRING and UniProt versioning
Always record the version used in `manifest.json`:
- UniProt: hit `https://rest.uniprot.org/configure/uniprotkb/release_summary` for current release
- STRING: hit `https://string-db.org/api/json/version`
- KEGG: KEGG release date is on https://www.kegg.jp/kegg/docs/relnote.html
- AlphaFold: pin to v4 (current as of skill creation)
