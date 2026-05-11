# Proteomics analysis: {{ORGANISM_NAME}} — {{CONDITION_A}} vs {{CONDITION_B}}

**Run ID**: {{RUN_ID}}
**Date**: {{DATE}}
**Input**: `{{INPUT_PATH}}` (sha256 `{{INPUT_SHA}}`, {{N_PROTEINS}} proteins)
**Notebook**: `notebook.ipynb`

---

## TL;DR

- Conditions: A = {{CONDITION_A}} vs B = {{CONDITION_B}} ({{ORGANISM_NAME}}, taxid {{TAXID}}).
- DE: **{{N_UP}} up**, **{{N_DOWN}} down** at p_adj<{{P_ADJ}}, |log2FC|>{{LOG2FC}}.
- Headline: {{HEADLINE_FINDING}}
- Strongest follow-up: {{TOP_RECOMMENDATION}}

## Methods

Pre-computed t-test xlsx loaded with pandas. DE called at p_adj<{{P_ADJ}} (BH), |log2FC|>{{LOG2FC}}. Annotation via UniProt REST ({{UNIPROT_VERSION}}). Enrichment via gseapy ({{GSEAPY_VERSION}}) with detected-protein background. Pathway maps via KEGG REST ({{KEGG_DATE}}). Interaction network via STRING DB v{{STRING_VERSION}} at confidence ≥0.7. Structural inspection via AlphaFold v4. Full package versions in `manifest.json`. Random seed: 42.

## Quality control

- Coverage: {{N_PROTEINS}} proteins, {{PCT_NAMED}}% with annotated names.
- Missingness: {{PCT_MISSING}}%.
- Single-replicate proteins (std=0): {{N_SINGLE_REP}} (kept but flagged).
- p-value distribution: {{P_DIST_NOTE}}.

## Biological story

{{LEAD_PARAGRAPH_HYPOTHESIS_VS_FINDINGS}}

### Pathway-level findings

{{PATHWAY_PARAGRAPH}}

### Network analysis

{{NETWORK_PARAGRAPH}}

### Structural & PTM observations

{{STRUCTURE_PTM_PARAGRAPH}}

### Surprises and contradictions

{{SURPRISES_PARAGRAPH}}

## Hub proteins / candidate regulators

{{HUB_TABLE}}

## Falsification check

Stage 2 falsification criteria: {{FALSIFICATION_CRITERIA}}
Outcome: {{FALSIFICATION_OUTCOME}}

## Recommended follow-up

### Wet-lab
{{WETLAB_RECS}}

### Computational
{{COMP_RECS}}

### Highest-priority candidates for biochemical characterization
{{PRIORITY_CANDIDATES}}

## Caveats

{{CAVEATS}}

## Appendix

- **Notebook**: `notebook.ipynb`
- **Tables**: `tables/`
  - `04_de_proteins.csv` — full DE results
  - `05_annotations.csv` — UniProt annotations
  - `06_enrichment.csv` — enriched GO/KEGG terms
  - `07_pathway_overlay.csv` — pathway coverage summary
  - `08_hubs.csv`, `08_communities.csv` — network analysis
  - `09_structural.csv` — AlphaFold metrics
  - `10_ptm_annotations.csv` — known PTM sites on DE proteins
- **Figures**: `figures/`
- **Manifest**: `manifest.json` (inputs, params, package versions, degraded stages)

### Citations

{{CITATIONS}}
