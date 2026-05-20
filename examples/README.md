# Reference examples

Real-world reusable scripts adapted from the LBNL P. putida ALE proteomics study (HGL1175 vs KT2440 across 3 media). Use these as starting templates when you have a similar dataset.

## `run_pipeline.py` — single-comparison end-to-end pipeline

Runs Stages 0, 1, 4, 5, 6, 7, 8 of the SKILL.md pipeline (intake, QC, DE call, UniProt annotation, GO enrichment, KEGG pathway enrichment, STRING network) on a single t-test xlsx file.

Skips Stages 2 (hypothesis), 3 (deep research), 9 (AlphaFold), 10 (PTM), 11 (synthesis) — those are interactive and best done by Claude in conversation.

**Usage:**
```bash
conda activate proteomics-agent
python run_pipeline.py --config <path-to-config.json>
```

**Config example:**
```json
{
  "input_file": "/path/to/t-test.xlsx",
  "a_id": "SN_001",
  "b_id": "SN_002",
  "a_label": "Strain X / glucose",
  "b_label": "Strain Y / glucose",
  "run_name": "X_vs_Y_glucose",
  "taxid": 160488,
  "organism": "Pseudomonas putida KT2440",
  "kegg_org": "ppu",
  "de_p": 0.05,
  "de_fc": 1.0,
  "expand_p": 0.05
}
```

Outputs go to `~/proteomics_runs/proteomics_run_<timestamp>_<run_name>/` with `manifest.json`, `tables/`, `figures/`, `annotations/api_cache.sqlite`.

**Customizations to consider for your project:**
- Change `taxid` and `kegg_org` for non-P-putida organisms (e.g., `83333` + `eco` for E. coli K-12)
- Adjust `de_p`/`de_fc` for sparser data
- Bump `expand_p` to 0.10 if your dataset has very few hits at p_adj < 0.05

## `phase3_meta_analysis.py` — multi-comparison cross-analysis

Cross-comparison meta-analysis across N pairwise comparisons. Builds:
- Master DE matrix (all proteins × all comparisons)
- Constitutive vs conditional decomposition (multi-significance intersection)
- Buffered-response 2D scatter (WT response vs ALE response per protein)
- Regulator co-expression heatmap
- Composite-score tolerance candidate ranking

Edit the `COMPARISONS` dict at the top to point at your set of files. Edit the `regulators` dict to query genes specific to your biology.

## `build_meta_notebook.py` — emit a fully-populated notebook

Generates a Jupyter notebook (.ipynb) that reproduces the entire Phase 1 + Phase 3 analysis from raw data → final figures. Useful for handing the full analysis to a collaborator.

## `PRT1296_end_to_end/` — full modular pipeline (2×3 factorial reference)

The most complete worked example in the repo. A 10-script modular pipeline (`scripts/00_…09_` + `run_all.sh`) that runs every stage of the SKILL.md pipeline on a real 2 (strain) × 3 (medium) factorial DIA-NN dataset, with config-driven sample mapping, contrast definitions, and gene-symbol tolerance modules.

Includes:
- Self-contained UpSet plotter (no `upsetplot` dependency)
- UniProt-GOA + KEGG REST custom hypergeometric ORA (works for any KEGG-indexed organism — *P. putida*, *B. subtilis*, etc.)
- Per-protein strain × medium two-factor ANOVA with interaction-term FDR
- Curated 15-module tolerance scoring + per-module ANOVA
- Interaction-protein categorizer (reciprocal / buffered / hyper-induced / constitutive)
- Per-cluster and per-DEP-set enrichment
- Stress-amplitude scatter + annotated volcano panels
- Full 13-section `data_analysis_plan.md` you can crib for your own project

Use this as the starting point for any new 2-factor (or higher) proteomics study. Drop your DIA-NN CSV in `data/`, edit `config/analysis.yaml`, run `bash scripts/run_all.sh`.

## When to use these vs the agent

- **Use the agent (`/proteomics-agent`)** when you need the full hypothesis + deep research + synthesis story — it does the interactive interpretation steps the scripts don't.
- **Use these scripts** when you have many similar files and want batch processing without re-running the interactive layers, or as starting code for an analysis you want to build on.

The agent uses a similar pipeline internally; these scripts are the externalized, reusable form.
