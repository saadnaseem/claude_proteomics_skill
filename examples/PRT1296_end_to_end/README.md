# PRT1296 — End-to-end 2×3 factorial proteomics pipeline

A complete, reproducible reference implementation of the proteomics-agent pipeline applied to a 2 (strain) × 3 (medium) factorial design. Drop your own data in `data/`, edit `config/analysis.yaml`, and run `bash scripts/run_all.sh`.

## What this example demonstrates

This is the most complete worked example in the skill repo. It runs the full SKILL.md pipeline on a real DIA-NN MaxLFQ output (1021 quantified proteins, 6 sample groups × 3 replicates, 18 LC-MS runs) and produces:

- QC: per-run correlation heatmap, intensity boxplot, dendrogram, CV table, missingness summary
- Exploratory: annotated PCA, hierarchical clustering, per-protein variance partitioning
- Differential expression: 11 pairwise t-test integrations + per-protein two-factor ANOVA (strain × medium interaction)
- DEP sets: UpSet plots (custom plotter — no `upsetplot` dependency), per-set membership tables
- Heatmap: z-scored, hierarchically clustered, with strain/medium annotation bars
- Enrichment: GO BP/MF/CC + KEGG ORA via UniProt-GOA + KEGG REST APIs (works for any organism KEGG knows, including those not supported by gProfiler — e.g. *Pseudomonas putida* KT2440)
- Tolerance-module scoring: 15 curated metabolic / stress / membrane modules, two-way ANOVA per module
- Interaction profiles: categorize interaction-significant proteins into constitutive / buffered / hyper-induced / reciprocal patterns
- Per-cluster + per-DEP-set enrichment (extended analysis)
- Stress-amplitude scatter (HGL vs KT log2FC per protein)
- Annotated volcano panels with gene labels

## Directory layout

```
PRT1296_end_to_end/
├── README.md                  this file
├── data_analysis_plan.md      the full plan that produced this pipeline (13 sections)
├── config/
│   ├── analysis.yaml          sample → strain × medium mapping, contrasts, thresholds, palette
│   └── tolerance_modules.yaml 15 curated gene-symbol modules
└── scripts/
    ├── _common.py             shared helpers (config loader, metadata, log2-safe)
    ├── 00_load_and_metadata.py
    ├── 01_qc.py
    ├── 02_exploratory.py
    ├── 03_two_factor_anova.py
    ├── 04_dep_sets.py         self-contained UpSet plotter (no upsetplot dep)
    ├── 05_heatmap.py
    ├── 06_enrichment.py       UniProt-GOA + KEGG REST + custom hypergeometric ORA
    ├── 07_modules.py
    ├── 08_interaction_profiles.py
    ├── 09_extended_analysis.py  per-cluster + per-set enrichment, stress-amplitude, volcanoes
    └── run_all.sh             reproducible end-to-end runner
```

## How to adapt to your data

1. **Place your DIA-NN protein-level CSV(s)** in a sibling `data/` directory and symlink the canonical inputs:
   - `data/proteins_long.csv` — long format with columns `Protein.Group, Protein, Protein.Description, Sample, Replicate, Counts_sum`
   - `data/proteins_summary.csv` — per-group means/std/CV (optional, used for sanity checks)
   - `data/top3_proteins.csv` — Top3 quant (optional, for absolute-abundance)
   - `data/t_test_xlsx/` — directory of pairwise t-test xlsx files with sheet `Full t-test output` containing columns `Protein, Protein.Group, Protein.Names, Protein.Description, log2_mean_<A>, log2_mean_<B>, log2_std_<A>, log2_std_<B>, t-test_stat, p-value, p_adjusted(BH), log2_Fold_change_A/B`

2. **Edit `config/analysis.yaml`**:
   - Set `project.organism`, `taxon_id`, `kegg_org` (KEGG organism code, e.g. `ppu` for P. putida, `eco` for E. coli K-12, `bsu` for B. subtilis)
   - Replace the `samples` block with your sample-id → strain × medium mapping
   - Update the `contrasts` block with your A-vs-B pairs (one entry per pairwise t-test xlsx file)
   - Adjust `thresholds.q_strict` / `log2fc_strict` if you want a different cutoff
   - Replace the `plot.palette_strain` and `palette_medium` colour maps to match your factor levels

3. **Edit `config/tolerance_modules.yaml`** to add or remove gene-symbol modules. Each module has either a `genes` list of exact symbols (lowercase-matched) or a `gene_prefixes` list (matches anything starting with the prefix — useful for ribosomal proteins, `rps`/`rpl`/`rpm`).

4. **Create the Python environment** once:
   ```
   python3 -m venv .venv
   source .venv/bin/activate
   pip install pandas numpy scipy statsmodels scikit-learn matplotlib seaborn \
               openpyxl requests pyyaml pingouin bioservices
   ```

5. **Run**:
   ```
   bash scripts/run_all.sh
   ```
   Outputs go to `outputs/figures/`, `outputs/tables/`, `outputs/enrichment/`.

## Design choices worth knowing

- **Background for ORA** = the actually-quantified proteins (whatever was in `proteins_long.csv`), *not* the full ORFeome. This corrects for detection-bias in single-shot DIA.
- **Enrichment goes through UniProt + KEGG REST** rather than gProfiler. Slower on first run but caches everything to `data/uniprot_annotation.tsv`, `data/kegg_pathway_map.tsv`, `data/go_ontology_map.tsv` — subsequent runs are fast. Works for any organism KEGG knows.
- **Variance partitioning** uses Type-II SS via residual subtraction (`SS_strain = r_strain² − r_inter²`) so it's correct for balanced designs and the per-protein percentages sum to 100 %.
- **UpSet plotter** is self-contained (gridspec + scatter) — no dependency on `upsetplot`, which has matplotlib-version compatibility issues.
- **Interaction categorizer** rule order: (1) sign-flip → reciprocal, (2) amplitude difference → buffered / hyper-induced, (3) baseline shift → constitutive. This order avoids over-calling constitutive-shift when one strain has a much bigger stress response.
- **DIA-NN protein-group accessions** that join multiple UniProt IDs with `;` (e.g. `Q877Q0;Q88FK3`) are cleaned to the first accession before UniProt batch queries.

## Provenance

Built 2026-05-04 / audited 2026-05-05 on the PRT1296 study (*P. putida* KT2440 wild-type vs HGL1175 ALE-evolved strain across glucose / M9+cholinium-lysinate / M9+butylamine, 35 % v/v hydrolysate). The cross-validation against the existing JBEI t-test pipeline confirmed all 11 contrast DEP counts matched exactly.

See `data_analysis_plan.md` for the full 13-section plan including biological hypotheses, statistical thresholds, interpretation framework, and risks/caveats. The plan is itself reusable — generalize section 1 (background) to your organism and stressor; the rest of the structure is data-agnostic.
