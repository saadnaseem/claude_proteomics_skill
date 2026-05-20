# PRT1296 — End-to-end 2×3 factorial proteomics pipeline

A complete, reproducible reference implementation of the proteomics-agent pipeline applied to a 2 (strain) × 3 (medium) factorial design. The cleanest starting point in the repo for a *new* multi-factor study: drop your own data in `data/`, edit `config/analysis.yaml`, run `bash scripts/run_all.sh`.

## What this example demonstrates

This is the most complete worked example in the skill repo. It runs the full SKILL.md pipeline on a real DIA-NN MaxLFQ output (1021 quantified proteins, 6 sample groups × 3 replicates, 18 LC-MS runs) and produces:

- **QC** — per-run correlation heatmap, intensity boxplot, dendrogram, CV table, missingness summary
- **Exploratory** — annotated PCA, hierarchical clustering, per-protein variance partitioning (strain / medium / interaction / residual)
- **Differential expression** — 11 pairwise t-test integrations + per-protein two-factor ANOVA with `strain × medium` interaction term, BH-adjusted
- **DEP sets** — UpSet plots (custom plotter — no `upsetplot` dependency), per-set membership tables
- **Heatmap** — z-scored, hierarchically clustered, with strain/medium annotation bars
- **Enrichment** — GO BP/MF/CC + KEGG ORA via UniProt-GOA + KEGG REST APIs (works for any organism KEGG indexes, including those not supported by gProfiler — e.g. *Pseudomonas putida* KT2440)
- **Tolerance-module scoring** — 15 curated metabolic / stress / membrane modules, two-way ANOVA per module
- **Interaction profiles** — categorize interaction-significant proteins into constitutive / buffered / hyper-induced / reciprocal patterns
- **Per-cluster + per-DEP-set enrichment** (extended analysis)
- **Stress-amplitude scatter** (HGL vs KT log2FC per protein), annotated volcano panels

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

When you run the pipeline, the following directories are created next to `scripts/` and `config/`:

```
data/                          # your inputs (or symlinks to them)
outputs/
├── figures/                   # F1_PCA.png … F13_stress_amplitude.png (~17–20 PNGs)
├── tables/                    # protein matrices, ANOVA, DEP master, modules, interactions
└── enrichment/                # one TSV per (contrast × ontology × direction); cluster + set summaries
```

---

## Setup (one-time, ~5 min)

The example reuses the skill's conda environment, which is created when you install the skill — see the top-level [README](../../README.md#installation). If you already installed the skill, you're done.

If you want to use this example standalone (without the skill), the minimal env is:

```bash
conda create -n proteomics-agent -c conda-forge python=3.11 \
    pandas numpy scipy statsmodels scikit-learn matplotlib seaborn \
    openpyxl pyyaml requests
conda activate proteomics-agent
```

Verify:

```bash
conda run -n proteomics-agent python -c "import pandas, scipy, statsmodels, \
    sklearn, matplotlib, seaborn, openpyxl, yaml, requests; print('OK')"
```

---

## Step-by-step setup for *your* data

### 1. Lay out a project directory

Copy this example folder somewhere outside the skill repo (so updates to the skill don't trample your work):

```bash
cp -r ~/.claude/skills/proteomics-agent/examples/PRT1296_end_to_end ~/proteomics_myproject
cd ~/proteomics_myproject
mkdir -p data
```

### 2. Prepare your input files

You need at minimum:

#### `data/proteins_long.csv` — long-format protein-level table

One row per (protein, sample, replicate). Required columns:

| Column | Type | Example |
|---|---|---|
| `Protein.Group` | str | `Q88GQ0` (UniProt accession; `;`-joined for protein groups OK) |
| `Protein` | str | `Katg` (gene symbol) |
| `Protein.Description` | str | `Catalase-peroxidase` |
| `Sample` | str | `SN_0725_80` (a *sample group* — strain × condition; not a single run) |
| `Replicate` | str | `R1`, `R2`, `R3` |
| `Counts_sum` | float | MaxLFQ-style intensity (linear scale, not log) |

First few rows look like:

```
Protein.Group,Protein,Protein.Description,Sample,Replicate,Counts_sum
Q88GQ0,Katg,Catalase-peroxidase,SN_0725_80,R1,8312345.0
Q88GQ0,Katg,Catalase-peroxidase,SN_0725_80,R2,7891234.0
Q88GQ0,Katg,Catalase-peroxidase,SN_0725_80,R3,8120987.0
Q88GQ0,Katg,Catalase-peroxidase,SN_0725_83,R1,6234567.0
...
```

#### `data/t_test_xlsx/` — directory of pairwise t-test xlsx files

The JBEI in-house DIA-NN post-processor produces one xlsx per pairwise comparison, with a sheet named `Full t-test output` containing:

| Column | Notes |
|---|---|
| `Protein`, `Protein.Group`, `Protein.Names`, `Protein.Description` | identifiers |
| `log2_mean_<A>`, `log2_mean_<B>` | group means (here A and B are sample IDs) |
| `log2_std_<A>`, `log2_std_<B>` | group SDs |
| `t-test_stat`, `p-value`, `p_adjusted(BH)` | stats |
| `log2_Fold_change_A/B` | log2(A/B) — sign matches contrast name |

Filename convention: `t-test_<contrast_name>_<timestamp>.xlsx`. Names should match the keys in `config/analysis.yaml:contrasts`.

#### Optional inputs

- `data/top3_proteins.csv` — Top3 quant (same schema as `proteins_long.csv` with `Top_3pep_counts_mean` instead of `Counts_sum`). Used by script `00` and `07` for absolute-abundance-ish numbers.
- `data/proteins_summary.csv` — per-group mean/std/CV summary (used as a sanity-check input by script `00`; auto-derivable if missing).

If you don't have the t-test xlsx files but you do have `proteins_long.csv`, the pipeline can still run script `03` (two-factor ANOVA) and most downstream steps. Scripts `04` (UpSet of pairwise t-test results) will be skipped.

### 3. Edit `config/analysis.yaml`

The shipped config encodes the PRT1296 study. Edit for your project:

```yaml
project:
  id: MYSTUDY                      # your experiment ID
  organism: <full Latin binomial>
  taxon_id: <NCBI taxid>           # 160488 for P. putida, 83333 for E. coli K-12, ...
  kegg_org: <3-letter KEGG code>   # ppu, eco, bsu, sce, hsa, etc.

samples:                           # sample-id → strain × medium × label
  SN_xxx_01: {strain: WT,  medium: ctrl,  label: "WT-ctrl"}
  SN_xxx_02: {strain: WT,  medium: stress, label: "WT-stress"}
  SN_xxx_03: {strain: KO,  medium: ctrl,  label: "KO-ctrl"}
  SN_xxx_04: {strain: KO,  medium: stress, label: "KO-stress"}

replicates: [R1, R2, R3]

contrasts:                         # one entry per pairwise t-test xlsx you have
  WT_stress_vs_ctrl:    [WT-stress, WT-ctrl]   # (A, B) → log2FC = mean(A) - mean(B)
  KO_stress_vs_ctrl:    [KO-stress, KO-ctrl]
  stress_KO_vs_WT:      [KO-stress, WT-stress]
  ctrl_KO_vs_WT:        [KO-ctrl,   WT-ctrl]

thresholds:
  q_strict: 0.05      # primary BH-FDR cutoff
  log2fc_strict: 1.0  # primary |log2FC| cutoff (2-fold)
  q_lenient: 0.10     # secondary, used for enrichment input
  log2fc_lenient: 0.585  # 1.5-fold

plot:
  dpi: 300
  palette_strain: {WT: "#4477AA", KO: "#EE6677"}
  palette_medium: {ctrl: "#BBBBBB", stress: "#228833"}
```

The shipped config is for a 2 (strain) × 3 (medium) design. The pipeline works for any number of factor levels — just add more entries.

### 4. Edit `config/tolerance_modules.yaml` (optional)

15 curated stress/metabolism modules tailored to *P. putida* biology (efflux, OMP porins, ROS, chaperones, choline/betaine, lysine catabolism, β-oxidation, ED pathway, etc.). For a different organism, swap the gene symbols. Each module:

```yaml
my_module_name:
  description: "What this module captures"
  genes: [gene1, gene2, gene3]       # exact gene symbols, case-insensitive
  gene_prefixes: [rps, rpl]          # optional — match any gene starting with these
```

Matching is done against the `Protein` field of `proteins_long.csv` and the UniProt `Gene Names` field of cached annotations.

### 5. Run

```bash
conda activate proteomics-agent
bash scripts/run_all.sh
```

Expected timing (Apple M-series or modern x86 laptop, ~1000 proteins × 18 runs):

| Stage | Time |
|---|---|
| `00_load_and_metadata.py` | 5 s |
| `01_qc.py` | 10 s |
| `02_exploratory.py` (per-protein ANOVA × 1021) | 30 s |
| `03_two_factor_anova.py` | 30 s |
| `04_dep_sets.py` | 15 s |
| `05_heatmap.py` | 10 s |
| `06_enrichment.py` first run (UniProt + KEGG fetch for ~1000 proteins) | **3–5 min** |
| `06_enrichment.py` cached re-runs | 30 s |
| `07_modules.py` | 5 s |
| `08_interaction_profiles.py` | 10 s |
| `09_extended_analysis.py` | 30 s |
| **Total first run** | **5–8 min** |
| **Total cached re-run** | **~2 min** |

### 6. Inspect outputs

After the run completes:

```bash
ls outputs/figures/      # ~17–20 PNGs; F1 PCA → F13 stress-amplitude
ls outputs/tables/       # ~27 TSVs; master DEP table, ANOVA, module scores, interactions
ls outputs/enrichment/   # ~60 TSVs; per-contrast × ontology × direction, plus cluster/set summaries
```

Key files to look at first:

| File | Use |
|---|---|
| `outputs/figures/F1_PCA.png` | Does PCA separate by your factors as expected? If not, batch-effect or labelling issue. |
| `outputs/figures/F2_QC_sample_correlation.png` | Within-group r should be > 0.90. |
| `outputs/figures/F8_heatmap_DEPs.png` | Headline figure for any presentation. |
| `outputs/figures/F10_module_heatmap.png` | Module-level summary across all groups. |
| `outputs/tables/two_factor_anova.tsv` | Search for proteins with `q_interaction <= 0.05`. |
| `outputs/enrichment/cluster_top_terms.tsv` | Pathway identity of each heatmap cluster. |

---

## Common gotchas

- **`HTTP 400` from UniProt** during script `06`: caused by `;`-joined protein-group accessions from DIA-NN. The `clean_accession()` helper splits on `;` and takes the first member. Already applied in this example.
- **Empty enrichment results**: usually a `kegg_org` mismatch. Make sure the 3-letter code in `analysis.yaml` matches a KEGG-indexed organism — run `curl https://rest.kegg.jp/list/organism | grep -i <your-organism>` to find it.
- **All proteins ns in two-factor ANOVA**: high replicate variance. Check `outputs/tables/qc_per_group_CV.tsv` — median CV should be < 35 %. If not, look at the QC PNGs for a bad run.
- **`Counts_sum` is already log2 in your data**: the pipeline expects *linear*-scale intensity and applies log2 itself. If yours is already log2, edit `scripts/_common.py:log2_safe` to be a pass-through, OR back-exponentiate your input column.
- **Different replicate naming** (e.g. `Rep1`, `Rep2` or `1`, `2`): edit `config/analysis.yaml:replicates` to match exactly.

## Design choices worth knowing

- **Background for ORA** = the actually-quantified proteins (whatever was in `proteins_long.csv`), *not* the full ORFeome. This corrects for detection-bias in single-shot DIA.
- **Enrichment goes through UniProt + KEGG REST** rather than gProfiler. Slower on first run but caches everything to `data/uniprot_annotation.tsv`, `data/kegg_pathway_map.tsv`, `data/go_ontology_map.tsv`. Works for any organism KEGG indexes.
- **Variance partitioning** uses Type-II SS via residual subtraction (`SS_strain = r_strain² − r_inter²`) so it's correct for balanced designs and the per-protein percentages sum to 100 %.
- **UpSet plotter** is self-contained (gridspec + scatter) — no dependency on `upsetplot`, which has matplotlib-version compatibility issues.
- **Interaction categorizer** rule order: (1) sign-flip → reciprocal, (2) amplitude difference → buffered / hyper-induced, (3) baseline shift → constitutive. This order avoids over-calling constitutive-shift when one strain has a much bigger stress response.

## Provenance

Built 2026-05-04 / audited 2026-05-05 on the PRT1296 study (*P. putida* KT2440 wild-type vs HGL1175 ALE-evolved strain across glucose / M9+cholinium-lysinate / M9+butylamine, 35 % v/v hydrolysate). Cross-validation against the existing JBEI t-test pipeline confirmed all 11 contrast DEP counts matched exactly. Variance partition verified to sum to 100 % per protein. Module gene matching audited — no false negatives.

See `data_analysis_plan.md` for the full 13-section plan including biological hypotheses, statistical thresholds, interpretation framework, and risks/caveats. The plan structure is data-agnostic; generalize §1 (background) to your organism and stressor and reuse the rest.
