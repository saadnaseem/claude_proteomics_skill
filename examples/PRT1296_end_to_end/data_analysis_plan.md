# Data Analysis Plan — PRT1296

**Project:** Proteomic and Metabolic Characterization of Wild-Type *P. putida* KT2440 and ALE-Evolved HGL1175 in Cholinium-Lysinate (Chlys) and Butylamine Hydrolysates
**Experiment ID:** PRT1296
**Author / PI lab:** Saad Naseem (LBL / JBEI)
**Plan drafted:** 2026-05-04
**Source data root:** `proteomics/ALE CHlys BAPRT1296_JBEI_20250820_SNaseem__with_full_ids.pr_matrix/`

---

## 1. Biological background and rationale

Lignocellulosic biomass pretreatment with ionic liquids and amines yields hydrolysates that are *fermentable but inhibitory*. Cholinium lysinate ([Ch][Lys]) and butylamine each leave residues with distinct toxicity modes:

| Pretreatment | Residual stressors | Expected stress biology |
|---|---|---|
| **Cholinium lysinate (Chlys)** | Cholinium cation (membrane intercalator, osmolyte mimic), free lysine, partially-degraded lignin phenolics, mild salt stress | Membrane fluidity perturbation; choline catabolism (BetA/BetB → glycine betaine); efflux of cation; lysine catabolism via cadaverine and δ-aminovalerate; phenolic detox via GST/SDR/oxidoreductases |
| **Butylamine hydrolysate** | Free butylamine (primary amine, protonophore, periplasmic oxidative substrate), residual lignin oligomers, possibly higher residual amine | Amine oxidase / periplasmic copper-dependent dehydrogenase (e.g., PedE/PedH/QedH); deamination feeding into central C; outer membrane porin remodelling; pH/charge stress; ROS from amine oxidation |
| **Glucose / M9 (control)** | none | Baseline central metabolism, ED pathway dominant in *P. putida* |

KT2440 already tolerates many lignocellulose-derived inhibitors (TtgABC, gst, ahpC, katG, oprH); HGL1175 is an ALE-evolved derivative selected for *enhanced* tolerance. Classic ALE outcomes in pseudomonads: (i) constitutive upregulation of efflux pumps, (ii) altered fatty-acid saturation/cyclopropanation, (iii) elevated ROS scavenging, (iv) rewired carbon entry, (v) loss-of-function in regulatory repressors (PA1226/TtgR-type) producing constitutive stress phenotypes.

The proteomic readout will tell us which of these mechanisms are operative, whether they are *constitutive* in HGL1175 (i.e., already on in glucose) versus *induced*, and which differ between Chlys and butylamine challenge.

---

## 2. Experimental design recap

### 2.1 Factor structure

A 2 (strain) × 3 (medium) factorial, 3 biological replicates per cell, 18 mass-spec runs total.

| | Glucose / M9 | M9 + Chlys | Butylamine hydrolysate |
|---|---|---|---|
| **KT2440** | SN_0725_80 (R1–R3) | SN_0725_83 (R1–R3) | SN_0725_86 (R1–R3) |
| **HGL1175** | SN_0725_89 (R1–R3) | SN_0725_92 (R1–R3) | SN_0725_95 (R1–R3) |

### 2.2 What's already in the pipeline (don't redo)

- DIA-NN protein quantification → `Full_list_proteins_20250902-195343.csv` (long format, 1021 proteins × 6 groups × 3 reps)
- Per-sample summaries with CV, z-score → `Full_list_proteins_summary_20250902-195343.csv`
- Top3 quantification → `Top3_quant_files/`
- PCA on 18 runs → `PCA_plot/`
- 11 pairwise t-tests with BH-adjusted p-values, volcano plots → `t-test_files/`, both p and q-value variants
- QC plots (per-run protein/peptide counts, CV violins, log2 boxplots) → `QC_files/`
- EDD-format upload table → `EDD_files/`
- Bar/strip charts for a curated pathway list → `Bar_Charts/`, `Strip_Charts/`

### 2.3 What's missing and worth building

1. Two-factor (strain × condition) linear model (currently only pairwise t-tests).
2. Functional enrichment (GO, KEGG, COG) per contrast.
3. UpSet / Venn intersection of DEP sets across contrasts.
4. Hierarchical-clustered z-score heatmap of significant proteins.
5. Curated "tolerance module" scoring (efflux, membrane, ROS, chaperone, amine/choline catabolism).
6. Strain-effect-in-glucose vs strain-effect-in-hydrolysate interaction analysis (the key ALE question).
7. Integration of proteomics with growth-curve / EDD metabolite measurements.

---

## 3. Hypotheses

### 3.1 Primary

**H1.** HGL1175 displays a *constitutive* stress-tolerance proteome — efflux pumps, ROS scavengers, membrane-remodelling enzymes are elevated relative to KT2440 *even in glucose*, indicating the ALE adaptation is regulatory, not inducible only in hydrolysate.

**H2.** Chlys and butylamine elicit *partially overlapping but mechanistically distinct* stress responses in KT2440. Chlys-specific signature is dominated by choline/glycine-betaine pathway proteins and cation efflux; butylamine-specific signature is dominated by periplasmic amine oxidation (PQQ-dependent dehydrogenases) and outer-membrane porin remodelling.

**H3.** HGL1175's tolerance advantage in hydrolysate manifests as *attenuated stress-response amplitude* (less induction needed because baseline is already elevated) plus *novel* hydrolysate-specific proteins not seen in KT2440.

### 3.2 Secondary

**H4.** Carbon entry routes shift between conditions: ED pathway (Edd, Eda, Zwf) dominant in glucose; β-oxidation, amino-acid catabolism, and GntZ/KGD shunts elevated in hydrolysates.

**H5.** A subset of DEPs in HGL1175 vs KT2440 (glucose) maps to genes annotated in the published HGL1175 ALE genome (mutations in regulatory loci, IS-element insertions). These give a falsifiable link between mutation and proteome.

**H6.** Translational and ribosomal proteins are *down* in hydrolysate-stressed KT2440 (general growth-rate reduction signature) but *less* down in HGL1175 (preserved growth).

---

## 4. Analysis questions (organized)

### 4.1 Quality / sanity
- Q1.1 — Do replicates cluster by condition in PCA? Are there outlier runs?
- Q1.2 — Are protein/peptide counts and CVs consistent across all 18 runs?
- Q1.3 — Is normalization adequate (boxplots, median-equalization)?

### 4.2 Strain effect (KT vs HGL1175)
- Q2.1 — In glucose only (lowest stress), what proteins differ between strains? This is the "ALE baseline shift."
- Q2.2 — Does the strain effect in glucose overlap with the strain effect in Chlys / butylamine? (i.e., interaction term: is HGL1175's advantage *strain-intrinsic* or *condition-revealed*?)
- Q2.3 — Are any proteins reciprocally regulated (up in HGL1175 only in hydrolysate, down in glucose, or vice versa)?

### 4.3 Condition effect (hydrolysate vs glucose)
- Q3.1 — What are the shared core hydrolysate-response proteins (KT and HGL1175 both regulate them in the same direction)?
- Q3.2 — Chlys-specific proteins (regulated in Chlys vs glucose but not in butylamine vs glucose).
- Q3.3 — Butylamine-specific proteins.
- Q3.4 — What metabolic pathways are reorganized? (Carbon entry, energy production, amino acid metab, fatty acid metab.)

### 4.4 Tolerance mechanisms (curated modules)
- Q4.1 — Efflux: TtgABC, TtgDEF, TtgGHI, MexEF, SrpABC, MdtABC. Levels in each cell?
- Q4.2 — Outer-membrane / porins: OprH, OprD, OprF, OprG, OmpA family.
- Q4.3 — Membrane lipid remodelling: cyclopropane fatty acid synthase (Cfa), FabA/FabB ratio, phospholipase, cardiolipin synthase.
- Q4.4 — ROS / oxidative stress: KatG, KatA, AhpC, AhpF, SodB, SodC, OxyR target genes, Gpx, Trx, Grx.
- Q4.5 — Chaperones / proteostasis: GroEL, GroES, DnaK, DnaJ, ClpB, ClpX, Lon, IbpA (HSP20).
- Q4.6 — General stress sigma factors: RpoS, RpoH (where detected).
- Q4.7 — Amine oxidation (butylamine-specific): PedE/PedH (PQQ-dependent), copper amine oxidase, BetA/BetB analogs.
- Q4.8 — Choline / glycine-betaine pathway (Chlys-specific): BetA, BetB, BetI, BetT.
- Q4.9 — Lysine catabolism (Chlys-specific): CadA, DavB, DavA, DavD, DavT (δ-aminovalerate route to glutarate).
- Q4.10 — Aromatic / phenolic detox: catechol meta-cleavage (XylE/CatA-like), GST, SDRs, FAD-monooxygenases.

### 4.5 Carbon and energy
- Q5.1 — ED pathway proteins (Zwf, Pgl, Edd, Eda) in each condition.
- Q5.2 — TCA cycle proteins.
- Q5.3 — Electron transport: cytochrome bo3 (Cyo), cbb3 oxidase, NADH dehydrogenases, F0F1-ATPase.
- Q5.4 — β-oxidation if present in HGL1175 (FadA, FadB, FadE).

### 4.6 Translation / growth
- Q6.1 — 30S/50S ribosomal protein abundance summed by condition. Surrogate for growth rate.
- Q6.2 — Translation factors (EF-Tu, EF-G, IF-3).

### 4.7 Integration
- Q7.1 — Do HGL1175 mutated loci (from genome resequencing) appear as DEPs?
- Q7.2 — Do growth-rate / OD600 differences track with ribosome/TCA proteome?
- Q7.3 — Do EDD-recorded extracellular metabolites (sugar uptake, organic acids) match up with the predicted metabolic shifts?

---

## 5. Step-by-step analysis breakdown

Each step lists: **input → method → output → success criterion**.

### Step 0 — Project bootstrap

- Create an analysis project directory under `proteomics/ALE_claude/` with subfolders: `data/`, `scripts/`, `outputs/figures/`, `outputs/tables/`, `outputs/enrichment/`, `notebooks/`, `memo/`.
- Symlink (or copy) the four authoritative inputs into `data/`:
  - `Full_list_proteins_20250902-195343.csv` (long-format counts)
  - `Top3_Full_list_proteins_20250902-195343.csv` (Top3 quant — preferred for absolute abundance)
  - `PCA-output_20250902-195343.csv`
  - All 11 `t-test_*.xlsx` files
- Write a `sample_metadata.tsv` with columns `sample_id, replicate, strain, medium, run_order, batch` — single source of truth for downstream code.

### Step 1 — QC and normalization sanity

**Input:** `Full_list_proteins_*.csv`, existing QC PNGs.
**Method:**
1. Read into pandas, pivot to a *protein × run* matrix (3 reps × 6 groups = 18 columns).
2. Verify median log2 intensities are aligned (they should be — DIA-NN MaxLFQ normalizes — but confirm).
3. Compute per-protein missing-value pattern. Plot fraction of proteins detected in (i) all 18 runs, (ii) all 3 reps of at least one group, (iii) only in one group (likely on/off proteins — important).
4. Pearson correlation heatmap of the 18 runs. Replicates within a group should be r>0.95.
5. CV per group; flag any group with median CV > 30 % (re-examine those runs).

**Output:**
- `outputs/figures/QC_corr_heatmap_18runs.png`
- `outputs/figures/QC_missingness_upset.png`
- `outputs/tables/per_group_CV.tsv`

**Success:** No outlier run; replicate correlation > 0.95; median group CV ≤ 25 %.

### Step 2 — Exploratory: PCA, hierarchical clustering, sample correlation

**Input:** protein × run matrix (log2-transformed, missing-value imputed by half-minimum or left as-is).
**Method:**
1. PCA on log2 intensities (already exists — but redo with explicit variance % and colour by strain *and* condition to see the dominant axis).
2. Hierarchical clustering of runs (Ward, 1-Pearson distance). Verify the dendrogram splits first by *strain* or by *condition* — this tells us which factor explains more variance.
3. Variance partitioning: fit a linear mixed model `~ strain + condition + strain:condition + (1|replicate)` per-protein, summarize variance contributions across the proteome.

**Output:**
- `outputs/figures/PCA_PC1PC2_strain_condition.png` (and PC1 vs PC3)
- `outputs/figures/dendrogram_18runs.png`
- `outputs/figures/variance_partition_violin.png`

**Success criterion:** Variance partitioning identifies which factor dominates; informs whether to discuss "strain-driven" vs "condition-driven" proteome.

### Step 3 — Differential abundance: re-use existing t-tests + add 2-way model

**Input:** the 11 pairwise t-test xlsx outputs.
**Re-use:**
- Load each `Full t-test output` sheet into a long-format DataFrame with columns `[contrast, protein, log2FC, p, q, t_stat]`.
- Apply consistent thresholds for downstream work: **|log2FC| ≥ 1 AND q ≤ 0.05** (primary stringent), and **|log2FC| ≥ 0.585 (1.5-fold) AND p ≤ 0.05** (lenient, for enrichment sensitivity).

**Add:** two-factor analysis (limma-style or `statsmodels` OLS on log2 protein).
- Per protein: `log2_intensity ~ strain + medium + strain:medium`
- Extract three F-tests per protein: main strain effect, main medium effect, interaction. BH-correct each.
- The **interaction term** is the most biologically interesting: proteins where the strain effect *depends* on medium = ALE-specific stress-response wiring.

**Output:**
- `outputs/tables/all_contrasts_long.tsv` (combined DEP table)
- `outputs/tables/two_factor_anova.tsv` (one row per protein, columns: F_strain, q_strain, F_medium, q_medium, F_interaction, q_interaction)
- `outputs/tables/interaction_significant_proteins.tsv` (q_interaction < 0.05)

### Step 4 — Set-based DEP analysis (Venn / UpSet)

**Input:** DEP lists from each contrast (q ≤ 0.05, |log2FC| ≥ 1).
**Method:**
1. Define the six "biology-direct" contrasts:
   - C1: KT_Chlys_vs_glucose
   - C2: KT_butylamine_vs_glucose
   - C3: HGL1175_Chlys_vs_glucose
   - C4: HGL1175_butylamine_vs_glucose
   - C5: glucose_HGL1175_vs_KT (strain baseline)
   - C6: Chlys_HGL1175_vs_KT, C7: butylamine_HGL1175_vs_KT (strain in stress)
2. UpSet plot of UP-regulated DEPs across C1–C4 (which proteins are core hydrolysate response vs medium-specific).
3. Separate UpSet for DOWN-regulated DEPs.
4. Strain-vs-strain UpSet of C5, C6, C7 — proteins differential between strains in *all* media = constitutive ALE shift; proteins differential only in stress = induced by stress and amplified by ALE.

**Output:**
- `outputs/figures/upset_hydrolysate_response_UP.png`
- `outputs/figures/upset_hydrolysate_response_DOWN.png`
- `outputs/figures/upset_strain_effect_across_media.png`
- `outputs/tables/dep_set_assignments.tsv`

**Success:** Clean partitioning into "core", "Chlys-only", "butylamine-only", "strain-baseline", "strain-amplified" categories.

### Step 5 — Heatmap of significant proteins

**Input:** union of proteins significant in any contrast (q ≤ 0.05, |log2FC| ≥ 1).
**Method:**
1. Z-score each protein across the 18 runs (or across 6 group means).
2. Hierarchical clustering on rows (proteins) and columns (groups).
3. Annotate column dendrogram with strain and medium colour bars.
4. Annotate row clusters with the dominant DEP-set assignment from Step 4.
5. Optionally: split heatmap into 4–6 row clusters and label them by inferred function (run enrichment per cluster, Step 6).

**Output:**
- `outputs/figures/heatmap_DEPs_clustered.png` (high-res, publication-grade)
- `outputs/tables/heatmap_cluster_membership.tsv`

### Step 6 — Functional enrichment

**Input:** DEP lists per contrast, plus per-cluster lists from Step 5.
**Method:**

1. Map UniProt IDs (e.g., Q88IC8) and gene symbols to:
   - GO terms (BP / MF / CC) via UniProt-GOA for *P. putida* KT2440 (taxon 160488).
   - KEGG pathways via `ppu` organism code.
   - COG categories via eggNOG-mapper or pre-annotated table.
2. Use `gprofiler-python` (organism `ppukt2440`) or `clusterProfiler`-equivalent for over-representation analysis (ORA) with the full quantified proteome (1021 proteins) as background — *not* the genome — to avoid bias from undetected proteins.
3. Run GSEA on the ranked log2FC list per contrast (sign-aware enrichment).
4. Multiple-testing within each ontology with BH.

**Output:**
- `outputs/enrichment/GO_BP_<contrast>.tsv` etc. for each contrast and each cluster
- `outputs/figures/enrichment_dotplot_<contrast>.png` (top 15 terms, dot size = count, colour = q)
- `outputs/figures/GSEA_enrichment_curves_top.png`

**Success:** Each contrast yields ≥ 5 significant pathways (q < 0.05), interpretable in light of hypotheses.

### Step 7 — Curated tolerance-module scoring

**Input:** Top3 quantification (absolute-ish abundance), curated gene lists.
**Method:**
1. Build a YAML/JSON file `config/tolerance_modules.yaml` with the gene lists from Section 4.4 (efflux, OMPs, membrane lipid, ROS, chaperone, amine oxidation, choline, lysine catabolism, aromatic detox, ribosome, TCA, ED).
2. For each module, compute (i) summed log2 intensity per sample, (ii) z-score each sample within module, (iii) plot heatmap (modules × 6 groups) or radar/spider plot (one polygon per group).
3. Test each module score for strain × medium effects with two-way ANOVA → BH-correct across modules.

**Output:**
- `outputs/figures/tolerance_module_heatmap.png`
- `outputs/figures/tolerance_module_radar.png`
- `outputs/tables/tolerance_module_anova.tsv`

**Why this matters:** moves the analysis from "list of proteins" to "tolerance phenotype score," easier to communicate.

### Step 8 — Strain × condition interaction deep-dive

**Input:** interaction-significant proteins from Step 3.
**Method:**
1. For each interaction-significant protein, plot the 4-point profile: (KT-glucose, KT-stress, HGL-glucose, HGL-stress) for Chlys and butylamine separately.
2. Categorize:
   - **Constitutively-elevated in HGL1175** (parallel curves, HGL > KT in all media)
   - **Hyper-induced in HGL1175** (HGL responds more strongly to stress)
   - **Buffered in HGL1175** (HGL responds less — stress is "ignored" because baseline already adequate)
   - **Reciprocal** (opposite direction)
3. Tabulate counts per category and per medium.

**Output:**
- `outputs/figures/interaction_profile_grid_top30.png`
- `outputs/figures/interaction_category_barplot.png`
- `outputs/tables/interaction_category_assignment.tsv`

### Step 9 — Pathway maps (visual)

**Input:** Per-protein log2FC for the four hydrolysate-vs-glucose contrasts.
**Method:**
1. Use `pathview`-style overlay on KEGG `ppu` central-carbon, oxidative-phosphorylation, amino-acid metabolism maps.
2. Colour each enzyme node by log2FC; render one map per contrast.
3. Custom-draw a stress-response cartoon (efflux + OMP + ROS + chaperone) with same colour scheme.

**Output:**
- `outputs/figures/pathway_KEGG_<map>_<contrast>.png` (e.g., `ppu00010_glycolysis_KT_Chlys.png`)
- `outputs/figures/pathway_stress_cartoon_<contrast>.svg` (editable)

### Step 10 — Integration with growth and EDD metabolite data

**Input:** EDD records (`EDD_files/EDDformat_*.csv`), and any growth-curve / OD600 data the user holds (request if absent).
**Method:**
1. Pull OD600 trajectories for the same 6 groups → estimate µ_max and lag.
2. Pull extracellular sugar / organic acid traces.
3. Correlate: ribosomal-protein module score vs µ_max (expect positive); ED-pathway score vs glucose uptake rate; amine-oxidase score vs butylamine disappearance rate.
4. Single integrative figure with growth, metabolite consumption, and proteome module scores stacked for the 6 groups.

**Output:**
- `outputs/figures/integration_growth_metabolome_proteome.png`
- `outputs/tables/integration_correlations.tsv`

### Step 11 — Cross-reference with HGL1175 genome (if resequencing exists)

**Input:** mutation list from HGL1175 vs KT2440 ALE genome resequencing (request from user).
**Method:** intersect mutated genes with DEPs from C5–C7. Flag concordant cases (mutation + protein change in same direction or in regulator → target).

**Output:** `outputs/tables/genotype_proteotype_concordance.tsv` and a small narrative `memo/genotype_proteotype.md`.

---

## 6. Plot inventory (final figure pack)

Numbered roughly in the order they would appear in a manuscript or report:

| # | Figure | Source step | Purpose |
|---|---|---|---|
| F1 | PCA (PC1 vs PC2) coloured by strain, shaped by medium | Step 2 | Show clean group separation, dominant variance axis |
| F2 | Sample correlation heatmap, 18 runs | Step 1 | QC, replicate consistency |
| F3 | Variance-partition violin (per-protein strain / medium / interaction) | Step 2 | Quantify what dominates the proteome |
| F4 | Volcano plots, 4-panel: KT-Chlys/glu, KT-buty/glu, HGL-Chlys/glu, HGL-buty/glu | existing | Per-condition response in each strain |
| F5 | Volcano plots, 3-panel: HGL vs KT in glucose / Chlys / butylamine | existing | Strain effect in each medium |
| F6 | UpSet of UP DEPs across 4 stress contrasts | Step 4 | Core vs medium-specific hydrolysate response |
| F7 | UpSet of DEPs across 3 strain contrasts | Step 4 | Constitutive vs induced ALE shift |
| F8 | Z-score heatmap of all DEPs, hierarchically clustered, with row-cluster annotations | Step 5 | Global proteome view |
| F9 | Enrichment dotplots, one per contrast (top 15 terms) | Step 6 | Pathway-level interpretation |
| F10 | Tolerance-module heatmap (modules × groups) | Step 7 | Stress-mechanism summary |
| F11 | Tolerance-module radar plot | Step 7 | Visual phenotype comparison |
| F12 | Interaction-protein profile grid (top 20–30) | Step 8 | The headline ALE-specific biology |
| F13 | KEGG pathway overlays (central carbon + amino-acid catabolism) | Step 9 | Metabolic rewiring picture |
| F14 | Custom stress-response cartoon, fold-change coloured | Step 9 | Communication figure for talks |
| F15 | Integrated growth × metabolite × proteome panel | Step 10 | Phenotype-to-mechanism story |
| F16 | Bar plots for curated gene panels (already drafted in `Bar_Charts/`) | existing | Targeted illustrations of key genes |

Supplementary figures: per-replicate boxplots, missing-value UpSet, full enrichment tables, all individual interaction-profile plots.

---

## 7. Output deliverables

1. **`outputs/tables/master_DEP_table.tsv`** — one row per protein × contrast, with all stats and category labels (DEP set, interaction category, module membership).
2. **`outputs/tables/two_factor_anova.tsv`** — per-protein strain / medium / interaction F-stats and BH q-values.
3. **`outputs/tables/tolerance_module_scores.tsv`** — per-sample score per module.
4. **`outputs/enrichment/`** — full ORA / GSEA tables for every contrast, every ontology.
5. **All figures** in `outputs/figures/` as PNG + SVG (editable).
6. **`memo/results_narrative.md`** — 3–5 page interpretation document, hypothesis-by-hypothesis verdict.
7. **`memo/figure_legends.md`** — manuscript-ready legends for F1–F16.

---

## 8. Implementation: scripts to write

| Script | Purpose | Key libs |
|---|---|---|
| `scripts/00_load_and_metadata.py` | Build sample_metadata.tsv, load long protein file, pivot to wide | pandas |
| `scripts/01_qc.py` | Correlation heatmap, missingness UpSet, CV table | pandas, seaborn, upsetplot |
| `scripts/02_exploratory.py` | PCA redo, dendrogram, variance partition | scikit-learn, scipy, statsmodels |
| `scripts/03_two_factor_model.py` | Per-protein OLS, ANOVA, interaction extraction, BH | statsmodels, pingouin |
| `scripts/04_dep_sets.py` | Aggregate t-test xlsx, UpSet, set assignments | upsetplot |
| `scripts/05_heatmap.py` | Z-score, cluster, plot, export cluster membership | scipy, seaborn |
| `scripts/06_enrichment.py` | gProfiler ORA + GSEApy on ranked log2FC | gprofiler-python, gseapy |
| `scripts/07_modules.py` | Curated module scores, ANOVA, heatmap, radar | pandas, matplotlib |
| `scripts/08_interaction_profiles.py` | Categorize interaction-significant proteins, profile grid | pandas, matplotlib |
| `scripts/09_pathway_maps.py` | KEGG pathway overlay via Bioservices/`pathview` (R) or pypath | bioservices, pypath |
| `scripts/10_integration_eddi.py` | Pull EDD measurements, correlate with proteome | pandas |
| `scripts/run_all.sh` | Reproducible end-to-end pipeline | bash |

A single `config/analysis.yaml` should hold thresholds (q-cutoff, log2FC cutoff), file paths, and module gene lists, so the pipeline is reproducible without code edits.

Recommended Python environment (pin in `environment.yml`):
`python>=3.11, pandas, numpy, scipy, statsmodels, scikit-learn, seaborn, matplotlib, upsetplot, gprofiler-official, gseapy, openpyxl, pyyaml, pingouin, bioservices`.

---

## 9. Statistical thresholds and reporting standards

- **Primary DEP cutoff:** BH-adjusted p ≤ 0.05 AND |log2FC| ≥ 1 (2-fold).
- **Lenient (enrichment input only):** raw p ≤ 0.05 AND |log2FC| ≥ 0.585 (1.5-fold).
- **Background for enrichment:** the 1021-protein quantified set, *not* the full KT2440 proteome (avoids detection bias).
- **Multiple testing across pathways:** BH per ontology (BP, MF, CC, KEGG separately).
- **Two-factor model:** report q_interaction; treat q_interaction ≤ 0.1 as suggestive (for an n=3 design, FDR-stringent interaction is hard).
- **Reporting:** every DEP in the narrative must cite the contrast, log2FC, q. Sentences like "X is upregulated" without numbers are not acceptable.

---

## 10. Interpretation framework — how each result maps to a hypothesis

| Observation | Conclusion about |
|---|---|
| Many DEPs in C5 (HGL vs KT in glucose) | H1: ALE caused a constitutive proteome shift |
| C5 DEPs largely overlap C6, C7 | H1 confirmed: shift is medium-independent |
| C6, C7 DEPs ⊃ C5 (extra DEPs in stress) | H3: ALE adds an induced layer on top of constitutive shift |
| Significant interaction (Step 3) for stress-response genes | H3 confirmed: ALE specifically rewires stress amplitude |
| Chlys-specific DEPs enriched in choline / glycine-betaine / lysine catabolism | H2 Chlys arm |
| Butylamine-specific DEPs enriched in PQQ / amine oxidase / OMP | H2 butylamine arm |
| Ribosomal module score: KT-stress < KT-glucose; HGL-stress ≈ HGL-glucose | H6 |
| Efflux module score elevated constitutively in HGL1175 | classic ALE outcome, supports H1/H3 |

Anti-hypothesis outcomes worth being explicit about:
- If C5 has very few DEPs but C6/C7 have many → ALE is *condition-revealed*, not constitutive (rejects H1).
- If Chlys and butylamine responses overlap > 80 % → toxicity is dominated by a shared component (lignin phenolics?), not by the named amine — rethink "Chlys-specific" framing.
- If HGL1175 responses are *larger* (not smaller) than KT in stress → tolerance is from sensing-and-responding faster, not from buffering.

---

## 11. Risks, caveats, and open questions

1. **n = 3 per cell** is the floor for statistics; interaction-term power is low. A "non-significant interaction" should not be over-interpreted.
2. **Coverage at 1021 proteins** (~19 % of the KT2440 ORFeome) means the most lowly-expressed regulators (RpoS, sigma factors, two-component sensors) may be missed. Cross-check with transcriptomics if available.
3. **Top3 quant** is semi-quantitative; absolute concentrations require iBAQ or label-mediated calibration. Statements about "absolute abundance" should be hedged.
4. **Hydrolysate composition variability** between batches — confirm the same Chlys and butylamine batches were used across the 6 cultures. If not, document it as a confound.
5. **Growth-rate confound**: DEPs in stress vs glucose conflate stress-response biology with reduced-growth-rate biology. The ribosome-module subtraction in Step 7 is a partial mitigation but not a full fix. To deconfound rigorously, would need matched-µ chemostat data.
6. **HGL1175 ALE history** — tolerance was selected on a specific medium; performance on a *different* hydrolysate (e.g., butylamine if it was selected on Chlys) is a generalization test, and a partial response is not a failure of ALE but an expected scope limit.

---

## 12. Suggested execution order and rough timing

Day 1 — Steps 0, 1, 2 (bootstrap, QC, exploratory).
Day 2 — Steps 3, 4 (two-factor model, DEP sets).
Day 3 — Steps 5, 6 (heatmap, enrichment).
Day 4 — Steps 7, 8 (modules, interaction deep-dive).
Day 5 — Step 9 (pathway maps); start narrative.
Day 6 — Step 10, 11 (integration, genotype overlay).
Day 7 — Polish figures, finish narrative memo, internal review.

---

## 13. Open items / questions for the user before starting

1. **Growth data:** are OD600 / µ_max / lag-time measurements available for the same 6 cultures? Critical for Step 10.
2. **EDD pulls:** which EDD study ID holds the metabolite time-courses? Confirm before Step 10.
3. **HGL1175 genome / mutation list:** is there a resequencing dataset for Step 11? If not, drop Step 11.
4. **Hydrolysate compositional data:** any HPLC characterization (sugar, residual amine) per batch?
5. **Sample SN_0725_xx → Strain×Medium mapping** confirmed via t-test filenames — does the user have an explicit metadata sheet to verify?
6. **Reporting target:** internal report, manuscript, conference talk? Drives figure polish vs analysis depth tradeoff.
7. **Significance thresholds:** confirm q ≤ 0.05, |log2FC| ≥ 1 are acceptable, or relax for an underpowered design (q ≤ 0.1)?

---

*End of plan.*
