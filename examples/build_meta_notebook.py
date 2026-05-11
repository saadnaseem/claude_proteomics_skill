"""Build a fully-populated Jupyter notebook that reproduces the entire Phase 1 + Phase 3 analysis."""
import nbformat as nbf
from pathlib import Path

NB = Path("/Users/snaseem/proteomics_runs/phase3_meta_analysis_20260510_212346/analysis.ipynb")
DATA_DIR = "/Users/snaseem/Library/CloudStorage/GoogleDrive-snaseem@lbl.gov/My Drive/Saad m-group/proteomics/ALE_claude/ALE CHlys BAPRT1296_JBEI_20250820_SNaseem__with_full_ids.pr_matrix/t-test_files/"
META_DIR = "/Users/snaseem/proteomics_runs/phase3_meta_analysis_20260510_212346"

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"name": "proteomics-agent", "display_name": "Python (proteomics-agent)", "language": "python"},
    "language_info": {"name": "python"},
}

cells = []

# ---- Title + intro ----
cells.append(nbf.v4.new_markdown_cell(f"""# Comparative proteomics: HGL1175 (ALE) vs KT2440 (WT) *Pseudomonas putida* across 3 media

**End-to-end reproducible analysis** — load all 11 t-test xlsx files, build master DE matrix, generate the buffered-response scatters, regulator co-expression heatmap, and tolerance candidate ranking.

This notebook reproduces the analysis described in:
- `STORY.md` — unified Phase 1 + Phase 3 narrative
- `STORY_groupB_addendum.md` — Group B (stress-response) detailed analysis

## Inputs
- 11 t-test xlsx files at `{DATA_DIR}` (sample IDs: SN_0725_80 KT/glu, _83 KT/Chlys, _86 KT/buta, _89 HGL/glu, _92 HGL/Chlys, _95 HGL/buta)

## Outputs (re-generated when this notebook is run)
- `tables/01_master_de_matrix.csv` — 1021 proteins × 9 comparisons
- `tables/02_constitutive_ALE_hits.csv`
- `tables/03_buffered_*.csv` — per-protein classification under each hydrolysate
- `tables/04_regulators_across_comparisons.csv`
- `tables/05_tolerance_candidates.csv`
- `figures/03_buffered_response_*.png`
- `figures/04_regulator_heatmap.png`

## Kernel
Run this with `Python (proteomics-agent)` kernel (already registered when the skill was installed).
"""))

# ---- Setup ----
cells.append(nbf.v4.new_markdown_cell("## Setup — imports, paths, plot defaults"))
cells.append(nbf.v4.new_code_cell(f'''import os, json, warnings
from pathlib import Path
from collections import defaultdict
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests, requests_cache
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
np.random.seed(42)

DATA_DIR = Path(r"{DATA_DIR}")
RUN_DIR  = Path(r"{META_DIR}")
FIG_DIR  = RUN_DIR / "figures"; FIG_DIR.mkdir(exist_ok=True)
TBL_DIR  = RUN_DIR / "tables"; TBL_DIR.mkdir(exist_ok=True)
ANN_DIR  = RUN_DIR / "annotations"; ANN_DIR.mkdir(exist_ok=True)

# Cache external API calls (UniProt, KEGG, STRING) for offline re-runs
requests_cache.install_cache(str(ANN_DIR / "api_cache"), backend="sqlite",
                              expire_after=86400 * 30)

plt.rcParams.update({{"figure.dpi": 100, "savefig.dpi": 200, "savefig.bbox": "tight",
                     "font.size": 10, "axes.spines.top": False, "axes.spines.right": False}})
COLOR_UP, COLOR_DOWN, COLOR_NS = "#D62728", "#1F77B4", "#AAAAAA"
print(f"Setup OK. RUN_DIR={{RUN_DIR}}")'''))

# ---- Section 1: load all comparisons ----
cells.append(nbf.v4.new_markdown_cell("""## 1. Load all 9 unique pairwise comparisons → master DE matrix

Each xlsx file has 5 sheets; we read only the `Full t-test output` sheet.

The comparisons are organized in three groups:
- **Group A** (strain effect): HGL1175 vs KT2440 under each medium — 3 files
- **Group B** (stress effect): each strain's response to each hydrolysate vs glucose — 4 files
- **Group C** (between hydrolysates): Chlys vs butamine within each strain — 2 files (the other 2 are sign-flipped duplicates)
"""))
cells.append(nbf.v4.new_code_cell('''SAMPLE = {
    "SN_0725_80": "KT/glu",  "SN_0725_83": "KT/Chlys", "SN_0725_86": "KT/buta",
    "SN_0725_89": "HGL/glu", "SN_0725_92": "HGL/Chlys","SN_0725_95": "HGL/buta",
}
COMPARISONS = {
    "A_glu_HGLvsKT":       ("t-test_glucose_HGL1175_vs_KT_20250902-195343.xlsx",  "SN_0725_89","SN_0725_80","A"),
    "A_Chlys_HGLvsKT":     ("t-test_Chlys_HGL1175_vs_KT_20250902-195343.xlsx",    "SN_0725_92","SN_0725_83","A"),
    "A_buta_HGLvsKT":      ("t-test_Butamine_HGL1175_vs_KT_20250902-195343.xlsx", "SN_0725_95","SN_0725_86","A"),
    "B_KT_Chlys_vs_glu":   ("t-test_KT_Chlys_vs_glucose_20250902-195343.xlsx",    "SN_0725_83","SN_0725_80","B"),
    "B_KT_buta_vs_glu":    ("t-test_KT_butamine_vs_glucose_20250902-195343.xlsx", "SN_0725_86","SN_0725_80","B"),
    "B_HGL_Chlys_vs_glu":  ("t-test_HGL1175_Chlys_vs_glucose_20250902-195343.xlsx","SN_0725_92","SN_0725_89","B"),
    "B_HGL_buta_vs_glu":   ("t-test_HGL1175_butamine_vs_glucose_20250902-195343.xlsx","SN_0725_95","SN_0725_89","B"),
    "C_KT_Chlys_vs_buta":  ("t-test_KT_Chlys_vs_butamine_20250902-195343.xlsx",   "SN_0725_83","SN_0725_86","C"),
    "C_HGL_Chlys_vs_buta": ("t-test_HGL1175_Chlys_vs_butamine_20250902-195343.xlsx","SN_0725_92","SN_0725_95","C"),
}

master = None
meta_dict = {}
for key, (fname, a_id, b_id, group) in COMPARISONS.items():
    df = pd.read_excel(DATA_DIR / fname, sheet_name="Full t-test output")
    df["acc"] = df["Protein.Group"].apply(lambda s: str(s).split(";")[0].strip())
    # Validate sign convention
    err = (df["log2_Fold_change_A/B"] - (df[f"log2_mean_{a_id}"] - df[f"log2_mean_{b_id}"])).abs().max()
    assert err < 0.01, f"Sign convention {key}: {err}"

    for _, r in df.iterrows():
        a = r["acc"]
        if a not in meta_dict:
            meta_dict[a] = {"protein": r["Protein"], "gene": r["Protein.Names"], "desc": r["Protein.Description"]}

    sub = df[["acc"]].copy()
    sub[f"{key}_log2fc"] = df["log2_Fold_change_A/B"]
    sub[f"{key}_padj"]   = df["p_adjusted(BH)"]
    master = sub if master is None else master.merge(sub, on="acc", how="outer")

    n_sig = ((df["p_adjusted(BH)"]<0.05) & (df["log2_Fold_change_A/B"].abs()>1)).sum()
    print(f"  {key:<22} A={SAMPLE[a_id]:<10} B={SAMPLE[b_id]:<10} group={group}  sig={n_sig}")

master["protein"] = master["acc"].map(lambda a: meta_dict.get(a,{}).get("protein",""))
master["gene"]    = master["acc"].map(lambda a: meta_dict.get(a,{}).get("gene",""))
master["desc"]    = master["acc"].map(lambda a: meta_dict.get(a,{}).get("desc",""))
master = master[["acc","protein","gene","desc"] + [c for c in master.columns if c not in ("acc","protein","gene","desc")]]
master.to_csv(TBL_DIR / "01_master_de_matrix.csv", index=False)
print(f"\\nMaster matrix: {master.shape}")
master.head()'''))

cells.append(nbf.v4.new_markdown_cell("""**Headline numbers** — note the asymmetry:
- WT (KT) mounts 17–27 sig changes when transferred to either hydrolysate
- ALE (HGL) mounts 0–3
- Strain difference under Chlys is only 4 hits (vs 22–29 under glucose/butamine) — convergence."""))

# ---- Section 2: Constitutive vs conditional ----
cells.append(nbf.v4.new_markdown_cell("""## 2. Constitutive vs conditional ALE signature

Decompose the Group A strain comparisons. Which proteins are differential in *all* three media (constitutive ALE signature) vs only some (conditional)?"""))
cells.append(nbf.v4.new_code_cell('''def sig_strict(prefix, p=0.05, fc=1.0):
    return (master[f"{prefix}_padj"]<p) & (master[f"{prefix}_log2fc"].abs()>fc)

in_glu  = sig_strict("A_glu_HGLvsKT")
in_buta = sig_strict("A_buta_HGLvsKT")
in_chl  = sig_strict("A_Chlys_HGLvsKT")

print(f"All-3 (constitutive, strict):     {(in_glu & in_buta & in_chl).sum()}")
print(f"glu+buta only:                    {(in_glu & in_buta & ~in_chl).sum()}")
print(f"glu only:                          {(in_glu & ~in_buta & ~in_chl).sum()}")
print(f"buta only:                         {(~in_glu & in_buta & ~in_chl).sum()}")
print(f"chlys only:                        {(~in_glu & ~in_buta & in_chl).sum()}")
print(f"Multi-sig (≥2 of 3):              {((in_glu.astype(int)+in_buta.astype(int)+in_chl.astype(int)) >= 2).sum()}")

multisig = master[in_glu | in_buta | in_chl].copy()
multisig["n_groupA_sig"] = (in_glu.astype(int)+in_buta.astype(int)+in_chl.astype(int))[multisig.index]
constitutive = multisig[multisig["n_groupA_sig"]>=2].sort_values("n_groupA_sig", ascending=False)
constitutive[["acc","gene","desc","n_groupA_sig",
              "A_glu_HGLvsKT_log2fc","A_buta_HGLvsKT_log2fc","A_Chlys_HGLvsKT_log2fc"]].to_csv(
    TBL_DIR / "02_constitutive_ALE_hits.csv", index=False)
constitutive[["acc","gene","desc","n_groupA_sig",
              "A_glu_HGLvsKT_log2fc","A_buta_HGLvsKT_log2fc","A_Chlys_HGLvsKT_log2fc"]]'''))

cells.append(nbf.v4.new_markdown_cell("""**Constitutive set (significant in ≥2 Group A comparisons):**
- **PhaG (PP_5008)** + **PhaF GA2 (PP_5007)** — PHA storage DOWN, all 3 media (most consistent finding)
- **Idh (NADP-isocitrate DH)** — UP, glu+buta
- **KatG (catalase-peroxidase)** — DOWN, glu+buta
- **YeaG (stationary-phase kinase)** — DOWN, glu+buta

The Chlys comparison adds little to the strict-cutoff list because of low power (only 4 sig hits total under Chlys), but the regulator heatmap below shows the patterns continue at sub-threshold magnitudes."""))

# ---- Section 3: Buffered response ----
cells.append(nbf.v4.new_markdown_cell("""## 3. Buffered-response analysis (Group B)

For each hydrolysate, plot:
- **x axis**: WT (KT2440) log2FC under stress (vs glucose)
- **y axis**: ALE (HGL1175) log2FC under stress (vs glucose)

Proteins on the **diagonal** = same response in both strains (general response to the substrate).
Proteins **off the diagonal toward the WT axis** = WT responds, ALE doesn't = "buffered" by ALE.

⚠️ **Important caveat (per `STORY_groupB_addendum.md`)**: many "buffered" classifications under Chlys are statistical artifacts — both strains activate the cholinium catabolic pathway at similar effect size, but HGL's p-value falls just above the 0.10 threshold for some genes. Under butamine, the buffered classification holds better (true magnitude differences)."""))
cells.append(nbf.v4.new_code_cell('''SIG_P = 0.10
buffered_results = {}
for hyd, kt_key, hgl_key in [("Chlys","B_KT_Chlys_vs_glu","B_HGL_Chlys_vs_glu"),
                              ("butamine","B_KT_buta_vs_glu","B_HGL_buta_vs_glu")]:
    sub = master[["acc","gene","desc",
                  f"{kt_key}_log2fc",f"{kt_key}_padj",
                  f"{hgl_key}_log2fc",f"{hgl_key}_padj"]].dropna()
    sub.columns = ["acc","gene","desc","kt_log2fc","kt_padj","hgl_log2fc","hgl_padj"]
    sub["category"] = "ns"
    sub.loc[(sub.kt_padj<SIG_P) & (sub.hgl_padj>=SIG_P), "category"] = "buffered_by_ALE"
    sub.loc[(sub.kt_padj>=SIG_P) & (sub.hgl_padj<SIG_P), "category"] = "ALE_specific"
    sub.loc[(sub.kt_padj<SIG_P) & (sub.hgl_padj<SIG_P), "category"] = "shared_response"
    sub.to_csv(TBL_DIR / f"03_buffered_{hyd}.csv", index=False)
    buffered_results[hyd] = sub
    print(f"{hyd}: {sub['category'].value_counts().to_dict()}")

fig, axes = plt.subplots(1, 2, figsize=(16, 7.5))
for ax, (hyd, sub) in zip(axes, buffered_results.items()):
    cmap = {"ns":"#CCCCCC","buffered_by_ALE":"#D62728","ALE_specific":"#9467BD","shared_response":"#2CA02C"}
    for cat, color in cmap.items():
        s = sub[sub.category==cat]
        ax.scatter(s["kt_log2fc"], s["hgl_log2fc"],
                   s=12 if cat=="ns" else 35,
                   c=color, alpha=0.25 if cat=="ns" else 0.8,
                   edgecolors="none", label=f"{cat} (n={len(s)})")
    lim = max(sub["kt_log2fc"].abs().max(), sub["hgl_log2fc"].abs().max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], color="grey", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(0, color="black", lw=0.5); ax.axvline(0, color="black", lw=0.5)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel(f"WT KT2440: log2FC ({hyd} vs glucose)")
    ax.set_ylabel(f"ALE HGL1175: log2FC ({hyd} vs glucose)")
    ax.set_title(f"{hyd} stress")
    ax.legend(loc="upper left", fontsize=9)

    # Label top buffered hits
    buf = sub[sub.category=="buffered_by_ALE"].copy()
    buf["abs"] = buf["kt_log2fc"].abs()
    for _, r in buf.nlargest(10, "abs").iterrows():
        ax.annotate(str(r.gene) if pd.notna(r.gene) and str(r.gene) else r.acc,
                    (r.kt_log2fc, r.hgl_log2fc), fontsize=7,
                    xytext=(3,3), textcoords="offset points")
plt.suptitle("Buffered-response: WT vs ALE under hydrolysate stress\\n"
             "Off-diagonal red = ALE doesn\\'t need to do what WT does", fontsize=12, y=1.0)
plt.tight_layout()
plt.savefig(FIG_DIR / "03_buffered_response_combined.png")
plt.show()'''))

# ---- Section 4: Regulator heatmap ----
cells.append(nbf.v4.new_markdown_cell("""## 4. Regulator co-expression heatmap (THE KEY FIGURE)

Direct view of the **multi-regulator suppression** model. TurA, TurB, RelA, CsrA, KatG should all be blue (DOWN in HGL) across the Group A columns. Idh should be red (UP). PhaG/F always blue."""))
cells.append(nbf.v4.new_code_cell('''regulators = {
    "TurA (PP_1366, MvaT P16)":   "Q88N50",
    "TurB (PP_3765, MvaT H-NS)":  "Q88GF9",
    "RelA (ppGpp synthase)":      "Q88MB8",
    "CsrA (RNA-binding repr.)":   "Q88G93",
    "KatG (catalase-peroxidase)": "Q88GQ0",
    "GroES (chaperonin)":         "Q88N56",
    "DnaK (chaperone)":           "Q88DU2",
    "Fur (iron uptake reg.)":     "Q88CR2",
    "PhaG (PHA granule)":         "Q88D20",
    "PhaF GA2 (PHA granule)":     "Q88D21",
    "Pp_4834 (SPFH/Band 7)":      "Q88DJ1",
    "OmpA family Pp_1502":        "Q88MR7",
    "Cyoc (cyt bo3 oxidase)":     "Q88PN5",
    "Idh (NADP-IDH, TCA)":        "Q88FS1",
    "Pp_0241 (myst transporter)": "Q88R92",
}

reg_rows = []
for name, acc in regulators.items():
    if acc in master["acc"].values:
        row = master[master["acc"]==acc].iloc[0]
        rdict = {"name":name,"acc":acc,"gene":row["gene"]}
        for k in COMPARISONS:
            rdict[k] = row[f"{k}_log2fc"]
            rdict[f"{k}_p"] = row[f"{k}_padj"]
        reg_rows.append(rdict)
reg_df = pd.DataFrame(reg_rows)
reg_df.to_csv(TBL_DIR / "04_regulators_across_comparisons.csv", index=False)

fig, ax = plt.subplots(figsize=(14, 9))
log2fc_cols = list(COMPARISONS.keys())
heat = reg_df[log2fc_cols].copy(); heat.index = reg_df["name"]
heat.columns = [c.replace("_log2fc","").replace("HGLvsKT","HGL/KT").replace("_vs_","→") for c in log2fc_cols]
annot = pd.DataFrame(index=heat.index, columns=heat.columns, dtype=str)
for i, name in enumerate(heat.index):
    for j, k in enumerate(log2fc_cols):
        val = reg_df.iloc[i][k]; p = reg_df.iloc[i][f"{k}_p"]
        s = ""
        if pd.notna(p):
            if p<0.001: s="***"
            elif p<0.01: s="**"
            elif p<0.05: s="*"
            elif p<0.10: s="·"
        annot.iloc[i,j] = f"{val:+.1f}{s}" if pd.notna(val) else ""
sns.heatmap(heat, annot=annot, fmt="", cmap="RdBu_r", center=0, vmin=-3, vmax=3,
            ax=ax, cbar_kws={"label":"log2FC"}, linewidths=0.3, linecolor="white")
ax.set_title("Key regulators + effectors: log2FC across all 9 comparisons\\n"
             "Look at the first 3 columns (Group A): TurA, TurB, RelA, CsrA, KatG all blue (DOWN in HGL); Idh red (UP)\\n"
             "*** p_adj<0.001  ** p_adj<0.01  * p_adj<0.05  · p_adj<0.10",
             fontsize=10)
plt.tight_layout()
plt.savefig(FIG_DIR / "04_regulator_heatmap.png")
plt.show()'''))

cells.append(nbf.v4.new_markdown_cell("""**The pattern to see**:
- Look at the first 3 columns (`A_glu_HGLvsKT`, `A_Chlys_HGLvsKT`, `A_buta_HGLvsKT`)
- All 4 growth-suppressor regulators (TurA, TurB, RelA, CsrA) are **blue (negative)** across all 3 columns
- KatG (catalase) is also blue (constitutively low stress baseline)
- PhaG, PhaF (storage) are blue
- Idh (TCA) is red (constitutively up)
- This pattern, especially the simultaneous reduction of multiple growth-suppressor regulators, is the molecular substrate of the constitutive growth-optimized phenotype."""))

# ---- Section 5: Tolerance candidate ranking ----
cells.append(nbf.v4.new_markdown_cell("""## 5. Tolerance candidate ranking

Composite score = `|HGL_vs_KT under stress|` × `|KT_stress_vs_glucose|`, averaged across the two hydrolysates. High score = ALE has changed AND WT cares about it under stress."""))
cells.append(nbf.v4.new_code_cell('''scoring = master[["acc","gene","desc"]].copy()
for hyd, a_key, b_key in [("Chlys","A_Chlys_HGLvsKT","B_KT_Chlys_vs_glu"),
                          ("buta","A_buta_HGLvsKT","B_KT_buta_vs_glu")]:
    a_lfc = master[f"{a_key}_log2fc"].abs().fillna(0)
    b_lfc = master[f"{b_key}_log2fc"].abs().fillna(0)
    scoring[f"{hyd}_strain_eff"] = master[f"{a_key}_log2fc"]
    scoring[f"{hyd}_WT_resp"] = master[f"{b_key}_log2fc"]
    scoring[f"{hyd}_score"] = a_lfc * b_lfc
scoring["composite_score"] = (scoring["Chlys_score"] + scoring["buta_score"]) / 2
scoring["glu_strain_eff"] = master["A_glu_HGLvsKT_log2fc"]
top = scoring.nlargest(20, "composite_score")
top.to_csv(TBL_DIR / "05_tolerance_candidates.csv", index=False)
top[["acc","gene","desc","composite_score","glu_strain_eff",
     "Chlys_strain_eff","Chlys_WT_resp","buta_strain_eff","buta_WT_resp"]]'''))

# ---- Section 6: Take-home ----
cells.append(nbf.v4.new_markdown_cell("""## 6. Take-home

**Three-layer model of ALE-acquired hydrolysate tolerance** (corrected via Group B):

1. **Layer 1 — constitutive growth program**: HGL1175 has constitutively reduced TurA, TurB, RelA, CsrA, KatG, PhaG, PhaF and elevated Idh, regardless of medium. Multi-regulator suppression releases brakes on growth.

2. **Layer 2 — substrate utilization** (induced equally in both strains): the cholinium catabolic pathway (BetA, BetB, sarcosine ox., etc.) is activated by both KT and HGL when cholinium is present. Not a tolerance mechanism — both strains have it.

3. **Layer 3 — generic stress response damping** (the actual tolerance mechanism): KT additionally induces Cyoc (high-O₂ respiration), OmpA (membrane remodeling), sugar-acid + aromatic-AA catabolism (under butamine), and catalase. **HGL skips most of this layer** because Layer 1 already covers the underlying need.

**Top engineering targets**:
- TurA + TurB + RelA + CsrA combined knockout in WT (test multi-regulator hypothesis directly)
- Cyoc-induction prevention (test Layer 3 damping)
- Sequence HGL1175 to find the cis/trans mutations driving regulator suppression

**Outputs from this notebook are reproducible** — all generated tables and figures are saved to `tables/` and `figures/` folders. External API calls (UniProt, KEGG, STRING) are cached in `annotations/api_cache.sqlite`.

For the narrative versions of these findings:
- `STORY.md` — Phase 1 + Phase 3 unified analysis
- `STORY_groupB_addendum.md` — Group B (stress responses) and corrigendum to the buffered-response interpretation"""))

nb.cells = cells
nbf.write(nb, str(NB))
print(f"Wrote: {NB}")
print(f"Cells: {len(nb.cells)}")
