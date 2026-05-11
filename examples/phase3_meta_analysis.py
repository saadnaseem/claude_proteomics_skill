"""Phase 3 — cross-comparison meta-analysis across all 9 unique strain×condition comparisons.
Builds master DE matrix, constitutive/conditional decomposition, buffered-response scatters,
regulator co-expression heatmap, tolerance candidate ranking."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json, os
from pathlib import Path
from collections import defaultdict

META_DIR = Path(open("/tmp/meta_dir.txt").read().strip())
FIG = META_DIR / "figures"; TBL = META_DIR / "tables"

DATA_DIR = "/Users/snaseem/Library/CloudStorage/GoogleDrive-snaseem@lbl.gov/My Drive/Saad m-group/proteomics/ALE_claude/ALE CHlys BAPRT1296_JBEI_20250820_SNaseem__with_full_ids.pr_matrix/t-test_files/"

# Sample IDs → human labels
SAMPLE = {
    "SN_0725_80": "KT/glu",  "SN_0725_83": "KT/Chlys", "SN_0725_86": "KT/buta",
    "SN_0725_89": "HGL/glu", "SN_0725_92": "HGL/Chlys","SN_0725_95": "HGL/buta",
}

# Map filename → (numerator_id, denominator_id, group)
# Group A = strain effect (HGL_vs_KT under each medium)
# Group B = stress effect (medium_vs_glucose per strain)
# Group C = hydrolysate vs hydrolysate per strain  (we have 2 unique per strain — one is sign-flip)
COMPARISONS = {
    # Group A
    "A_glu_HGLvsKT":        ("t-test_glucose_HGL1175_vs_KT_20250902-195343.xlsx", "SN_0725_89","SN_0725_80","A"),
    "A_Chlys_HGLvsKT":      ("t-test_Chlys_HGL1175_vs_KT_20250902-195343.xlsx",   "SN_0725_92","SN_0725_83","A"),
    "A_buta_HGLvsKT":       ("t-test_Butamine_HGL1175_vs_KT_20250902-195343.xlsx","SN_0725_95","SN_0725_86","A"),
    # Group B
    "B_KT_Chlys_vs_glu":    ("t-test_KT_Chlys_vs_glucose_20250902-195343.xlsx",    "SN_0725_83","SN_0725_80","B"),
    "B_KT_buta_vs_glu":     ("t-test_KT_butamine_vs_glucose_20250902-195343.xlsx", "SN_0725_86","SN_0725_80","B"),
    "B_HGL_Chlys_vs_glu":   ("t-test_HGL1175_Chlys_vs_glucose_20250902-195343.xlsx","SN_0725_92","SN_0725_89","B"),
    "B_HGL_buta_vs_glu":    ("t-test_HGL1175_butamine_vs_glucose_20250902-195343.xlsx","SN_0725_95","SN_0725_89","B"),
    # Group C — pick one orientation per pair (Chlys_vs_butamine; the other is sign-flip)
    "C_KT_Chlys_vs_buta":   ("t-test_KT_Chlys_vs_butamine_20250902-195343.xlsx",   "SN_0725_83","SN_0725_86","C"),
    "C_HGL_Chlys_vs_buta":  ("t-test_HGL1175_Chlys_vs_butamine_20250902-195343.xlsx","SN_0725_92","SN_0725_95","C"),
}

# === 1. Build master DE matrix ===
print("[1] Building master DE matrix")
master = None
meta_dict = {}  # gene metadata

for key, (fname, a_id, b_id, group) in COMPARISONS.items():
    fpath = os.path.join(DATA_DIR, fname)
    df = pd.read_excel(fpath, sheet_name="Full t-test output")
    df["uniprot_acc"] = df["Protein.Group"].apply(lambda s: str(s).split(";")[0].strip())

    # Validate sign convention
    implied = df[f"log2_mean_{a_id}"] - df[f"log2_mean_{b_id}"]
    err = (df["log2_Fold_change_A/B"] - implied).abs().max()
    if err > 0.01:
        print(f"  WARN sign convention {key}: err={err}")

    # Build per-comparison columns
    cmp_df = df[["Protein","uniprot_acc","Protein.Names","Protein.Description"]].copy()
    cmp_df.columns = ["protein","acc","gene","desc"]
    cmp_df[f"{key}_log2fc"] = df["log2_Fold_change_A/B"]
    cmp_df[f"{key}_padj"]   = df["p_adjusted(BH)"]
    cmp_df[f"{key}_p"]      = df["p-value"]

    # Update gene metadata (first non-null wins)
    for _, r in cmp_df.iterrows():
        a = r["acc"]
        if a not in meta_dict:
            meta_dict[a] = {"protein": r["protein"], "gene": r["gene"], "desc": r["desc"]}

    # Merge into master
    cmp_subset = cmp_df[["acc", f"{key}_log2fc", f"{key}_padj", f"{key}_p"]]
    if master is None:
        master = cmp_subset
    else:
        master = master.merge(cmp_subset, on="acc", how="outer")

    n_sig_strict = ((df["p_adjusted(BH)"]<0.05) & (df["log2_Fold_change_A/B"].abs()>1)).sum()
    print(f"  {key:<24} A={SAMPLE[a_id]:<10} B={SAMPLE[b_id]:<10} group={group}  sig_strict={n_sig_strict}")

# Add metadata columns
master["protein"] = master["acc"].map(lambda a: meta_dict.get(a,{}).get("protein",""))
master["gene"]    = master["acc"].map(lambda a: meta_dict.get(a,{}).get("gene",""))
master["desc"]    = master["acc"].map(lambda a: meta_dict.get(a,{}).get("desc",""))

cols_order = ["acc","protein","gene","desc"] + [c for c in master.columns if c not in ("acc","protein","gene","desc")]
master = master[cols_order]
master.to_csv(TBL / "01_master_de_matrix.csv", index=False)
print(f"  Master matrix: {master.shape} (proteins x columns)")

# === 2. Constitutive vs conditional ALE signature (Group A decomposition) ===
print("\n[2] Constitutive vs conditional ALE signature (Group A)")
# Significance: p_adj<0.05 & |log2FC|>1 (strict) or p_adj<0.05 (expanded for sparse Chlys)
def sig_strict(df, prefix):
    return (df[f"{prefix}_padj"]<0.05) & (df[f"{prefix}_log2fc"].abs() > 1)
def sig_lenient(df, prefix, padj=0.10):
    return df[f"{prefix}_padj"] < padj

a_keys = ["A_glu_HGLvsKT","A_Chlys_HGLvsKT","A_buta_HGLvsKT"]
sigA = {k: sig_strict(master, k) for k in a_keys}
sigA_l = {k: sig_lenient(master, k, 0.10) for k in a_keys}

# Decompose
all_sig = sigA["A_glu_HGLvsKT"] | sigA["A_Chlys_HGLvsKT"] | sigA["A_buta_HGLvsKT"]
in_glu  = sigA["A_glu_HGLvsKT"]
in_buta = sigA["A_buta_HGLvsKT"]
in_chl  = sigA["A_Chlys_HGLvsKT"]

print(f"  Total ever-sig (any Group A, strict): {all_sig.sum()}")

categories = {
    "all3 (strict, constitutive)": (in_glu & in_buta & in_chl).sum(),
    "glu+buta only": (in_glu & in_buta & ~in_chl).sum(),
    "glu+chlys only": (in_glu & in_chl & ~in_buta).sum(),
    "buta+chlys only": (in_buta & in_chl & ~in_glu).sum(),
    "glu only": (in_glu & ~in_buta & ~in_chl).sum(),
    "buta only": (in_buta & ~in_glu & ~in_chl).sum(),
    "chlys only": (in_chl & ~in_glu & ~in_buta).sum(),
}
for k, v in categories.items():
    print(f"    {k:<35} {v}")

# Save the multi-significant proteins (in ≥2 Group A comparisons)
multisig = master[in_glu.astype(int)+in_buta.astype(int)+in_chl.astype(int) >= 2].copy()
multisig["n_groupA_sig"] = (in_glu.astype(int)+in_buta.astype(int)+in_chl.astype(int))[multisig.index]
multisig_cols = ["acc","gene","protein","desc","n_groupA_sig",
                 "A_glu_HGLvsKT_log2fc","A_glu_HGLvsKT_padj",
                 "A_buta_HGLvsKT_log2fc","A_buta_HGLvsKT_padj",
                 "A_Chlys_HGLvsKT_log2fc","A_Chlys_HGLvsKT_padj"]
multisig[multisig_cols].sort_values("n_groupA_sig", ascending=False).to_csv(TBL/"02_constitutive_ALE_hits.csv", index=False)
print(f"  Multi-sig in Group A (≥2 of 3): {len(multisig)} proteins → tables/02_constitutive_ALE_hits.csv")
print(f"\n  All-3 constitutive (strict in glu+buta+chlys):")
all3 = master[in_glu & in_buta & in_chl]
for _, r in all3.iterrows():
    print(f"    {r['acc']:<10} {str(r['gene'])[:14]:<14}  glu={r['A_glu_HGLvsKT_log2fc']:+.2f}  buta={r['A_buta_HGLvsKT_log2fc']:+.2f}  Chlys={r['A_Chlys_HGLvsKT_log2fc']:+.2f}  | {str(r['desc'])[:50]}")

# Also report glu+buta intersection (likely the strongest constitutive set given Chlys sparsity)
print(f"\n  glu+buta intersection (strict in both, Chlys may be subthreshold):")
glubuta = master[in_glu & in_buta].sort_values("A_glu_HGLvsKT_log2fc", ascending=False)
for _, r in glubuta.iterrows():
    print(f"    {r['acc']:<10} {str(r['gene'])[:14]:<14}  glu={r['A_glu_HGLvsKT_log2fc']:+.2f}  buta={r['A_buta_HGLvsKT_log2fc']:+.2f}  | {str(r['desc'])[:55]}")

# === 3. Buffered-response 2D scatters ===
print("\n[3] Buffered-response analysis (WT vs ALE stress responses)")
for hyd, kt_key, hgl_key in [("Chlys", "B_KT_Chlys_vs_glu", "B_HGL_Chlys_vs_glu"),
                              ("butamine","B_KT_buta_vs_glu","B_HGL_buta_vs_glu")]:
    sub = master[[ "acc","gene","desc",
                   f"{kt_key}_log2fc", f"{kt_key}_padj",
                   f"{hgl_key}_log2fc", f"{hgl_key}_padj"]].dropna()
    sub.columns = ["acc","gene","desc","kt_log2fc","kt_padj","hgl_log2fc","hgl_padj"]

    # Classify
    SIG_P = 0.10
    sub["kt_sig"]  = sub["kt_padj"] < SIG_P
    sub["hgl_sig"] = sub["hgl_padj"] < SIG_P
    sub["category"] = "ns"
    sub.loc[sub["kt_sig"] & ~sub["hgl_sig"], "category"]   = "buffered_by_ALE"  # WT responds, ALE doesn't
    sub.loc[~sub["kt_sig"] & sub["hgl_sig"], "category"]   = "ALE_specific"     # ALE responds, WT doesn't
    sub.loc[sub["kt_sig"] & sub["hgl_sig"], "category"]    = "shared_response"

    sub["delta"] = sub["kt_log2fc"] - sub["hgl_log2fc"]  # how much MORE WT changes than ALE
    sub.to_csv(TBL/f"03_buffered_{hyd}.csv", index=False)

    # Counts
    counts = sub["category"].value_counts()
    print(f"  {hyd}:")
    for c, n in counts.items():
        print(f"    {c}: {n}")

    # Top buffered hits (most strongly KT-responsive but ALE-flat)
    buf = sub[sub["category"]=="buffered_by_ALE"].copy()
    buf["abs_kt"] = buf["kt_log2fc"].abs()
    top_buf = buf.nlargest(15, "abs_kt")
    print(f"  Top 10 buffered by ALE under {hyd} (WT changes a lot, ALE barely):")
    for _, r in top_buf.head(10).iterrows():
        direction = "↑" if r["kt_log2fc"]>0 else "↓"
        print(f"    {r['acc']:<10} {str(r['gene'])[:14]:<14}  KT={r['kt_log2fc']:+.2f} {direction}, HGL={r['hgl_log2fc']:+.2f}  | {str(r['desc'])[:55]}")

    # Plot
    fig, ax = plt.subplots(figsize=(9, 8))
    color_map = {"ns":"#CCCCCC", "buffered_by_ALE":"#D62728",
                 "ALE_specific":"#9467BD", "shared_response":"#2CA02C"}
    for cat, color in color_map.items():
        sm = sub[sub["category"]==cat]
        size = 12 if cat=="ns" else 35
        alpha = 0.25 if cat=="ns" else 0.8
        ax.scatter(sm["kt_log2fc"], sm["hgl_log2fc"], s=size, c=color, alpha=alpha,
                   edgecolors="none", label=f"{cat} (n={len(sm)})")

    # Diagonal & axes
    lim = max(sub["kt_log2fc"].abs().max(), sub["hgl_log2fc"].abs().max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], color="grey", ls="--", lw=0.8, alpha=0.5, label="x=y (same response)")
    ax.axhline(0, color="black", lw=0.5); ax.axvline(0, color="black", lw=0.5)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

    # Label the top buffered hits
    for _, r in top_buf.head(15).iterrows():
        ax.annotate(str(r["gene"]) if pd.notna(r["gene"]) and str(r["gene"]) else r["acc"],
                    (r["kt_log2fc"], r["hgl_log2fc"]),
                    fontsize=7, xytext=(3,3), textcoords="offset points")

    ax.set_xlabel(f"WT KT2440: log2FC ({hyd} vs glucose)")
    ax.set_ylabel(f"ALE HGL1175: log2FC ({hyd} vs glucose)")
    ax.set_title(f"Buffered-response scatter: {hyd} stress\n"
                 f"Off-diagonal in red box = WT responds, ALE doesn't = TOLERANCE CANDIDATES")
    ax.legend(loc="upper left", fontsize=9)
    plt.tight_layout()
    plt.savefig(FIG/f"03_buffered_response_{hyd}.png")
    plt.close()
    print(f"  Saved: figures/03_buffered_response_{hyd}.png")

# === 4. Regulator co-expression heatmap ===
print("\n[4] Regulator co-expression across all 9 comparisons")
regulators = {
    "TurA (PP_1366, MvaT P16)":   "Q88N50",
    "TurB (PP_3765, MvaT H-NS)":  "Q88GF9",
    "RelA (ppGpp synthase)":      "Q88MB8",
    "CsrA":                       None,  # find by gene name
    "KatG (catalase-peroxidase)": None,
    "GroES":                      None,
    "DnaK":                       None,
    "Lon protease":               None,
    "RpoS (sigma 38)":            None,
    "Fur (iron uptake reg.)":     None,
    "PhaG (PHA granule)":         "Q88D20",
    "PhaF GA2":                   "Q88D21",
    "Pp_4834 (SPFH)":             "Q88DJ1",
    "OmpA family Pp_1502":        "Q88MR7",
    "Cyoc (cyt bo3)":             "Q88PN5",
    "Idh (NADP-IDH)":             "Q88FS1",
    "Pp_0241 (myst transporter)": "Q88R92",
}
def find_acc(gene_substr):
    hits = master[master["gene"].astype(str).str.contains(gene_substr, case=False, na=False) |
                  master["desc"].astype(str).str.contains(gene_substr, case=False, na=False)]
    return hits["acc"].tolist()

for name, acc in regulators.items():
    if acc is None:
        gn = name.split()[0]
        accs = find_acc(gn)
        if accs:
            regulators[name] = accs[0]

reg_rows = []
for name, acc in regulators.items():
    if acc and acc in master["acc"].values:
        row = master[master["acc"]==acc].iloc[0]
        rdict = {"name": name, "acc": acc, "gene": row["gene"], "desc": str(row["desc"])[:50]}
        for k in COMPARISONS:
            rdict[k] = row[f"{k}_log2fc"] if pd.notna(row[f"{k}_log2fc"]) else np.nan
            rdict[f"{k}_p"] = row[f"{k}_padj"] if pd.notna(row[f"{k}_padj"]) else np.nan
        reg_rows.append(rdict)

reg_df = pd.DataFrame(reg_rows)
reg_df.to_csv(TBL/"04_regulators_across_comparisons.csv", index=False)

# Heatmap
fig, ax = plt.subplots(figsize=(13, 9))
log2fc_cols = list(COMPARISONS.keys())
heat_data = reg_df[log2fc_cols].copy()
heat_data.index = reg_df["name"]
heat_data.columns = [c.replace("_log2fc","").replace("HGLvsKT","HGL/KT").replace("_vs_","→")
                     for c in log2fc_cols]
# Annotate with significance
annot = pd.DataFrame(index=heat_data.index, columns=heat_data.columns, dtype=str)
for i, name in enumerate(heat_data.index):
    for j, k in enumerate(log2fc_cols):
        val = reg_df.iloc[i][k]
        p = reg_df.iloc[i][f"{k}_p"]
        s = ""
        if pd.notna(p):
            if p<0.001: s = "***"
            elif p<0.01: s = "**"
            elif p<0.05: s = "*"
            elif p<0.10: s = "·"
        annot.iloc[i,j] = f"{val:+.1f}{s}" if pd.notna(val) else ""
sns.heatmap(heat_data, annot=annot, fmt="", cmap="RdBu_r", center=0,
            vmin=-3, vmax=3, ax=ax, cbar_kws={"label": "log2FC"},
            linewidths=0.3, linecolor="white")
ax.set_title("Regulator + key-gene log2FC across all 9 comparisons\n"
             "*** p_adj<0.001  ** p_adj<0.01  * p_adj<0.05  · p_adj<0.10",
             fontsize=11)
plt.tight_layout()
plt.savefig(FIG/"04_regulator_heatmap.png")
plt.close()
print(f"  Saved: figures/04_regulator_heatmap.png ({len(reg_df)} genes x {len(log2fc_cols)} comparisons)")

# === 5. Tolerance candidate ranking ===
print("\n[5] Tolerance candidate ranking")
# Score: (|HGL-vs-KT under hydrolysate|) × (|KT_hydrolysate_vs_glucose|)
# = how much ALE has changed AND how much WT actively responds to that protein under stress

scoring = master[["acc","gene","desc"]].copy()
for hyd, a_key, b_key in [("Chlys", "A_Chlys_HGLvsKT", "B_KT_Chlys_vs_glu"),
                          ("buta",  "A_buta_HGLvsKT",  "B_KT_buta_vs_glu")]:
    a_lfc = master[f"{a_key}_log2fc"].abs().fillna(0)
    b_lfc = master[f"{b_key}_log2fc"].abs().fillna(0)
    scoring[f"{hyd}_strain_eff"] = master[f"{a_key}_log2fc"]
    scoring[f"{hyd}_WT_resp"] = master[f"{b_key}_log2fc"]
    scoring[f"{hyd}_score"] = a_lfc * b_lfc

# Composite score = mean of the two hydrolysates
scoring["composite_score"] = (scoring["Chlys_score"] + scoring["buta_score"]) / 2

# Also bring in the glucose strain effect (constitutive baseline)
scoring["glu_strain_eff"] = master["A_glu_HGLvsKT_log2fc"]

top_candidates = scoring.nlargest(30, "composite_score")
top_candidates.to_csv(TBL/"05_tolerance_candidates.csv", index=False)
print("  Top 15 tolerance-candidate proteins (high ALE-effect × WT-response):")
print(f"  {'Acc':<10} {'Gene':<14} {'Glu strain':<12} {'Chlys strain':<14} {'Chlys WT-resp':<14} {'Buta strain':<14} {'Buta WT-resp':<14} {'Score':<8} {'Description'}")
for _, r in top_candidates.head(15).iterrows():
    print(f"  {r['acc']:<10} {str(r['gene'])[:14]:<14} {r['glu_strain_eff']:>+11.2f} {r['Chlys_strain_eff']:>+13.2f} {r['Chlys_WT_resp']:>+13.2f} {r['buta_strain_eff']:>+13.2f} {r['buta_WT_resp']:>+13.2f} {r['composite_score']:>7.2f}  {str(r['desc'])[:45]}")

print("\n--- Phase 3 done ---")
