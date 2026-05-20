"""Curated tolerance/metabolism module scoring + radar/heatmap visualization."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (load_config, load_modules, metadata_df, ensure_dirs,
                     ROOT, palette, run_order)


def gene_matches(gene_field: str, symbols: set, prefixes: list) -> bool:
    if not isinstance(gene_field, str): return False
    tokens = [t.strip().lower() for t in gene_field.replace(",", " ").split()]
    if any(t in symbols for t in tokens): return True
    for t in tokens:
        for p in prefixes:
            if t.startswith(p): return True
    return False


def main() -> None:
    cfg = load_config()
    mods = load_modules()
    ensure_dirs(cfg)
    meta = metadata_df(cfg).set_index("run_id").loc[run_order(cfg)]

    log2 = pd.read_csv(ROOT / "outputs" / "tables" / "protein_log2_wide.tsv",
                       sep="\t", index_col=[0, 1, 2])
    annot = pd.read_csv(ROOT / "data" / "uniprot_annotation.tsv", sep="\t")
    # Map first-accession -> gene_names string
    acc_to_gene = dict(zip(annot["Entry"], annot["Gene Names"].fillna("")))

    # Assign each detected protein its gene field (first accession in group)
    proteins_meta = pd.DataFrame({
        "Protein.Group": log2.index.get_level_values(0),
        "Protein":       log2.index.get_level_values(1),
        "Description":   log2.index.get_level_values(2),
    })
    proteins_meta["first_acc"] = proteins_meta["Protein.Group"].str.split(";").str[0]
    proteins_meta["genes"]     = proteins_meta["first_acc"].map(acc_to_gene).fillna("")

    # Module membership
    membership_rows = []
    module_proteins = {}
    for mod_name, info in mods.items():
        symbols  = set(s.lower() for s in info.get("genes", []))
        prefixes = [p.lower() for p in info.get("gene_prefixes", [])]
        sel = proteins_meta[proteins_meta["genes"].apply(
            lambda g: gene_matches(g, symbols, prefixes))]
        module_proteins[mod_name] = sel
        for _, row in sel.iterrows():
            membership_rows.append({"module": mod_name,
                                    "Protein.Group": row["Protein.Group"],
                                    "Protein": row["Protein"],
                                    "Description": row["Description"],
                                    "genes": row["genes"]})
        print(f"  module {mod_name:30s}: {len(sel):3d} proteins matched")

    pd.DataFrame(membership_rows).to_csv(
        ROOT / "outputs" / "tables" / "module_membership.tsv", sep="\t", index=False)

    # Compute per-run module score = mean log2 intensity over module proteins,
    # then z-score across the 18 runs.
    run_scores = {}
    for mod_name, sel in module_proteins.items():
        if len(sel) < 2:
            continue
        idx = list(zip(sel["Protein.Group"], sel["Protein"], sel["Description"]))
        sub = log2.loc[idx].dropna(how="all")
        run_scores[mod_name] = sub.mean(axis=0)
    score_df = pd.DataFrame(run_scores).T
    score_df = score_df[run_order(cfg)]
    z = score_df.sub(score_df.mean(axis=1), axis=0).div(score_df.std(axis=1), axis=0)
    z.to_csv(ROOT / "outputs" / "tables" / "module_zscores_per_run.tsv", sep="\t")

    # Group means
    group_z = pd.DataFrame({
        cfg["samples"][s]["label"]:
            z[meta.loc[meta["sample"] == s].index.tolist()].mean(axis=1)
        for s in cfg["samples"]
    })
    col_order = ["KT-glu", "KT-Ch", "KT-Bu", "HGL-glu", "HGL-Ch", "HGL-Bu"]
    group_z = group_z[col_order]
    group_z.to_csv(ROOT / "outputs" / "tables" / "module_zscores_per_group.tsv", sep="\t")

    # Heatmap of module × group
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(group_z, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                annot=True, fmt=".1f", cbar_kws={"label": "z-score (log2)"},
                ax=ax, linewidths=0.4, linecolor="white")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Tolerance / metabolism module scores (z-score across 6 groups)")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F10_module_heatmap.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    # Two-way ANOVA on each module score: sample-level z-score ~ strain * medium
    rows = []
    fit_meta = meta.reset_index()
    for mod, scores in z.iterrows():
        df = fit_meta.copy()
        df["y"] = scores.values
        df = df.dropna(subset=["y"])
        if df["y"].nunique() < 4: continue
        try:
            anv = sm.stats.anova_lm(smf.ols("y ~ C(strain) * C(medium)", data=df).fit(), typ=2)
            rows.append({
                "module": mod,
                "p_strain":      anv.loc["C(strain)",           "PR(>F)"],
                "p_medium":      anv.loc["C(medium)",           "PR(>F)"],
                "p_interaction": anv.loc["C(strain):C(medium)", "PR(>F)"],
            })
        except Exception: continue
    mod_anova = pd.DataFrame(rows)
    for term in ["strain", "medium", "interaction"]:
        valid = mod_anova[f"p_{term}"].notna()
        q = np.full(len(mod_anova), np.nan)
        if valid.sum():
            _, q_v, _, _ = multipletests(mod_anova.loc[valid, f"p_{term}"].values,
                                         method="fdr_bh")
            q[valid.values] = q_v
        mod_anova[f"q_{term}"] = q
    mod_anova.to_csv(ROOT / "outputs" / "tables" / "module_anova.tsv",
                     sep="\t", index=False)
    print("\n[module ANOVA] modules with q≤0.10:")
    print(mod_anova[(mod_anova["q_strain"]<=0.1)|(mod_anova["q_medium"]<=0.1)|
                   (mod_anova["q_interaction"]<=0.1)].to_string(index=False))

    # Radar plot (one polygon per group)
    radar(group_z, cfg)
    print("\n[done]")


def radar(group_z: pd.DataFrame, cfg: dict) -> None:
    modules = group_z.index.tolist()
    angles = np.linspace(0, 2 * np.pi, len(modules), endpoint=False).tolist()
    angles += angles[:1]

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5),
                             subplot_kw=dict(projection="polar"))
    palette_strain = cfg["plot"]["palette_strain"]
    medium_styles = {"glucose": "-", "Chlys": "--", "butylamine": ":"}

    for ax, strain in zip(axes, ["KT2440", "HGL1175"]):
        prefix = "KT" if strain == "KT2440" else "HGL"
        for med, style in medium_styles.items():
            col = {"glucose": f"{prefix}-glu",
                   "Chlys":   f"{prefix}-Ch",
                   "butylamine": f"{prefix}-Bu"}[med]
            vals = group_z[col].tolist() + [group_z[col].iloc[0]]
            ax.plot(angles, vals, style, label=med, linewidth=2,
                    color=palette_strain[strain])
            ax.fill(angles, vals, alpha=0.08, color=palette_strain[strain])
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels([m.replace("_", "\n") for m in modules], fontsize=7)
        ax.set_ylim(-2, 2)
        ax.set_title(strain, color=palette_strain[strain], fontsize=12, pad=18)
        ax.legend(loc="upper right", bbox_to_anchor=(1.25, 1.1), fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("Module scores per group (z across all 6 groups)", fontsize=12)
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F11_module_radar.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
