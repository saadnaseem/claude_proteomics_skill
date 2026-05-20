"""Categorize interaction-significant proteins into ALE-response patterns and plot top profiles."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import load_config, metadata_df, ensure_dirs, ROOT, palette, run_order


def categorize(profile: pd.Series) -> str:
    """profile indexed by KT-glu, KT-Ch, KT-Bu, HGL-glu, HGL-Ch, HGL-Bu (log2 means).

    Rule priority (after audit, 2026-05-05):
      1. reciprocal     — sign-flip in either medium with magnitude difference ≥ 1
      2. buffered/hyper — when KT and HGL respond in same direction but with different amplitudes
      3. constitutive   — when amplitudes are similar but baselines differ by ≥ 1
      4. other          — none of the above
    """
    kt_baseline  = profile["KT-glu"]
    hgl_baseline = profile["HGL-glu"]
    delta_kt_ch  = profile["KT-Ch"]  - kt_baseline
    delta_hgl_ch = profile["HGL-Ch"] - hgl_baseline
    delta_kt_bu  = profile["KT-Bu"]  - kt_baseline
    delta_hgl_bu = profile["HGL-Bu"] - hgl_baseline
    delta_baseline = hgl_baseline - kt_baseline

    # 1. Sign-flip / reciprocal
    if (delta_kt_ch * delta_hgl_ch < 0 and abs(delta_kt_ch - delta_hgl_ch) >= 1) or \
       (delta_kt_bu * delta_hgl_bu < 0 and abs(delta_kt_bu - delta_hgl_bu) >= 1):
        return "reciprocal"

    # 2. Same-sign amplitude differences
    avg_kt  = 0.5 * (abs(delta_kt_ch)  + abs(delta_kt_bu))
    avg_hgl = 0.5 * (abs(delta_hgl_ch) + abs(delta_hgl_bu))
    if avg_hgl > avg_kt + 0.5:
        return "hyper_induced_in_HGL"
    if avg_kt > avg_hgl + 0.5:
        return "buffered_in_HGL"

    # 3. Parallel-curve constitutive shift
    if abs(delta_baseline) >= 1.0:
        return "constitutive_HGL_shift"

    return "other"


def main() -> None:
    cfg = load_config()
    ensure_dirs(cfg)
    meta = metadata_df(cfg).set_index("run_id").loc[run_order(cfg)]
    log2 = pd.read_csv(ROOT / "outputs" / "tables" / "protein_log2_wide.tsv",
                       sep="\t", index_col=[0, 1, 2])
    inter = pd.read_csv(ROOT / "outputs" / "tables" / "interaction_significant_proteins.tsv",
                        sep="\t")
    print(f"[interaction] {len(inter)} interaction-significant proteins")
    if inter.empty:
        print("No proteins with q_interaction <= 0.05; nothing to do.")
        return

    # Compute group means for each interaction protein
    group_cols = ["KT-glu", "KT-Ch", "KT-Bu", "HGL-glu", "HGL-Ch", "HGL-Bu"]
    group_means = {}
    for sid, info in cfg["samples"].items():
        runs = meta.loc[meta["sample"] == sid].index.tolist()
        group_means[info["label"]] = log2[runs].mean(axis=1)
    group_means = pd.DataFrame(group_means)[group_cols]

    inter_idx = inter[["Protein.Group", "Protein", "Description"]].apply(tuple, axis=1).tolist()
    available = [k for k in inter_idx if k in group_means.index]
    profiles = group_means.loc[available]

    # Assign categories
    cats = profiles.apply(categorize, axis=1)
    inter_set = inter.set_index(["Protein.Group", "Protein", "Description"]).loc[available].copy()
    inter_set["category"] = cats.values
    inter_set[group_cols] = profiles[group_cols].values
    inter_set.to_csv(ROOT / "outputs" / "tables" / "interaction_categorized.tsv", sep="\t")
    counts = inter_set["category"].value_counts()
    print("\n[category counts]")
    print(counts.to_string())

    # Bar plot of category counts
    fig, ax = plt.subplots(figsize=(7, 4))
    cat_order = ["constitutive_HGL_shift", "buffered_in_HGL", "hyper_induced_in_HGL",
                 "reciprocal", "other"]
    counts = counts.reindex(cat_order, fill_value=0)
    colors = {"constitutive_HGL_shift": "#4477AA",
              "buffered_in_HGL": "#228833",
              "hyper_induced_in_HGL": "#EE6677",
              "reciprocal": "#CCBB44",
              "other": "#888888"}
    ax.bar(range(len(counts)), counts.values,
           color=[colors[c] for c in counts.index])
    for i, v in enumerate(counts.values):
        ax.text(i, v + 0.2, str(v), ha="center", va="bottom", fontsize=10)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("# interaction-significant proteins")
    ax.set_title(f"Interaction-significant proteins (n={len(inter_set)}) by ALE-response category")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F12b_interaction_categories.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    # Profile grid for top interaction proteins
    inter_top = inter_set.sort_values("q_interaction").head(24)
    n = len(inter_top)
    cols = 4; rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.4, rows * 2.5),
                             sharex=True)
    axes = np.array(axes).flatten()
    x = np.arange(3)
    media = ["glucose", "Chlys", "butylamine"]
    for ax, ((pg, pname, desc), row) in zip(axes, inter_top.iterrows()):
        kt_y  = [row["KT-glu"],  row["KT-Ch"],  row["KT-Bu"]]
        hgl_y = [row["HGL-glu"], row["HGL-Ch"], row["HGL-Bu"]]
        ax.plot(x, kt_y,  "o-", color=palette(cfg, "strain")["KT2440"],
                label="KT2440", lw=1.5, markersize=6)
        ax.plot(x, hgl_y, "s--", color=palette(cfg, "strain")["HGL1175"],
                label="HGL1175", lw=1.5, markersize=6)
        ax.set_xticks(x); ax.set_xticklabels(media, fontsize=7)
        title = f"{pname}\nq_int={row['q_interaction']:.2g} | {row['category'].replace('_',' ')}"
        ax.set_title(title, fontsize=7.5)
        ax.set_ylabel("log2 mean")
        ax.grid(alpha=0.2)
    for ax in axes[n:]: ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, -0.01), frameon=False, fontsize=10)
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F12_interaction_profiles_top24.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()
    print("[done]")


if __name__ == "__main__":
    main()
