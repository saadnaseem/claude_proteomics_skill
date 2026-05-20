"""Exploratory: PCA (annotated), hierarchical dendrogram, per-protein variance partition."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (load_config, metadata_df, ensure_dirs, ROOT, palette, run_order)


def main() -> None:
    cfg = load_config()
    ensure_dirs(cfg)
    meta = metadata_df(cfg).set_index("run_id").loc[run_order(cfg)]
    log2 = pd.read_csv(ROOT / "outputs" / "tables" / "protein_log2_wide.tsv",
                       sep="\t", index_col=[0, 1, 2])

    # Use only complete-case proteins for PCA / clustering (avoids imputation choice)
    complete = log2.dropna()
    print(f"[PCA input] {complete.shape[0]} complete-case proteins (of {log2.shape[0]})")

    # 1. PCA
    X = complete.T.values  # runs × proteins
    X_centered = X - X.mean(axis=0)
    pca = PCA(n_components=4)
    scores = pca.fit_transform(X_centered)
    pcs = pd.DataFrame(scores, columns=[f"PC{i+1}" for i in range(4)],
                       index=complete.columns)
    pcs = pcs.join(meta[["sample", "label", "strain", "medium"]])
    pcs.to_csv(ROOT / "outputs" / "tables" / "pca_scores.tsv", sep="\t")
    var = pca.explained_variance_ratio_ * 100
    print(f"[PCA] explained variance %: {var.round(1).tolist()}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, (xpc, ypc) in zip(axes, [("PC1", "PC2"), ("PC1", "PC3")]):
        for med in pcs["medium"].unique():
            for st in pcs["strain"].unique():
                sel = pcs[(pcs["medium"] == med) & (pcs["strain"] == st)]
                ax.scatter(sel[xpc], sel[ypc],
                           color=palette(cfg, "strain")[st],
                           edgecolor="black",
                           marker={"glucose": "o", "Chlys": "s", "butylamine": "^"}[med],
                           s=140, label=f"{st} / {med}", alpha=0.85)
        ix = "PC1 PC2 PC3 PC4".split().index(xpc)
        iy = "PC1 PC2 PC3 PC4".split().index(ypc)
        ax.set_xlabel(f"{xpc} ({var[ix]:.1f} %)")
        ax.set_ylabel(f"{ypc} ({var[iy]:.1f} %)")
        ax.axhline(0, ls=":", lw=0.5, color="gray")
        ax.axvline(0, ls=":", lw=0.5, color="gray")
        ax.grid(alpha=0.2)
    handles, labels = axes[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc="lower center",
               ncol=6, bbox_to_anchor=(0.5, -0.04), frameon=False, fontsize=9)
    fig.suptitle("PCA of 18 runs (complete-case log2 intensity)", fontsize=12)
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F1_PCA.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    # Scree
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(range(1, 5), var, color="#4477AA")
    ax.plot(range(1, 5), var.cumsum(), "ko-", label="cumulative")
    ax.set_xlabel("PC"); ax.set_ylabel("% variance"); ax.legend()
    ax.set_title("Scree plot")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F1b_PCA_scree.png",
                dpi=cfg["plot"]["dpi"])
    plt.close()

    # 2. Hierarchical clustering dendrogram
    dist = pdist(X_centered, metric="correlation")
    Z = hierarchy.linkage(dist, method="average")
    fig, ax = plt.subplots(figsize=(11, 4.5))
    labels = [f"{r}\n{meta.loc[r,'label']}" for r in complete.columns]
    label_colors = {r: palette(cfg, "strain")[meta.loc[r, "strain"]] for r in complete.columns}
    hierarchy.dendrogram(Z, labels=labels, leaf_rotation=90, ax=ax,
                         color_threshold=0)
    for lbl in ax.get_xmajorticklabels():
        run = lbl.get_text().split("\n")[0]
        lbl.set_color(palette(cfg, "strain")[meta.loc[run, "strain"]])
    ax.set_title("Hierarchical clustering of runs (1 - Pearson correlation)")
    ax.set_ylabel("distance")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F2c_dendrogram.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    # 3. Per-protein variance partition
    # For each protein with no missing values across the 18 runs, fit log2 ~ C(strain) + C(medium) + C(strain):C(medium)
    # then partition variance via type-II sums of squares.
    df_var = []
    meta_for_fit = meta.reset_index()
    for prot_key, row in complete.iterrows():
        df_fit = meta_for_fit.copy()
        df_fit["y"] = row.values
        try:
            full = smf.ols("y ~ C(strain) * C(medium)", data=df_fit).fit()
            ss_total = ((df_fit["y"] - df_fit["y"].mean()) ** 2).sum()
            ss_resid = (full.resid ** 2).sum()
            r_strain  = smf.ols("y ~ C(medium)", data=df_fit).fit().resid
            r_medium  = smf.ols("y ~ C(strain)", data=df_fit).fit().resid
            r_inter   = smf.ols("y ~ C(strain) + C(medium)", data=df_fit).fit().resid
            ss_strain = (r_strain ** 2).sum() - ss_resid
            ss_medium = (r_medium ** 2).sum() - ss_resid
            ss_inter  = (r_inter  ** 2).sum() - ss_resid
            ss_strain = max(ss_strain - ss_inter, 0)
            ss_medium = max(ss_medium - ss_inter, 0)
            df_var.append({
                "protein": prot_key[1],
                "strain_pct":      100 * ss_strain / ss_total if ss_total > 0 else np.nan,
                "medium_pct":      100 * ss_medium / ss_total if ss_total > 0 else np.nan,
                "interaction_pct": 100 * ss_inter  / ss_total if ss_total > 0 else np.nan,
                "residual_pct":    100 * ss_resid  / ss_total if ss_total > 0 else np.nan,
            })
        except Exception:
            continue
    var_df = pd.DataFrame(df_var)
    var_df.to_csv(ROOT / "outputs" / "tables" / "variance_partition.tsv",
                  sep="\t", index=False)
    print("\n[variance partition] median % per term:")
    print(var_df[["strain_pct", "medium_pct", "interaction_pct", "residual_pct"]].median())

    fig, ax = plt.subplots(figsize=(7, 4.5))
    melt = var_df.melt(id_vars="protein", value_vars=["strain_pct", "medium_pct",
                                                       "interaction_pct", "residual_pct"],
                       var_name="term", value_name="pct")
    sns.violinplot(data=melt, x="term", y="pct", inner="quartile", ax=ax,
                   palette=["#4477AA", "#228833", "#EE6677", "#BBBBBB"], cut=0)
    ax.set_ylabel("% variance explained per protein")
    ax.set_xlabel("")
    ax.set_xticklabels(["strain", "medium", "strain×medium", "residual"])
    ax.set_title("Per-protein variance partitioning")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F3_variance_partition.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    print("\n[done] exploratory outputs in outputs/figures/ and outputs/tables/")


if __name__ == "__main__":
    main()
