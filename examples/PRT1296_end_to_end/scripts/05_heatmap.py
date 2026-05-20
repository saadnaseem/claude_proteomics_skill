"""Z-score heatmap of all DEPs (union across contrasts), hierarchically clustered."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import load_config, metadata_df, ensure_dirs, ROOT, palette, run_order


def main() -> None:
    cfg = load_config()
    ensure_dirs(cfg)
    meta = metadata_df(cfg).set_index("run_id").loc[run_order(cfg)]
    log2 = pd.read_csv(ROOT / "outputs" / "tables" / "protein_log2_wide.tsv",
                       sep="\t", index_col=[0, 1, 2])

    # Union of DEPs from any biology-relevant contrast at strict threshold
    master = pd.read_csv(ROOT / "outputs" / "tables" / "all_contrasts_long.tsv", sep="\t")
    deps = master.loc[master["sig_dir"].isin(["UP", "DOWN"]), "Protein.Group"].unique()
    print(f"[heatmap] {len(deps)} unique DEPs across all contrasts (strict q≤0.05, |log2FC|≥1)")

    log2_deps = log2[log2.index.get_level_values(0).isin(deps)].dropna()
    print(f"[heatmap] {log2_deps.shape[0]} DEPs with complete log2 values across 18 runs")

    # Collapse replicates → group-mean log2
    grp_means = pd.DataFrame({
        meta.loc[meta["sample"] == sid, "label"].iloc[0]:
            log2_deps[meta.loc[meta["sample"] == sid].index].mean(axis=1)
        for sid in cfg["samples"]
    })
    z = grp_means.sub(grp_means.mean(axis=1), axis=0).div(grp_means.std(axis=1), axis=0)
    z = z.dropna()

    # Cluster rows
    row_dist = pdist(z.values, metric="correlation")
    row_link = hierarchy.linkage(row_dist, method="average")
    n_clust = 6
    row_clusters = hierarchy.fcluster(row_link, t=n_clust, criterion="maxclust")

    # Save cluster membership
    cluster_df = pd.DataFrame({
        "Protein.Group":  z.index.get_level_values(0),
        "Protein":        z.index.get_level_values(1),
        "Description":    z.index.get_level_values(2),
        "cluster":        row_clusters,
    })
    cluster_df.to_csv(ROOT / "outputs" / "tables" / "heatmap_cluster_membership.tsv",
                      sep="\t", index=False)
    print("[clusters]")
    print(cluster_df["cluster"].value_counts().sort_index().to_string())

    # Plot with clustermap
    col_order = ["KT-glu", "KT-Ch", "KT-Bu", "HGL-glu", "HGL-Ch", "HGL-Bu"]
    z = z[col_order]
    col_strain = [c.split("-")[0] for c in col_order]
    col_medium = ["glucose" if "glu" in c else ("Chlys" if "Ch" in c else "butylamine")
                  for c in col_order]
    col_strain_full = ["KT2440" if s == "KT" else "HGL1175" for s in col_strain]
    col_strain_colors = [palette(cfg, "strain")[s] for s in col_strain_full]
    col_medium_colors = [palette(cfg, "medium")[m] for m in col_medium]
    col_colors = pd.DataFrame({"strain": col_strain_colors,
                                "medium": col_medium_colors}, index=col_order)
    cluster_palette = sns.color_palette("tab10", n_clust)
    row_colors = pd.Series([cluster_palette[c - 1] for c in row_clusters],
                           index=z.index, name="cluster")

    g = sns.clustermap(
        z,
        row_linkage=row_link,
        col_cluster=False,
        cmap="RdBu_r",
        center=0, vmin=-2, vmax=2,
        col_colors=col_colors,
        row_colors=row_colors,
        yticklabels=False,
        figsize=(7.5, 10),
        dendrogram_ratio=(0.12, 0.05),
        cbar_kws={"label": "row z-score (log2)"},
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel(f"{z.shape[0]} DEPs (clustered)")
    plt.suptitle("DEP heatmap (union across all contrasts)", y=1.01, fontsize=12)
    plt.savefig(ROOT / "outputs" / "figures" / "F8_heatmap_DEPs.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    # Also write z-scored matrix for downstream use
    z.to_csv(ROOT / "outputs" / "tables" / "heatmap_zscore_matrix.tsv", sep="\t")
    print("[done]")


if __name__ == "__main__":
    main()
