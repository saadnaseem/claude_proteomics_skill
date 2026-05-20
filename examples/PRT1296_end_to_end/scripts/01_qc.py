"""QC: sample correlation heatmap, missingness, per-group CV, intensity boxplots."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (load_config, metadata_df, ensure_dirs, ROOT, palette, run_order)


def main() -> None:
    cfg = load_config()
    ensure_dirs(cfg)
    meta = metadata_df(cfg).set_index("run_id").loc[run_order(cfg)]
    log2 = pd.read_csv(ROOT / "outputs" / "tables" / "protein_log2_wide.tsv",
                       sep="\t", index_col=[0, 1, 2])

    # 1. Sample-correlation heatmap
    corr = log2.corr(method="pearson", min_periods=200)
    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(corr, vmin=0.85, vmax=1.0, cmap="viridis",
                xticklabels=True, yticklabels=True, ax=ax,
                cbar_kws={"label": "Pearson r"})
    strain_colors = [palette(cfg, "strain")[meta.loc[r, "strain"]] for r in corr.columns]
    medium_colors = [palette(cfg, "medium")[meta.loc[r, "medium"]] for r in corr.columns]
    for i, c in enumerate(strain_colors):
        ax.add_patch(plt.Rectangle((i, len(corr)), 1, 0.6, color=c, transform=ax.transData,
                                    clip_on=False))
    for i, c in enumerate(medium_colors):
        ax.add_patch(plt.Rectangle((i, len(corr) + 0.7), 1, 0.6, color=c, transform=ax.transData,
                                    clip_on=False))
    ax.set_title("Sample correlation (Pearson on log2 intensity)")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F2_QC_sample_correlation.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    corr.to_csv(ROOT / "outputs" / "tables" / "qc_sample_correlation.tsv", sep="\t")
    print(f"[corr] median within-group r={corr_within(corr, meta):.3f}")

    # 2. Per-group CV
    cv_rows = []
    for sample, sub in meta.groupby("sample", sort=False):
        runs = sub.index.tolist()
        intensity = 2 ** log2[runs]
        cv = (intensity.std(axis=1) / intensity.mean(axis=1)) * 100
        cv = cv.dropna()
        cv_rows.append({"sample": sample, "label": sub["label"].iloc[0],
                        "n_proteins": len(cv),
                        "median_CV%": float(cv.median()),
                        "mean_CV%": float(cv.mean()),
                        "p90_CV%": float(cv.quantile(0.9))})
    cv_df = pd.DataFrame(cv_rows)
    cv_df.to_csv(ROOT / "outputs" / "tables" / "qc_per_group_CV.tsv", sep="\t", index=False)
    print("\n[per-group CV%]")
    print(cv_df.to_string(index=False))

    # 3. Per-run intensity boxplot
    fig, ax = plt.subplots(figsize=(10, 5))
    log2_long = log2.melt(var_name="run_id", value_name="log2_intensity").dropna()
    log2_long = log2_long.merge(meta.reset_index()[["run_id", "label", "strain"]], on="run_id")
    order = run_order(cfg)
    box_palette = [palette(cfg, "strain")[meta.loc[r, "strain"]] for r in order]
    sns.boxplot(data=log2_long, x="run_id", y="log2_intensity", order=order,
                palette=box_palette, ax=ax, fliersize=1.5, linewidth=0.6)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=7)
    ax.set_xlabel("")
    ax.set_ylabel("log2 intensity")
    ax.set_title("Per-run log2 intensity distribution")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F2b_QC_intensity_boxplot.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    # 4. Missingness summary
    detected = (~log2.isna()).astype(int)
    missing_per_protein = 18 - detected.sum(axis=1)
    miss_summary = missing_per_protein.value_counts().sort_index()
    miss_summary.to_csv(ROOT / "outputs" / "tables" / "qc_missingness_per_protein.tsv",
                        sep="\t", header=["n_proteins"])
    print("\n[missingness] proteins missing in N runs:")
    print(miss_summary.to_string())

    # 5. Per-run protein/peptide counts
    counts_per_run = detected.sum(axis=0).to_frame("n_proteins_detected")
    counts_per_run["label"] = [meta.loc[r, "label"] for r in counts_per_run.index]
    counts_per_run.to_csv(ROOT / "outputs" / "tables" / "qc_per_run_detection.tsv", sep="\t")

    print("\n[done] QC outputs in outputs/figures/ and outputs/tables/")


def corr_within(corr: pd.DataFrame, meta: pd.DataFrame) -> float:
    rs = []
    for s, sub in meta.groupby("sample"):
        rids = sub.index.tolist()
        sub_corr = corr.loc[rids, rids].values
        iu = np.triu_indices(len(rids), k=1)
        rs.extend(sub_corr[iu].tolist())
    return float(np.median(rs))


if __name__ == "__main__":
    main()
