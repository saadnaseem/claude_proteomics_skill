"""Aggregate the 11 existing pairwise t-tests into a master DEP table and draw UpSet plots."""
from __future__ import annotations
import sys
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import load_config, ensure_dirs, ROOT


def parse_t_test_xlsx(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Full t-test output")
    cols = list(df.columns)
    fc_col = [c for c in cols if "log2_Fold_change" in c][0]
    p_col  = "p-value"
    q_col  = "p_adjusted(BH)"
    out = df[["Protein", "Protein.Group", "Protein.Names", "Protein.Description",
              fc_col, p_col, q_col]].copy()
    out.columns = ["Protein", "Protein.Group", "Protein.Names", "Protein.Description",
                   "log2FC", "p", "q"]
    return out


def main() -> None:
    cfg = load_config()
    ensure_dirs(cfg)
    q_cut  = cfg["thresholds"]["q_strict"]
    fc_cut = cfg["thresholds"]["log2fc_strict"]

    t_dir = ROOT / "data" / "t_test_xlsx"
    files = sorted(t_dir.glob("t-test_*.xlsx"))
    long_rows = []
    for f in files:
        name = re.sub(r"^t-test_|_\d{8}-\d+\.xlsx$", "", f.name)
        df = parse_t_test_xlsx(f).assign(contrast=name)
        df["sig"]      = (df["q"] <= q_cut) & (df["log2FC"].abs() >= fc_cut)
        df["sig_dir"]  = np.where(df["sig"], np.where(df["log2FC"] > 0, "UP", "DOWN"), "ns")
        long_rows.append(df)

    master = pd.concat(long_rows, ignore_index=True)
    master.to_csv(ROOT / "outputs" / "tables" / "all_contrasts_long.tsv",
                  sep="\t", index=False)
    summary = (master.groupby(["contrast", "sig_dir"])
                     .size()
                     .unstack(fill_value=0)
                     .rename_axis(columns=None))
    summary["total_DEPs"] = summary.get("UP", 0) + summary.get("DOWN", 0)
    summary.to_csv(ROOT / "outputs" / "tables" / "dep_counts_per_contrast.tsv", sep="\t")
    print("\n[DEP counts per contrast]")
    print(summary.to_string())

    # UpSet for hydrolysate-vs-glucose contrasts
    hydro_contrasts = ["KT_Chlys_vs_glucose", "KT_butamine_vs_glucose",
                       "HGL1175_Chlys_vs_glucose", "HGL1175_butamine_vs_glucose"]
    for direction in ["UP", "DOWN"]:
        ind = build_indicator(master, hydro_contrasts, direction)
        if ind.empty or ind[hydro_contrasts].sum().sum() == 0:
            print(f"[upset] no {direction} DEPs across hydrolysate contrasts")
            continue
        plot_upset(ind, hydro_contrasts,
                   ROOT / "outputs" / "figures" / f"F6_upset_hydrolysate_{direction}.png",
                   title=f"Hydrolysate vs glucose response — {direction} DEPs")
        ind.to_csv(ROOT / "outputs" / "tables" / f"upset_hydrolysate_{direction}.tsv",
                   sep="\t", index=False)

    # UpSet for strain-effect contrasts
    strain_contrasts = ["glucose_HGL1175_vs_KT", "Chlys_HGL1175_vs_KT", "Butamine_HGL1175_vs_KT"]
    for direction in ["UP", "DOWN"]:
        ind = build_indicator(master, strain_contrasts, direction)
        if ind.empty or ind[strain_contrasts].sum().sum() == 0:
            continue
        plot_upset(ind, strain_contrasts,
                   ROOT / "outputs" / "figures" / f"F7_upset_strain_{direction}.png",
                   title=f"HGL1175 vs KT2440 across media — {direction} DEPs")
        ind.to_csv(ROOT / "outputs" / "tables" / f"upset_strain_{direction}.tsv",
                   sep="\t", index=False)

    # Set assignment table: per-protein membership across the 4 hydrolysate contrasts
    assign = (master[master["contrast"].isin(hydro_contrasts)]
                    .pivot_table(index=["Protein.Group", "Protein", "Protein.Description"],
                                 columns="contrast", values="sig_dir",
                                 aggfunc="first"))
    assign.to_csv(ROOT / "outputs" / "tables" / "dep_set_assignment_hydrolysate.tsv",
                  sep="\t")

    print("\n[done]")


def build_indicator(master: pd.DataFrame, contrasts: list, direction: str) -> pd.DataFrame:
    sub = master[master["contrast"].isin(contrasts)].copy()
    sub["hit"] = sub["sig_dir"] == direction
    pivot = sub.pivot_table(index=["Protein.Group", "Protein", "Protein.Description"],
                            columns="contrast", values="hit", aggfunc="any", fill_value=False)
    pivot = pivot.reindex(columns=contrasts, fill_value=False)
    return pivot[pivot.any(axis=1)].reset_index()


def plot_upset(ind: pd.DataFrame, contrasts: list, out: Path, title: str) -> None:
    """Self-contained UpSet plot: intersection bar chart + dot matrix + per-set bars."""
    bool_only = ind[contrasts].astype(bool)
    if bool_only.shape[0] < 2:
        return

    # Group rows by their boolean signature; count members of each intersection
    sig = bool_only.apply(lambda r: tuple(r.values), axis=1)
    counts = sig.value_counts()
    counts = counts[[s for s in counts.index if any(s)]]
    counts = counts.sort_values(ascending=False)
    if counts.empty:
        return

    n_sets = len(contrasts)
    n_intersect = len(counts)
    set_totals = bool_only.sum(axis=0)

    fig = plt.figure(figsize=(max(8, 0.6 * n_intersect + 4), 5.8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1.4, max(2.5, 0.4 * n_intersect)],
                           height_ratios=[3, 2], hspace=0.05, wspace=0.05)

    # Top-right: intersection size bars
    ax_top = fig.add_subplot(gs[0, 1])
    x = np.arange(n_intersect)
    ax_top.bar(x, counts.values, color="#333333", width=0.7)
    for xi, v in zip(x, counts.values):
        ax_top.text(xi, v + 0.02 * counts.max(), str(int(v)),
                    ha="center", va="bottom", fontsize=8)
    ax_top.set_xticks([]); ax_top.set_xlim(-0.5, n_intersect - 0.5)
    ax_top.set_ylabel("intersection size", fontsize=9)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    # Bottom-right: dot matrix
    ax_dot = fig.add_subplot(gs[1, 1], sharex=ax_top)
    for j, set_name in enumerate(contrasts):
        for i, sig_tuple in enumerate(counts.index):
            on = sig_tuple[j]
            ax_dot.scatter(i, n_sets - 1 - j,
                           s=140, color=("#222222" if on else "#dddddd"),
                           edgecolor="none", zorder=3)
        # connecting line for "on" cells in each intersection
    for i, sig_tuple in enumerate(counts.index):
        ys = [n_sets - 1 - j for j, on in enumerate(sig_tuple) if on]
        if len(ys) >= 2:
            ax_dot.plot([i, i], [min(ys), max(ys)], color="#222222", lw=2, zorder=2)
    ax_dot.set_xticks(x); ax_dot.set_xticklabels([])
    ax_dot.set_yticks(range(n_sets))
    ax_dot.set_yticklabels([contrasts[n_sets - 1 - i] for i in range(n_sets)], fontsize=9)
    ax_dot.set_xlim(-0.5, n_intersect - 0.5)
    ax_dot.set_ylim(-0.6, n_sets - 0.4)
    for s in ["top", "right", "bottom"]:
        ax_dot.spines[s].set_visible(False)
    ax_dot.tick_params(axis="x", which="both", bottom=False)

    # Bottom-left: per-set total bars (horizontal, mirrored)
    ax_left = fig.add_subplot(gs[1, 0], sharey=ax_dot)
    ys = [n_sets - 1 - j for j in range(n_sets)]
    totals = [set_totals[contrasts[j]] for j in range(n_sets)]
    ax_left.barh(ys, totals, color="#888888", height=0.6)
    for y, v in zip(ys, totals):
        ax_left.text(-0.02 * max(totals), y, str(int(v)),
                     ha="right", va="center", fontsize=8)
    ax_left.invert_xaxis()
    ax_left.set_xlabel("set size", fontsize=9)
    ax_left.set_yticks([]); ax_left.set_ylim(-0.6, n_sets - 0.4)
    for s in ["top", "right", "left"]:
        ax_left.spines[s].set_visible(False)

    fig.suptitle(title, fontsize=11)
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
