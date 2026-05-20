"""Per-protein two-factor ANOVA (strain x medium) with BH-correction on each term."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import load_config, metadata_df, ensure_dirs, ROOT, run_order


def main() -> None:
    cfg = load_config()
    ensure_dirs(cfg)
    meta = metadata_df(cfg).set_index("run_id").loc[run_order(cfg)]
    log2 = pd.read_csv(ROOT / "outputs" / "tables" / "protein_log2_wide.tsv",
                       sep="\t", index_col=[0, 1, 2])

    rows = []
    runs = log2.columns.tolist()
    fit_meta = meta.reset_index()

    for prot_key, intensities in log2.iterrows():
        df = fit_meta.copy()
        df["y"] = intensities.values
        df = df.dropna(subset=["y"])
        if df["y"].nunique() < 4 or len(df) < 8:
            continue
        try:
            model = smf.ols("y ~ C(strain) * C(medium)", data=df).fit()
            anova = sm.stats.anova_lm(model, typ=2)
        except Exception:
            continue

        f_strain = anova.loc["C(strain)",            "F"]
        p_strain = anova.loc["C(strain)",            "PR(>F)"]
        f_medium = anova.loc["C(medium)",            "F"]
        p_medium = anova.loc["C(medium)",            "PR(>F)"]
        f_inter  = anova.loc["C(strain):C(medium)",  "F"]
        p_inter  = anova.loc["C(strain):C(medium)",  "PR(>F)"]

        rows.append({
            "Protein.Group": prot_key[0],
            "Protein":       prot_key[1],
            "Description":   prot_key[2],
            "F_strain":      f_strain, "p_strain":      p_strain,
            "F_medium":      f_medium, "p_medium":      p_medium,
            "F_interaction": f_inter,  "p_interaction": p_inter,
        })

    out = pd.DataFrame(rows)
    for term in ["strain", "medium", "interaction"]:
        valid = out[f"p_{term}"].notna()
        q = np.full(len(out), np.nan)
        if valid.sum() > 0:
            _, q_valid, _, _ = multipletests(out.loc[valid, f"p_{term}"].values, method="fdr_bh")
            q[valid.values] = q_valid
        out[f"q_{term}"] = q

    out.to_csv(ROOT / "outputs" / "tables" / "two_factor_anova.tsv", sep="\t", index=False)
    print(f"[ANOVA] {len(out)} proteins fit")
    for term in ["strain", "medium", "interaction"]:
        n_sig = (out[f"q_{term}"] <= 0.05).sum()
        print(f"  q_{term} ≤ 0.05: {n_sig} proteins")

    # Volcano-style plots: -log10(q) vs F-statistic for each term
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, term, color in zip(axes, ["strain", "medium", "interaction"],
                                ["#4477AA", "#228833", "#EE6677"]):
        sub = out.dropna(subset=[f"q_{term}"]).copy()
        sub["nl10q"] = -np.log10(sub[f"q_{term}"].clip(lower=1e-12))
        sub["log10F"] = np.log10(sub[f"F_{term}"].clip(lower=1e-3))
        sig = sub[f"q_{term}"] <= 0.05
        ax.scatter(sub.loc[~sig, "log10F"], sub.loc[~sig, "nl10q"],
                   s=8, color="#cccccc", alpha=0.5)
        ax.scatter(sub.loc[ sig, "log10F"], sub.loc[ sig, "nl10q"],
                   s=12, color=color, alpha=0.9)
        ax.axhline(-np.log10(0.05), ls=":", color="gray")
        ax.set_xlabel("log10 F-statistic")
        ax.set_ylabel("-log10 q (BH)")
        ax.set_title(f"{term} (q≤0.05: {sig.sum()})")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F3b_anova_terms.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    # Save interaction-significant subset
    inter_sig = out[(out["q_interaction"] <= 0.05)].sort_values("q_interaction")
    inter_sig.to_csv(ROOT / "outputs" / "tables" / "interaction_significant_proteins.tsv",
                     sep="\t", index=False)
    print(f"\n[interaction q≤0.05] {len(inter_sig)} proteins → interaction_significant_proteins.tsv")

    print("[done]")


if __name__ == "__main__":
    main()
