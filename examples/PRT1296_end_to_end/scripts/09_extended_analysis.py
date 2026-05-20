"""Extended analyses (Phase 2):
  (a) per-cluster enrichment for the 6 heatmap clusters
  (b) per-DEP-set enrichment (KT-Chlys-only, KT-butyl-only, KT-shared, HGL-only)
  (c) stress-amplitude scatter (HGL amplitude vs KT amplitude per protein)
  (d) annotated volcano panels with gene labels for top hits
"""
from __future__ import annotations
import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import load_config, ensure_dirs, ROOT, palette

# Reuse helpers from 06_enrichment
import importlib.util
spec = importlib.util.spec_from_file_location("enrichment_module",
                                              ROOT / "scripts" / "06_enrichment.py")
em = importlib.util.module_from_spec(spec); spec.loader.exec_module(em)


def main() -> None:
    cfg = load_config()
    ensure_dirs(cfg)

    annot = pd.read_csv(ROOT / "data" / "uniprot_annotation.tsv", sep="\t")
    log2  = pd.read_csv(ROOT / "outputs" / "tables" / "protein_log2_wide.tsv",
                         sep="\t", index_col=[0, 1, 2])
    raw_pg = log2.index.get_level_values(0).unique().tolist()
    background = sorted({em.clean_accession(a) for a in raw_pg if a})

    # GO term sets
    go_to_acc = defaultdict(set)
    go_names = {}
    for _, row in annot.iterrows():
        for go_id, name, _ in em.parse_go_field(row["Gene Ontology (GO)"]):
            go_to_acc[go_id].add(row["Entry"])
            go_names[go_id] = name
    ont_map = em.fetch_go_ontology_map(set(go_to_acc.keys()))
    go_BP = {gid: accs for gid, accs in go_to_acc.items() if ont_map.get(gid) == "BP"}
    go_MF = {gid: accs for gid, accs in go_to_acc.items() if ont_map.get(gid) == "MF"}

    kegg_map, kegg_names = em.fetch_kegg_pathway_map()
    kegg_id_for_acc = {row["Entry"]: str(row["KEGG"]).rstrip(";")
                       for _, row in annot.iterrows()
                       if str(row.get("KEGG", "")).startswith("ppu:")}
    pathway_to_acc = defaultdict(set)
    for acc, gid in kegg_id_for_acc.items():
        for pw in kegg_map.get(gid, []):
            pathway_to_acc[pw].add(acc)

    bg_set = set(background)

    # ----- (a) per-cluster enrichment -----
    cluster_df = pd.read_csv(ROOT / "outputs" / "tables" / "heatmap_cluster_membership.tsv",
                              sep="\t")
    cluster_summary = []
    for cl in sorted(cluster_df["cluster"].unique()):
        members = cluster_df[cluster_df.cluster == cl]
        hits = {em.clean_accession(a) for a in members["Protein.Group"]}
        for ont_label, term_dict in [("GO_BP", go_BP), ("GO_MF", go_MF), ("KEGG", pathway_to_acc)]:
            term_names = go_names if ont_label.startswith("GO") else kegg_names
            df = em.hypergeom_ora(hits, bg_set, term_dict, term_names, ont_label)
            if df.empty: continue
            df["cluster"] = cl
            df.to_csv(ROOT / "outputs" / "enrichment" /
                      f"cluster{cl}_{ont_label}.tsv", sep="\t", index=False)
            top = df.head(3)
            for _, r in top.iterrows():
                cluster_summary.append({
                    "cluster": cl, "n_members": len(members),
                    "ontology": ont_label,
                    "term": r["term"], "name": r["name"],
                    "k": int(r["k"]), "K": int(r["K"]),
                    "fold": round(float(r["fold_enrichment"]), 2),
                    "q": float(r["q"])})
    pd.DataFrame(cluster_summary).to_csv(
        ROOT / "outputs" / "enrichment" / "cluster_top_terms.tsv", sep="\t", index=False)
    print(f"[clusters] {cluster_df['cluster'].nunique()} clusters annotated")

    # ----- (b) per-DEP-set enrichment -----
    master = pd.read_csv(ROOT / "outputs" / "tables" / "all_contrasts_long.tsv", sep="\t")
    kt_chlys = set(master[(master.contrast == "KT_Chlys_vs_glucose") &
                          (master.sig_dir == "UP")]["Protein.Group"])
    kt_buty  = set(master[(master.contrast == "KT_butamine_vs_glucose") &
                          (master.sig_dir == "UP")]["Protein.Group"])
    sets = {
        "KT_chlys_only":   kt_chlys - kt_buty,
        "KT_butyl_only":   kt_buty - kt_chlys,
        "KT_shared":       kt_chlys & kt_buty,
        "HGL_const_up":    set(master[(master.contrast == "glucose_HGL1175_vs_KT") &
                                      (master.sig_dir == "UP")]["Protein.Group"]),
        "HGL_const_down":  set(master[(master.contrast == "glucose_HGL1175_vs_KT") &
                                      (master.sig_dir == "DOWN")]["Protein.Group"]),
    }
    set_summary = []
    for set_name, members in sets.items():
        hits = {em.clean_accession(a) for a in members}
        if len(hits) < 3: continue
        for ont_label, term_dict in [("GO_BP", go_BP), ("KEGG", pathway_to_acc)]:
            term_names = go_names if ont_label.startswith("GO") else kegg_names
            df = em.hypergeom_ora(hits, bg_set, term_dict, term_names, ont_label)
            if df.empty: continue
            df["set"] = set_name
            df.to_csv(ROOT / "outputs" / "enrichment" /
                      f"set_{set_name}_{ont_label}.tsv", sep="\t", index=False)
            for _, r in df.head(3).iterrows():
                set_summary.append({
                    "set": set_name, "n_members": len(hits),
                    "ontology": ont_label,
                    "term": r["term"], "name": r["name"],
                    "k": int(r["k"]), "K": int(r["K"]),
                    "fold": round(float(r["fold_enrichment"]), 2),
                    "q": float(r["q"])})
    pd.DataFrame(set_summary).to_csv(
        ROOT / "outputs" / "enrichment" / "set_top_terms.tsv", sep="\t", index=False)
    print(f"[sets] {len(sets)} DEP sets annotated")

    # ----- (c) stress-amplitude scatter -----
    plot_stress_amplitude(master, cfg)

    # ----- (d) annotated volcano panels -----
    plot_volcanoes(master, cfg)
    print("[done]")


def plot_stress_amplitude(master: pd.DataFrame, cfg: dict) -> None:
    """Per-protein |log2FC| in HGL stress-vs-glucose vs |log2FC| in KT stress-vs-glucose."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4))
    for ax, medium in zip(axes, ["Chlys", "butamine"]):
        kt = master[master.contrast == f"KT_{medium}_vs_glucose"][
                    ["Protein", "Protein.Description", "log2FC", "q"]].set_index("Protein")
        hg = master[master.contrast == f"HGL1175_{medium}_vs_glucose"][
                    ["Protein", "log2FC", "q"]].set_index("Protein")
        merged = kt.join(hg, lsuffix="_KT", rsuffix="_HGL", how="inner")
        merged["sig_KT"]  = merged["q_KT"]  <= 0.05
        merged["sig_HGL"] = merged["q_HGL"] <= 0.05

        ax.scatter(merged.loc[~merged.sig_KT & ~merged.sig_HGL, "log2FC_KT"],
                   merged.loc[~merged.sig_KT & ~merged.sig_HGL, "log2FC_HGL"],
                   s=6, color="#bbbbbb", alpha=0.5, label="ns")
        ax.scatter(merged.loc[merged.sig_KT & ~merged.sig_HGL, "log2FC_KT"],
                   merged.loc[merged.sig_KT & ~merged.sig_HGL, "log2FC_HGL"],
                   s=20, color=palette(cfg, "strain")["KT2440"], alpha=0.85,
                   label="KT-only DEP")
        ax.scatter(merged.loc[~merged.sig_KT & merged.sig_HGL, "log2FC_KT"],
                   merged.loc[~merged.sig_KT & merged.sig_HGL, "log2FC_HGL"],
                   s=20, color=palette(cfg, "strain")["HGL1175"], alpha=0.85,
                   label="HGL-only DEP")
        ax.scatter(merged.loc[merged.sig_KT & merged.sig_HGL, "log2FC_KT"],
                   merged.loc[merged.sig_KT & merged.sig_HGL, "log2FC_HGL"],
                   s=22, color="black", alpha=0.85, label="both")

        # y=x reference
        lim = max(merged["log2FC_KT"].abs().max(), merged["log2FC_HGL"].abs().max())
        lim = float(np.ceil(lim))
        ax.plot([-lim, lim], [-lim, lim], ls="--", color="#cccccc", lw=1, zorder=0)
        ax.axhline(0, lw=0.5, color="gray"); ax.axvline(0, lw=0.5, color="gray")
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xlabel(f"log2FC  KT   ({medium} / glucose)")
        ax.set_ylabel(f"log2FC  HGL  ({medium} / glucose)")
        ax.set_title(f"Stress-response amplitude: HGL vs KT in {medium}")
        ax.grid(alpha=0.2)

        # Label top KT-DEPs
        kt_top = merged[merged.sig_KT].nlargest(8, "log2FC_KT")
        for prot, r in kt_top.iterrows():
            ax.annotate(prot, (r["log2FC_KT"], r["log2FC_HGL"]),
                        fontsize=7, alpha=0.85,
                        xytext=(3, 3), textcoords="offset points")
    axes[0].legend(loc="upper left", fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F13_stress_amplitude.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()


def plot_volcanoes(master: pd.DataFrame, cfg: dict) -> None:
    """4-panel volcano: KT-Chlys, KT-Bu, HGL-Chlys, HGL-Bu  vs glucose. Annotate top genes."""
    contrasts = [
        ("KT_Chlys_vs_glucose",          "KT2440 / Chlys"),
        ("KT_butamine_vs_glucose",       "KT2440 / butylamine"),
        ("HGL1175_Chlys_vs_glucose",     "HGL1175 / Chlys"),
        ("HGL1175_butamine_vs_glucose",  "HGL1175 / butylamine"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9))
    for ax, (name, title) in zip(axes.flat, contrasts):
        sub = master[master.contrast == name].copy()
        sub["nl10q"] = -np.log10(sub["q"].clip(lower=1e-12))
        ax.scatter(sub.loc[sub.sig_dir == "ns", "log2FC"],
                   sub.loc[sub.sig_dir == "ns", "nl10q"],
                   s=6, color="#cccccc", alpha=0.6)
        ax.scatter(sub.loc[sub.sig_dir == "UP", "log2FC"],
                   sub.loc[sub.sig_dir == "UP", "nl10q"],
                   s=22, color="#EE6677", alpha=0.85, label="UP")
        ax.scatter(sub.loc[sub.sig_dir == "DOWN", "log2FC"],
                   sub.loc[sub.sig_dir == "DOWN", "nl10q"],
                   s=22, color="#4477AA", alpha=0.85, label="DOWN")
        ax.axhline(-np.log10(0.05), ls=":", color="gray")
        ax.axvline( 1, ls=":", color="gray"); ax.axvline(-1, ls=":", color="gray")
        ax.set_xlabel("log2FC"); ax.set_ylabel("-log10 q")
        n_up = (sub.sig_dir == "UP").sum(); n_dn = (sub.sig_dir == "DOWN").sum()
        ax.set_title(f"{title}  (UP={n_up}, DOWN={n_dn})", fontsize=10)
        # Label top 6 by q
        top = sub[sub.sig_dir != "ns"].nsmallest(6, "q")
        for _, r in top.iterrows():
            ax.annotate(r["Protein"], (r["log2FC"], r["nl10q"]),
                        fontsize=7.5, alpha=0.95,
                        xytext=(3, 3), textcoords="offset points")
        ax.grid(alpha=0.2); ax.legend(fontsize=8, loc="upper left")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F4_volcanoes_hydrolysate_response.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()

    # Strain-effect volcano
    contrasts2 = [
        ("glucose_HGL1175_vs_KT",      "HGL1175 vs KT in glucose"),
        ("Chlys_HGL1175_vs_KT",        "HGL1175 vs KT in Chlys"),
        ("Butamine_HGL1175_vs_KT",     "HGL1175 vs KT in butylamine"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, (name, title) in zip(axes, contrasts2):
        sub = master[master.contrast == name].copy()
        sub["nl10q"] = -np.log10(sub["q"].clip(lower=1e-12))
        ax.scatter(sub.loc[sub.sig_dir == "ns", "log2FC"],
                   sub.loc[sub.sig_dir == "ns", "nl10q"],
                   s=6, color="#cccccc", alpha=0.6)
        ax.scatter(sub.loc[sub.sig_dir == "UP", "log2FC"],
                   sub.loc[sub.sig_dir == "UP", "nl10q"],
                   s=22, color="#EE6677", alpha=0.85, label="UP in HGL")
        ax.scatter(sub.loc[sub.sig_dir == "DOWN", "log2FC"],
                   sub.loc[sub.sig_dir == "DOWN", "nl10q"],
                   s=22, color="#4477AA", alpha=0.85, label="DOWN in HGL")
        ax.axhline(-np.log10(0.05), ls=":", color="gray")
        ax.axvline( 1, ls=":", color="gray"); ax.axvline(-1, ls=":", color="gray")
        ax.set_xlabel("log2FC (HGL / KT)"); ax.set_ylabel("-log10 q")
        n_up = (sub.sig_dir == "UP").sum(); n_dn = (sub.sig_dir == "DOWN").sum()
        ax.set_title(f"{title}\n(UP={n_up}, DOWN={n_dn})", fontsize=10)
        top = sub[sub.sig_dir != "ns"].nsmallest(8, "q")
        for _, r in top.iterrows():
            ax.annotate(r["Protein"], (r["log2FC"], r["nl10q"]),
                        fontsize=7.5, alpha=0.95,
                        xytext=(3, 3), textcoords="offset points")
        ax.grid(alpha=0.2); ax.legend(fontsize=8, loc="upper left")
    plt.tight_layout()
    plt.savefig(ROOT / "outputs" / "figures" / "F5_volcanoes_strain_effect.png",
                dpi=cfg["plot"]["dpi"], bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
