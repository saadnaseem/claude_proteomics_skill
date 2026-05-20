"""Functional enrichment via UniProt-GO + KEGG annotation, hypergeometric ORA.

Caches annotation in data/uniprot_annotation.tsv so subsequent runs are fast.
Background = the 1021 quantified proteins (corrects for detection bias).
"""
from __future__ import annotations
import sys
import time
import re
from pathlib import Path
from collections import defaultdict
import requests
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import load_config, ensure_dirs, ROOT


CACHE = ROOT / "data" / "uniprot_annotation.tsv"
KEGG_CACHE = ROOT / "data" / "kegg_pathway_map.tsv"


def clean_accession(acc: str) -> str:
    """DIA-NN concatenates protein-group members with ';'. Use first member."""
    return acc.split(";")[0].strip() if isinstance(acc, str) else acc


def fetch_uniprot(accessions: list[str]) -> pd.DataFrame:
    if CACHE.exists():
        return pd.read_csv(CACHE, sep="\t")
    cleaned = sorted(set(clean_accession(a) for a in accessions if a))
    rows = []
    url = "https://rest.uniprot.org/uniprotkb/accessions"
    fields = "accession,id,gene_names,go,xref_kegg,protein_name,xref_eggnog"
    import io
    for i in range(0, len(cleaned), 100):
        chunk = cleaned[i:i + 100]
        r = requests.get(url, params={"accessions": ",".join(chunk),
                                      "fields": fields, "format": "tsv"})
        r.raise_for_status()
        df = pd.read_csv(io.StringIO(r.text), sep="\t")
        rows.append(df)
        time.sleep(0.3)
        print(f"  fetched {min(i+100,len(cleaned))}/{len(cleaned)}")
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(CACHE, sep="\t", index=False)
    return out


def fetch_kegg_pathway_map() -> dict[str, list[str]]:
    """Map ppu KEGG gene IDs (e.g. ppu:PP_2679) -> list of pathway IDs."""
    if KEGG_CACHE.exists():
        df = pd.read_csv(KEGG_CACHE, sep="\t")
    else:
        r = requests.get("https://rest.kegg.jp/link/pathway/ppu")
        r.raise_for_status()
        rows = [line.split("\t") for line in r.text.strip().splitlines()]
        df = pd.DataFrame(rows, columns=["gene", "pathway"])
        df.to_csv(KEGG_CACHE, sep="\t", index=False)
    m = defaultdict(list)
    for _, row in df.iterrows():
        m[row["gene"]].append(row["pathway"])

    # Pathway name lookup
    nm_path = ROOT / "data" / "kegg_pathway_names.tsv"
    if nm_path.exists():
        names = pd.read_csv(nm_path, sep="\t")
    else:
        r = requests.get("https://rest.kegg.jp/list/pathway/ppu")
        rows = [line.split("\t") for line in r.text.strip().splitlines()]
        names = pd.DataFrame(rows, columns=["pathway", "name"])
        names.to_csv(nm_path, sep="\t", index=False)
    # KEGG list returns ids like "ppu01100"; link returns "path:ppu01100" → normalize
    name_map = {}
    for pw, nm in zip(names["pathway"], names["name"]):
        nm_clean = nm.split(" - Pseudomonas putida")[0] if isinstance(nm, str) else nm
        key1 = pw if pw.startswith("path:") else f"path:{pw}"
        key2 = pw.replace("path:", "")
        name_map[key1] = nm_clean
        name_map[key2] = nm_clean
    return dict(m), name_map


def parse_go_field(go_text: str) -> list[tuple[str, str, str]]:
    """Return list of (go_id, name, ontology_guess) — ontology_guess inferred."""
    if not isinstance(go_text, str) or not go_text.strip():
        return []
    out = []
    for term in go_text.split(";"):
        term = term.strip()
        if not term:
            continue
        m = re.match(r"(.*?) \[(GO:\d+)\]\s*$", term)
        if m:
            out.append((m.group(2), m.group(1), ""))
    return out


def fetch_go_ontology_map(go_ids: set[str]) -> dict[str, str]:
    """Fetch ontology (BP/MF/CC) for GO IDs via QuickGO (cached)."""
    cache = ROOT / "data" / "go_ontology_map.tsv"
    if cache.exists():
        df = pd.read_csv(cache, sep="\t")
        existing = dict(zip(df["go_id"], df["ontology"]))
        missing = [g for g in go_ids if g not in existing]
        if not missing:
            return existing
    else:
        existing = {}
        missing = list(go_ids)

    print(f"  fetching ontology for {len(missing)} new GO IDs from QuickGO ...")
    url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{ids}"
    rows = []
    for i in range(0, len(missing), 200):
        chunk = missing[i:i + 200]
        r = requests.get(url.format(ids=",".join(chunk)),
                         headers={"Accept": "application/json"})
        if r.ok:
            data = r.json()
            for entry in data.get("results", []):
                aspect = entry.get("aspect")
                ont = {"biological_process": "BP",
                       "molecular_function": "MF",
                       "cellular_component": "CC"}.get(aspect, "?")
                rows.append((entry["id"], ont))
        time.sleep(0.3)
    new_df = pd.DataFrame(rows, columns=["go_id", "ontology"])
    if cache.exists():
        old = pd.read_csv(cache, sep="\t")
        new_df = pd.concat([old, new_df]).drop_duplicates("go_id")
    new_df.to_csv(cache, sep="\t", index=False)
    return dict(zip(new_df["go_id"], new_df["ontology"]))


def hypergeom_ora(hits: set[str], background: set[str],
                  term_to_genes: dict[str, set[str]],
                  term_names: dict[str, str],
                  ontology_label: str = "GO") -> pd.DataFrame:
    M = len(background)
    N = len(hits & background)
    rows = []
    for term, genes in term_to_genes.items():
        K = len(genes & background)
        x = len(genes & hits & background)
        if x < 2 or K < 2:
            continue
        # Probability of seeing >= x hits (sf gives P(X > x-1))
        p = hypergeom.sf(x - 1, M, K, N)
        rows.append({"term": term, "name": term_names.get(term, ""),
                     "ontology": ontology_label,
                     "k": x, "K": K, "n": N, "M": M,
                     "fold_enrichment": (x / N) / (K / M) if N and K else np.nan,
                     "p": p})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    _, q, _, _ = multipletests(df["p"].values, method="fdr_bh")
    df["q"] = q
    return df.sort_values("q")


def main() -> None:
    cfg = load_config()
    ensure_dirs(cfg)

    log2 = pd.read_csv(ROOT / "outputs" / "tables" / "protein_log2_wide.tsv",
                       sep="\t", index_col=[0, 1, 2])
    raw_pg = log2.index.get_level_values(0).unique().tolist()
    background_acc = sorted({clean_accession(a) for a in raw_pg if a})
    # raw_to_clean: for mapping master DEPs (which use Protein.Group as raw) onto cleaned accessions
    raw_to_clean = {a: clean_accession(a) for a in raw_pg}
    print(f"[enrichment] background = {len(background_acc)} quantified proteins")

    print("[enrichment] fetching UniProt annotations …")
    annot = fetch_uniprot(background_acc)

    # Build GO id → set(accession) and GO id → name (parse_go_field gives name in [text])
    go_to_acc = defaultdict(set)
    go_names = {}
    for _, row in annot.iterrows():
        for go_id, name, _ in parse_go_field(row["Gene Ontology (GO)"]):
            go_to_acc[go_id].add(row["Entry"])
            go_names[go_id] = name
    print(f"[enrichment] {len(go_to_acc)} GO terms in background")

    # Determine ontology (BP/MF/CC) per GO id
    ont_map = fetch_go_ontology_map(set(go_to_acc.keys()))
    by_ont = {"BP": {}, "MF": {}, "CC": {}}
    for go_id, accs in go_to_acc.items():
        ont = ont_map.get(go_id, "?")
        if ont in by_ont:
            by_ont[ont][go_id] = accs

    # KEGG pathway annotation
    print("[enrichment] fetching KEGG pathway map …")
    kegg_map, kegg_names = fetch_kegg_pathway_map()
    # Build accession -> KEGG gene id map from UniProt 'KEGG' column ("ppu:PP_xxxx;")
    kegg_id_for_acc = {}
    for _, row in annot.iterrows():
        v = str(row.get("KEGG", "")).rstrip(";")
        if v.startswith("ppu:"):
            kegg_id_for_acc[row["Entry"]] = v
    pathway_to_acc = defaultdict(set)
    for acc, gid in kegg_id_for_acc.items():
        for pw in kegg_map.get(gid, []):
            pathway_to_acc[pw].add(acc)
    print(f"[enrichment] {len(pathway_to_acc)} KEGG pathways in background")

    # COG categories from eggNOG (COGxxxx → category letters)
    # Skip for now — would need a separate COG-to-category lookup. Note in plan.

    # Run ORA for each contrast (UP and DOWN separately)
    master = pd.read_csv(ROOT / "outputs" / "tables" / "all_contrasts_long.tsv", sep="\t")
    background_set = set(background_acc)
    summary_rows = []
    for contrast in master["contrast"].unique():
        sub = master[master["contrast"] == contrast]
        for direction in ["UP", "DOWN"]:
            hits = {clean_accession(a) for a in sub.loc[sub["sig_dir"] == direction, "Protein.Group"]}
            if len(hits) < 3:
                continue
            for ont_label, term_dict in by_ont.items():
                df = hypergeom_ora(hits, background_set, term_dict, go_names,
                                   f"GO_{ont_label}")
                if not df.empty:
                    df["contrast"] = contrast
                    df["direction"] = direction
                    df.to_csv(ROOT / "outputs" / "enrichment" /
                              f"GO{ont_label}_{contrast}_{direction}.tsv",
                              sep="\t", index=False)
                    sig = (df["q"] <= 0.10).sum()
                    summary_rows.append({"contrast": contrast, "direction": direction,
                                         "ontology": f"GO_{ont_label}",
                                         "n_terms_q≤0.10": int(sig)})

            kdf = hypergeom_ora(hits, background_set, pathway_to_acc, kegg_names,
                                "KEGG")
            if not kdf.empty:
                kdf["contrast"] = contrast
                kdf["direction"] = direction
                kdf.to_csv(ROOT / "outputs" / "enrichment" /
                           f"KEGG_{contrast}_{direction}.tsv",
                           sep="\t", index=False)
                sig = (kdf["q"] <= 0.10).sum()
                summary_rows.append({"contrast": contrast, "direction": direction,
                                     "ontology": "KEGG",
                                     "n_terms_q≤0.10": int(sig)})

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(ROOT / "outputs" / "enrichment" / "enrichment_summary.tsv",
                   sep="\t", index=False)
    print("\n[enrichment summary, top per contrast]")
    print(summary.to_string(index=False))

    # Dotplots: top 10 BP terms per contrast, faceted UP vs DOWN
    plot_dotplots(by_ont, master, background_set, go_names, kegg_names,
                  pathway_to_acc, cfg)
    print("\n[done]")


def plot_dotplots(by_ont, master, background_set, go_names, kegg_names,
                  pathway_to_acc, cfg):
    contrasts = ["KT_Chlys_vs_glucose", "KT_butamine_vs_glucose",
                 "HGL1175_Chlys_vs_glucose", "HGL1175_butamine_vs_glucose",
                 "glucose_HGL1175_vs_KT", "Chlys_HGL1175_vs_KT",
                 "Butamine_HGL1175_vs_KT"]
    for direction in ["UP", "DOWN"]:
        rows_all = []
        for c in contrasts:
            hits = {clean_accession(a) for a in master.loc[(master["contrast"] == c) &
                                  (master["sig_dir"] == direction),
                                  "Protein.Group"]}
            if len(hits) < 3: continue
            df = hypergeom_ora(hits, background_set, by_ont["BP"], go_names, "GO_BP")
            if df.empty: continue
            df = df.head(8).assign(contrast=c)
            rows_all.append(df)
        if not rows_all: continue
        all_df = pd.concat(rows_all)
        if all_df.empty: continue
        fig, ax = plt.subplots(figsize=(11, max(4, 0.35 * all_df["name"].nunique())))
        terms = all_df.groupby("name")["q"].min().sort_values().index.tolist()
        cs = [c for c in contrasts if c in all_df["contrast"].values]
        for c in cs:
            sub = all_df[all_df["contrast"] == c]
            xs = [cs.index(c)] * len(sub)
            ys = [terms.index(t) for t in sub["name"]]
            sizes = sub["k"].values * 30
            colors = -np.log10(sub["q"].clip(lower=1e-12).values)
            sc = ax.scatter(xs, ys, s=sizes, c=colors, cmap="viridis_r",
                            vmin=1, vmax=4, edgecolor="black", lw=0.4)
        cb = plt.colorbar(sc, ax=ax)
        cb.set_label("-log10 q")
        ax.set_xticks(range(len(cs))); ax.set_xticklabels(cs, rotation=30, ha="right")
        ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms, fontsize=8)
        ax.set_title(f"GO_BP enrichment — {direction} DEPs (top 8 per contrast)")
        ax.grid(alpha=0.2)
        plt.tight_layout()
        plt.savefig(ROOT / "outputs" / "figures" /
                    f"F9_enrichment_GOBP_{direction}.png",
                    dpi=cfg["plot"]["dpi"], bbox_inches="tight")
        plt.close()


if __name__ == "__main__":
    main()
