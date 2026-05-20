"""Build sample metadata, load proteins, write authoritative wide matrix."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (load_config, metadata_df, load_protein_long, load_top3,
                     long_to_wide, log2_safe, ensure_dirs, ROOT, run_order)


def main() -> None:
    cfg = load_config()
    ensure_dirs(cfg)

    meta = metadata_df(cfg)
    meta_path = ROOT / "outputs" / "tables" / "sample_metadata.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)

    long = load_protein_long(cfg)
    print(f"[long] {len(long)} rows, {long['Protein.Group'].nunique()} unique proteins")

    wide = long_to_wide(long, "Counts_sum")
    wide = wide.reindex(columns=run_order(cfg))
    log2_wide = log2_safe(wide)

    wide.to_csv(ROOT / "outputs" / "tables" / "protein_counts_wide.tsv", sep="\t")
    log2_wide.to_csv(ROOT / "outputs" / "tables" / "protein_log2_wide.tsv", sep="\t")
    print(f"[wide] {wide.shape[0]} proteins × {wide.shape[1]} runs")

    top3 = load_top3(cfg)
    top3["run_id"] = top3["Sample"] + "-" + top3["Replicate"]
    top3_wide = top3.pivot_table(index=["Protein.Group", "Protein", "Protein.Description"],
                                  columns="run_id", values="Top_3pep_counts_mean", aggfunc="first")
    top3_wide = top3_wide.reindex(columns=run_order(cfg))
    top3_wide.to_csv(ROOT / "outputs" / "tables" / "top3_counts_wide.tsv", sep="\t")
    log2_safe(top3_wide).to_csv(ROOT / "outputs" / "tables" / "top3_log2_wide.tsv", sep="\t")
    print(f"[top3] {top3_wide.shape[0]} proteins × {top3_wide.shape[1]} runs")

    detected_per_run = (~log2_wide.isna()).sum(axis=0)
    detected_per_protein = (~log2_wide.isna()).sum(axis=1)
    print("\nProteins detected per run (min/median/max):",
          int(detected_per_run.min()), int(detected_per_run.median()), int(detected_per_run.max()))
    print("Runs detecting each protein (min/median/max):",
          int(detected_per_protein.min()), int(detected_per_protein.median()), int(detected_per_protein.max()))

    print(f"\n[done] metadata + wide matrices written to {ROOT/'outputs'/'tables'}/")


if __name__ == "__main__":
    main()
