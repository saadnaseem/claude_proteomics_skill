"""Shared helpers for the PRT1296 analysis pipeline."""
from __future__ import annotations

from pathlib import Path
import yaml
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "analysis.yaml"
MODULES_PATH = ROOT / "config" / "tolerance_modules.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def load_modules() -> dict:
    with open(MODULES_PATH) as f:
        return yaml.safe_load(f)


def metadata_df(cfg: dict) -> pd.DataFrame:
    rows = []
    for sid, info in cfg["samples"].items():
        for rep in cfg["replicates"]:
            rows.append({
                "run_id":   f"{sid}-{rep}",
                "sample":   sid,
                "replicate": rep,
                "strain":   info["strain"],
                "medium":   info["medium"],
                "label":    info["label"],
                "group":    info["label"],
            })
    return pd.DataFrame(rows)


def load_protein_long(cfg: dict) -> pd.DataFrame:
    """Long-format protein counts: one row per (protein, sample, replicate)."""
    df = pd.read_csv(ROOT / "data" / "proteins_long.csv")
    return df


def load_top3(cfg: dict) -> pd.DataFrame:
    df = pd.read_csv(ROOT / "data" / "top3_proteins.csv")
    return df


def long_to_wide(df: pd.DataFrame, value_col: str = "Counts_sum") -> pd.DataFrame:
    """Pivot long format to protein × run wide matrix."""
    df = df.copy()
    df["run_id"] = df["Sample"] + "-" + df["Replicate"]
    wide = df.pivot_table(index=["Protein.Group", "Protein", "Protein.Description"],
                          columns="run_id", values=value_col, aggfunc="first")
    return wide


def log2_safe(x):
    """Log2 transform, treating <=0 and NaN as missing. Preserves pandas index/columns."""
    if isinstance(x, pd.DataFrame):
        arr = np.where((x.values > 0) & np.isfinite(x.values), x.values, np.nan)
        return pd.DataFrame(np.log2(arr, where=~np.isnan(arr), out=np.full_like(arr, np.nan)),
                            index=x.index, columns=x.columns)
    arr = np.where((x > 0) & np.isfinite(x), x, np.nan)
    return np.log2(arr, where=~np.isnan(arr), out=np.full_like(arr, np.nan, dtype=float))


def ensure_dirs(cfg: dict) -> None:
    for k in ("figures", "tables", "enrichment", "memo"):
        (ROOT / cfg["paths"][k]).mkdir(parents=True, exist_ok=True)


def palette(cfg: dict, kind: str) -> dict:
    return cfg["plot"][f"palette_{kind}"]


def group_order(cfg: dict) -> list[str]:
    return [info["label"] for info in cfg["samples"].values()]


def run_order(cfg: dict) -> list[str]:
    return [f"{sid}-{rep}" for sid in cfg["samples"] for rep in cfg["replicates"]]
