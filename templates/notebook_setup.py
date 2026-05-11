"""
Notebook setup cell — paste as the first code cell of every proteomics-agent notebook.
Run-specific values (paths, conditions, organism) are filled in by the agent at Stage 0.
"""
import os
import sys
import json
import hashlib
import random
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import requests_cache
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---- Reproducibility ----
SEED = 42
np.random.seed(SEED)
random.seed(SEED)

# ---- Run config (filled at Stage 0) ----
RUN_DIR = Path(os.path.abspath("."))               # set to proteomics_run_<ts>/ at Stage 0
INPUT_FILES = []                                    # list of dicts: {path, sha256, n_proteins}
CONDITIONS = {"A": None, "B": None}                 # e.g. {"A": "SN_0725_92", "B": "SN_0725_83"}
ORGANISM = {"name": None, "taxid": None, "kegg_code": None, "uniprot_proteome_id": None}

# ---- Output dirs ----
FIG_DIR = RUN_DIR / "figures"
TBL_DIR = RUN_DIR / "tables"
ANN_DIR = RUN_DIR / "annotations"
for d in (FIG_DIR, TBL_DIR, ANN_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ---- API cache (mandatory) ----
requests_cache.install_cache(
    str(ANN_DIR / "api_cache"),
    backend="sqlite",
    expire_after=86400 * 30,
    allowable_methods=["GET", "POST"],
)

# ---- Plot defaults ----
plt.rcParams.update({
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})
COLOR_UP = "#D62728"
COLOR_DOWN = "#1F77B4"
COLOR_NS = "#AAAAAA"

# ---- Helpers ----
def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()

def manifest_path() -> Path:
    return RUN_DIR / "manifest.json"

def update_manifest(updates: dict) -> None:
    p = manifest_path()
    data = {}
    if p.exists():
        data = json.loads(p.read_text())
    data.update(updates)
    data["updated_at"] = datetime.utcnow().isoformat() + "Z"
    p.write_text(json.dumps(data, indent=2, default=str))

def degraded(stage: str, reason: str) -> None:
    p = manifest_path()
    data = json.loads(p.read_text()) if p.exists() else {}
    data.setdefault("degraded_stages", []).append(
        {"stage": stage, "reason": reason, "at": datetime.utcnow().isoformat() + "Z"}
    )
    p.write_text(json.dumps(data, indent=2, default=str))

print(f"Setup OK. Run dir: {RUN_DIR}")
print(f"Cache: {ANN_DIR}/api_cache.sqlite")
