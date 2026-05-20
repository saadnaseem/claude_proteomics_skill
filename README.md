# proteomics-agent

A Claude Code skill that runs an end-to-end microbial proteomics analysis from a pre-computed DIA-NN/Spectronaut/MaxQuant t-test xlsx — differential expression, UniProt annotation, GO/KEGG enrichment, KEGG pathway maps, STRING interaction networks, AlphaFold structural inspection, UniProt PTM annotation, and a literature-grounded synthesis report.

The skill ships as a folder of instructions for [Claude Code](https://claude.com/claude-code) to execute. There is no separate process to launch and **no API keys to manage** — Claude Code is the runtime.

## What it does

Eleven stages, with human-in-the-loop approval gates at the critical points:

| # | Stage | Pause for approval? |
|---|---|---|
| 0 | Intake & schema validation, organism → NCBI taxid | ✋ |
| 1 | QC summary (coverage, missingness, p-value distribution, volcano preview) | — |
| 2 | Hypothesis + analysis plan, with self-critique | ✋ |
| 3 | Deep Research (organism + perturbation literature) | — |
| 4 | Differential expression call (default: p_adj<0.05, \|log2FC\|>1) | ✋ |
| 5 | UniProt annotation for DE set (cached) | — |
| 6 | GO / KEGG over-representation analysis | — |
| 7 | KEGG pathway maps with up/down overlay | — |
| 8 | STRING interaction network + hub & community analysis | ✋ |
| 9 | AlphaFold structural inspection for top hits | ✋ |
| 10 | UniProt PTM site annotation | — |
| 11 | Synthesis report with literature cross-refs | ✋ |

## Expected input

An xlsx with these columns (case-sensitive):

```
Protein, Protein.Group, Protein.Names, Protein.Description,
log2_mean_<conditionA>, log2_mean_<conditionB>,
log2_std_<conditionA>,  log2_std_<conditionB>,
t-test_stat, p-value, p_adjusted(BH), log2_Fold_change_A/B
```

This is downstream-only — raw mzML/RAW files and PTM-site-localized data (e.g., FragPipe `Phospho (STY)Sites.txt`) are out of scope for v0.1.

## Outputs

Everything lands in your current working directory under a timestamped folder:

```
proteomics_run_YYYYMMDD_HHMMSS/
├── notebook.ipynb            # the live notebook
├── report.md                 # standalone biological synthesis
├── figures/                  # volcano, MA, pathway maps, network HTML, AF structures
├── tables/                   # de_proteins.csv, enriched_pathways.csv, hubs.csv, ...
├── annotations/              # cached UniProt/STRING/KEGG/AlphaFold responses
└── manifest.json             # inputs, params, package versions, random seed, degraded stages
```

External API calls are cached to disk via `requests-cache` — re-runs are free and offline-capable.

## Prerequisites

- macOS or Linux (Windows untested; should work via WSL2)
- [Claude Code](https://claude.com/claude-code) installed and authenticated (`claude --version`)
- [Miniconda](https://docs.conda.io/projects/miniconda/) or Anaconda (Python 3.11 is pinned in `env/environment.yml`)
- git
- ~2 GB free disk for the conda environment, plus space for your run outputs and API response cache

You do **not** need any API keys. UniProt, KEGG, STRING, and AlphaFold are all free public APIs.

## Installation

```bash
# 1. Clone into your Claude Code skills directory
mkdir -p ~/.claude/skills
git clone https://github.com/saadnaseem/claude_proteomics_skill.git \
          ~/.claude/skills/proteomics-agent

# 2. Create the conda environment (takes ~5 min)
conda env create -f ~/.claude/skills/proteomics-agent/env/environment.yml

# 3. Register the Jupyter kernel (so .ipynb cells use the right env)
conda run -n proteomics-agent python -m ipykernel install --user \
  --name proteomics-agent --display-name "Python (proteomics-agent)"

# 4. Verify the install
conda run -n proteomics-agent python -c "import pandas, scipy, statsmodels, \
    sklearn, matplotlib, seaborn, bioservices, gseapy, goatools, networkx; \
    print('OK')"
```

If step 4 prints `OK`, the skill is ready. It is now available in every Claude Code session on this machine — no per-project setup needed.

### Updating later

```bash
cd ~/.claude/skills/proteomics-agent
git pull
conda env update -f env/environment.yml --prune
```

### Uninstalling

```bash
conda env remove -n proteomics-agent
rm -rf ~/.claude/skills/proteomics-agent
```

## Usage

### Option 1 — interactive agent (recommended for new analyses)

In any Claude Code session:

```bash
cd ~/wherever-your-data-lives
claude
```

Then either:

```
/proteomics-agent path/to/data.xlsx
```

or natural language:

> *Run a proteomics analysis on `~/data/exp042.xlsx` — conditions A=SN_0725_92 vs B=SN_0725_83, organism Pseudomonas putida KT2440.*

Claude Code will pick up the skill and walk through stages 0–11, pausing at the ✋ gates for your approval. Expect 1–3 hours for a first run (most of which is API calls — subsequent runs are much faster thanks to caching).

### Option 2 — non-interactive scripted pipeline

If you have many similar files (batch mode) or want a reproducible analysis without the interactive interpretation layers, run the example scripts directly. The cleanest starting point is the modular 2-factor reference in [`examples/PRT1296_end_to_end/`](examples/PRT1296_end_to_end/) — see its README for how to adapt to your data. Older simpler examples (`run_pipeline.py`, `phase3_meta_analysis.py`, `build_meta_notebook.py`) are also in `examples/`.

### Quick smoke test

To verify the agent triggers correctly without running real data:

```bash
cd /tmp && claude
```

Then type: `/proteomics-agent`. Claude should respond by asking for the xlsx path, condition labels, and organism — that's the stage-0 intake check.

## Repo layout

```
proteomics-agent/
├── SKILL.md                          # entry point Claude reads on invocation
├── env/environment.yml               # conda spec
├── prompts/
│   ├── 00_intake.md … 11_synthesis.md  # one prompt per pipeline stage
│   ├── coding_guidelines.md            # reproducibility, caching, error handling
│   └── critic.md                       # self-critique template
├── reference/
│   ├── api_endpoints.md              # UniProt, STRING, KEGG, AlphaFold, kgMicrobe URLs
│   ├── microbial_proteomics_playbook.md  # domain knowledge / patterns to look for
│   └── interpretation_heuristics.md  # confidence calibration
├── templates/
│   ├── notebook_setup.py             # standard imports + cache + plotting defaults
│   └── report_template.md            # synthesis report skeleton
├── examples/
│   ├── README.md                     # guide to which example to use when
│   ├── PRT1296_end_to_end/           # full modular 2×3 factorial pipeline (recommended starting point)
│   ├── run_pipeline.py               # single-comparison batch script
│   ├── phase3_meta_analysis.py       # multi-comparison cross-analysis
│   └── build_meta_notebook.py        # emit a populated notebook for handoff
└── README.md
```

## Troubleshooting

| Symptom | Fix |
|---|---|
| `conda: command not found` | Install Miniconda first: <https://docs.conda.io/projects/miniconda/> |
| `claude: command not found` | Install Claude Code: <https://claude.com/claude-code> |
| `ModuleNotFoundError: bioservices` (or any other) inside a run | The active env is wrong. Either `conda activate proteomics-agent` first, or prefix commands with `conda run -n proteomics-agent`. |
| Agent does not trigger on `/proteomics-agent` | Confirm the skill is at `~/.claude/skills/proteomics-agent/SKILL.md` and that you started `claude` (not a separate process). Run `/skills` inside Claude Code to list registered skills. |
| KEGG REST returns 403 / connection refused | KEGG enforces ~3 req/s. The example pipeline already throttles; if you wrote custom code, add `time.sleep(0.3)` between calls. Re-run with cache enabled. |
| `HTTP 400` from UniProt batch lookup with `Q877Q0;Q88FK3`-style accessions | DIA-NN joins protein-group members with `;`. Split on `;` and use the first accession (see `examples/PRT1296_end_to_end/scripts/06_enrichment.py:clean_accession` for the canonical fix). |
| `upsetplot` crashes with `Invalid RGBA argument: nan` on matplotlib ≥ 3.10 | Known version mismatch. The PRT1296 example uses a self-contained UpSet plotter (`scripts/04_dep_sets.py`) that avoids the dep entirely. |
| Re-runs feel as slow as the first run | The API cache lives in your *output* directory, not the skill repo. To share a cache across runs, set `--cache-dir` or symlink `annotations/api_cache.sqlite` between run directories. |

## External resources used at runtime

- **UniProt REST** — protein annotation, taxonomy resolution, PTM sites
- **STRING DB** — protein-protein interactions
- **KEGG REST** — organism pathways, gene mapping, pathway diagrams
- **AlphaFold DB** — predicted structures + pLDDT
- **kgMicrobe / Monarch** — microbe-specific traits and metabolism (with fallback if unavailable)
- **PubMed / WebSearch** — literature for the Deep Research stage

## Acknowledgments

Architecture inspired by [CellVoyager](https://github.com/zou-group/CellVoyager) (Zou Lab) — the plan-and-execute pattern with structured outputs, self-critique, and human-in-the-loop approval gates is adapted from their single-cell RNA-seq agent. This implementation differs by being a Claude Code skill rather than a standalone Python application — no API keys, no MCP server, no Streamlit GUI.

## License

MIT — see [LICENSE](LICENSE).
