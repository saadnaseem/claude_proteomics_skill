#!/usr/bin/env bash
# Reproducible end-to-end PRT1296 analysis pipeline.
# Run from project root (the directory containing this scripts/ folder).

set -euo pipefail
cd "$(dirname "$0")/.."

# Activate venv if present
if [ -d .venv ]; then
  source .venv/bin/activate
fi

echo "=========================================="
echo "PRT1296 analysis pipeline"
echo "started:  $(date)"
echo "python:   $(python --version)"
echo "=========================================="

for s in scripts/00_load_and_metadata.py \
         scripts/01_qc.py \
         scripts/02_exploratory.py \
         scripts/03_two_factor_anova.py \
         scripts/04_dep_sets.py \
         scripts/05_heatmap.py \
         scripts/06_enrichment.py \
         scripts/07_modules.py \
         scripts/08_interaction_profiles.py \
         scripts/09_extended_analysis.py; do
  echo ""
  echo "--- running $s ---"
  python "$s"
done

echo ""
echo "=========================================="
echo "done.  $(date)"
echo "outputs are in outputs/"
echo "=========================================="
