# Stage 0 — Intake & schema validation

## Goal
Load the xlsx, verify the column structure matches expectations, infer condition labels A/B, resolve the microbe to an NCBI taxonomy ID, and write an initial `manifest.json`.

## Expected schema
Required columns (case-sensitive):
- `Protein` — primary accession (UniProt format, may be semicolon-delimited for protein groups)
- `Protein.Group` — full group string (semicolon-delimited accessions)
- `Protein.Names` — gene/protein name (may be empty)
- `Protein.Description` — free text annotation (may be empty)
- `log2_mean_<conditionA_id>` — mean log2 intensity in condition A
- `log2_mean_<conditionB_id>` — mean log2 intensity in condition B
- `log2_std_<conditionA_id>` — std dev in condition A
- `log2_std_<conditionB_id>` — std dev in condition B
- `t-test_stat` — t-statistic
- `p-value` — raw p-value
- `p_adjusted(BH)` — Benjamini-Hochberg adjusted p
- `log2_Fold_change_A/B` — log2(mean_A / mean_B); positive = up in A

## What to do

1. Create the output dir: `proteomics_run_$(date +%Y%m%d_%H%M%S)/` in cwd, with `figures/`, `tables/`, `annotations/` subdirs.
2. Read the xlsx with `pd.read_excel(path, engine='openpyxl')`. Report number of rows, number of columns, sheet name(s).
3. Auto-detect condition IDs by regex on column names: `^log2_mean_(.+)$` → capture group is the condition ID. Verify there are exactly two and they have matching `log2_std_*` columns.
4. Validate the fold-change column name encodes the same A/B order as the user-provided labels. If mismatched, **stop and ask** rather than silently flipping signs.
5. Resolve organism → NCBI taxid:
   - Use UniProt's taxonomy REST endpoint: `https://rest.uniprot.org/taxonomy/search?query=<name>&format=json&size=5`
   - If multiple candidates, list them and ask user to pick.
   - Cache the result to `annotations/taxonomy_cache.json`.
6. Write `manifest.json` with:
   ```json
   {
     "run_id": "<timestamp>",
     "input_files": [{"path": "...", "sha256": "...", "n_proteins": ...}],
     "conditions": {"A": "SN_0725_92", "B": "SN_0725_83"},
     "organism": {"name": "...", "taxid": "...", "uniprot_proteome_id": "..."},
     "skill_version": "0.1.0",
     "created_at": "<iso8601>"
   }
   ```
7. Build the notebook setup cell from `~/.claude/skills/proteomics-agent/templates/notebook_setup.py`. Use `NotebookEdit` to create `notebook.ipynb` with the setup cell.

## ✋ Gate
After stage 0, output to the user:
```
Stage 0 complete:
  File: <path>
  Proteins: <n>
  Condition A: <label> (<n_replicates> via std presence)
  Condition B: <label>
  Organism: <name> (taxid <id>, UniProt proteome <id>)
  Output dir: <abspath>
Proceed to QC? (yes / change anything?)
```

## Failure modes
- **Missing required column**: list which, and ask the user to either rename in the xlsx or confirm we should map a different column.
- **Multiple sheets**: list them, ask which to use.
- **Organism ambiguous**: show top 5 taxonomy hits, ask which is correct.
- **Schema partially matches** (e.g. uses `log2FC` instead of `log2_Fold_change_A/B`): propose a column mapping and ask user to confirm before proceeding.
