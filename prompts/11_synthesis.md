# Stage 11 — Synthesis report

## Goal
Write `report.md` — the standalone biological story a wet-lab collaborator can read in 10 minutes and act on. Cross-references the notebook for technical detail.

## Structure (use `templates/report_template.md` as the skeleton)

1. **TL;DR** (3-5 bullets, numbers-first): conditions, organism, n_DE_up/down, the headline finding, the strongest follow-up suggestion.
2. **Methods** (≤1 paragraph): input file, cutoffs, packages + versions (from manifest.json), DBs queried (UniProt vYY.MM, STRING v12.0, KEGG <date>, AlphaFold v4).
3. **QC** (3-4 bullets from Stage 1).
4. **The biological story** (the central section, 4-8 paragraphs):
   - Lead with the hypothesis (Stage 2) and whether enrichment / network / structural data support or refute it.
   - Walk through the strongest enriched pathways with the actual gene names — not just "TCA cycle up" but "*aceA, gltA, sucA, sdhA* are all upregulated 2-4× (p_adj<1e-3)".
   - Pull in 2-4 papers from Stage 3 Deep Research that directly contextualize the finding ("Smith et al. 2024 showed similar TCA upregulation in *P. putida* under acetate growth").
   - Highlight regulatory PTM candidates from Stage 10 if they tie into the story.
   - Discuss surprises: hits that contradict the hypothesis or are novel.
5. **Hub proteins / key regulators** (from Stage 8 + structural notes from Stage 9).
6. **Falsification check**: did the data match the Stage 2 falsification criteria?
7. **Recommended follow-up** (the most important section for the user):
   - Specific wet-lab experiments (knockouts, growth assays, targeted MS for PTMs)
   - Specific computational extensions (e.g., transcriptomic correlation if RNA-seq exists, integration with metabolomics)
   - Specific candidate proteins for biochemical characterization, ranked
8. **Caveats**:
   - Single-replicate proteins flagged in Stage 1
   - Annotation gaps (Trembl-only proteins)
   - Pathway/network coverage limits
   - Any degraded_stages from manifest.json
9. **Appendix**:
   - Link to notebook
   - Link to all tables
   - Link to all figures
   - Citations (DOIs, accession numbers, DB versions)

## Self-critique pass

After drafting, apply `prompts/critic.md` once:
- Does the report actually answer the original biological question?
- Is every claim backed by a number, a figure, or a citation?
- Are alternative explanations addressed?
- Would a wet-lab collaborator know what to do tomorrow morning?

Revise as needed.

## ✋ Final gate

```
Stage 11 — Report draft complete:
  Length: ~1,800 words
  Citations: 7 papers, 12 DBs, 247 proteins
  Top recommendation: <one-sentence>
  Saved: report.md

Review and tell me what to revise, or accept as final?
```

If user requests revisions, edit and re-gate. When accepted, write a final entry to `manifest.json` with `final_report_at: <iso8601>`.
