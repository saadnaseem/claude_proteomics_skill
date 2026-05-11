# Stage 3 — Deep Research

## Goal
Build a literature-grounded biological context for the microbe + perturbation. Outputs feed into stages 6 (enrichment interpretation) and 11 (synthesis).

## Tools
Use `WebSearch` for broad queries, `WebFetch` for specific URLs (UniProt proteome page, KEGG organism page, recent papers).

## Searches to run (parallelize with multiple WebSearch calls)

1. **Organism baseline biology**:
   - `<organism name> proteome characterization review`
   - `<organism name> transcriptome <perturbation> proteome`
   - Pull the organism's KEGG entry: `https://www.kegg.jp/entry/<kegg_org_code>` if known
2. **Perturbation literature**:
   - `<organism name> <perturbation, e.g., nitrogen limitation> proteomics differential expression`
   - `<organism name> <perturbation> regulon`
3. **Recent papers (last 3 years)**:
   - Filter results to 2023–2026
   - Prioritize Nature/eLife/mSystems/PLOS Biol/J Proteome Res
4. **Known regulators**:
   - For each top-5 hit identified in Stage 2 (by gene name): `<gene name> <organism> regulation function`
5. **kgMicrobe** (if organism is in their coverage):
   - REST query to `https://kg-hub.berkeleybop.io/kg-microbe/` or use the API at `https://api.monarchinitiative.org/api/`
   - Look for traits / metabolism / growth conditions of the organism

## What to capture

For each useful source, write to a notebook markdown cell:
- 2-3 sentence summary in your own words (do not copy text >15 words verbatim)
- The URL
- One reason this is relevant to the current hypothesis

End the stage with a "Deep Research summary" cell (~10 bullets) that distills what's known and what's novel about the current dataset relative to literature.

Save the raw research notes to `tables/03_deepresearch_notes.md` for the synthesis report to consume.

## Failure handling
If WebSearch returns nothing useful for the organism (rare microbe), fall back to closest characterized relative (e.g., for *Rhodopseudomonas palustris* → check *Rhodobacter sphaeroides* literature) and flag the substitution explicitly.

## No ✋ gate — proceed automatically to Stage 4.
