# Stage 2 — Hypothesis & analysis plan

## Goal
Generate a focused biological hypothesis and a step-by-step analysis plan **before** writing any biology code. Then self-critique and refine.

## Inputs available to you
- The QC summary from Stage 1
- The microbe + study context the user provided
- The dataset summary (gene/protein names of top 50 hits each direction)
- Reference: `~/.claude/skills/proteomics-agent/reference/microbial_proteomics_playbook.md`
- Reference: `~/.claude/skills/proteomics-agent/reference/interpretation_heuristics.md`

## Step 1 — Pre-look at top hits

In the notebook, display:
- Top 25 proteins by p_adj (most significant), with log2FC and Protein.Names
- Top 15 most upregulated (largest positive log2FC, p_adj<0.05) with names
- Top 15 most downregulated with names

This grounds the hypothesis in actual gene names, not abstractions.

## Step 2 — Generate hypothesis

Write a markdown cell with:
- **Hypothesis**: one sentence. Specific. Names a pathway, regulon, or functional class. Example: *"In condition A, P. putida KT2440 upregulates the TCA cycle and downregulates ribosomal biogenesis, consistent with carbon-source switching from glucose to a slower carbon source."*
- **Why**: 2-3 bullets pointing to specific top hits that motivate the hypothesis.
- **Plan**: numbered steps mapping to stages 3-11. Be specific about *which* GO terms / KEGG pathways / proteins to focus on.
- **Falsification**: what would we expect to see in the enrichment / network / structural stages if the hypothesis is right? What would falsify it?

## Step 3 — Self-critique

Apply `prompts/critic.md` to your hypothesis. Specifically:
- Is the hypothesis falsifiable, or just a description?
- Could the top hits be explained by a simpler alternative (e.g., growth-rate effects, stress response)?
- Are you ignoring downregulated hits that contradict the story?
- For this organism, are there published priors that conflict?

If critique surfaces issues, revise the hypothesis. Iterate up to 2 times.

## ✋ Gate

Output to user:
```
Stage 2 — Proposed hypothesis:

  H: <one-sentence hypothesis>

  Motivated by:
    - <top hit name> (log2FC=X, p_adj=Y)
    - <top hit name> (log2FC=X, p_adj=Y)
    - <top hit name> (log2FC=X, p_adj=Y)

  Plan:
    1. Deep Research: <specific search>
    2. DE call: defaults p_adj<0.05, |log2FC|>1
    3. UniProt annotation for ~250 DE proteins (~5min, cached)
    4. GO + KEGG enrichment, focus on <X, Y, Z categories>
    5. KEGG pathway maps for <pathway_id_1>, <pathway_id_2>
    6. STRING network for top 100 DE proteins
    7. AlphaFold for top 10 by |log2FC|
    8. UniProt PTM annotation
    9. Synthesis report

  Falsification: <what would change the story>

Proceed? (yes / change hypothesis / change plan)
```

If user says "change hypothesis" or "focus on X instead", regenerate and re-gate.
