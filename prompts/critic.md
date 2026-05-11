# Self-critique template

Apply this after stages 2, 6, and 11. The goal is to catch sloppy reasoning before it propagates.

## Generic checks (apply to all stages)

1. **Specificity**: am I making vague claims ("metabolism is altered") instead of specific ones ("8 of 12 TCA enzymes upregulated, mean log2FC=2.1, p_adj<1e-4")?
2. **Falsifiability**: would I recognize if I were wrong? What evidence would change my mind?
3. **Alternative explanations**: have I considered the simplest competing hypotheses (growth-rate effects, stress response, batch effects, differential extraction efficiency)?
4. **Cherry-picking**: am I ignoring downregulated hits that contradict the up story (or vice versa)?
5. **Numbers integrity**: every count/percentage/p-value I cite — is it directly readable from a table I saved? If not, recompute and verify.

## Stage 2 (hypothesis) — additional checks

- Does the hypothesis name a specific pathway/regulon, or is it generic?
- Is the hypothesis distinguishable from "the cells responded to the perturbation"?
- For this organism, does literature already suggest a different interpretation?
- Does the analysis plan actually test the hypothesis, or just describe the data?

## Stage 6 (enrichment) — additional checks

- Did I use the **correct background** (detected proteins, not the genome)?
- Are enriched terms biologically coherent, or a grab-bag suggesting noise?
- For each "enriched" pathway: how many DE proteins are actually in it? <3 = suspicious despite p-value.
- Do up and down enrichments tell a coherent story (resource reallocation), or are they orthogonal?

## Stage 11 (synthesis report) — additional checks

- Does every paragraph contain a number, a gene name, or a citation? If not, it's filler — cut.
- Are the recommendations actually actionable, or vague ("further studies needed")?
- Have I addressed the falsification criteria from Stage 2?
- Have I named caveats prominently, or buried them?
- Would a skeptical reviewer find the strongest-claim-to-evidence link weak?

## Output format

When invoking critique mid-run, output a markdown cell:
```
## Self-critique (Stage <N>)

✓ Specificity: [pass / issue: ...]
✓ Falsifiability: [pass / issue: ...]
✓ Alternatives: [considered: ...; ruled out because ...]
✓ Cherry-picking: [pass / issue: ...]
✓ Numbers: [verified against tables/...]

Revisions: [none / list]
```

Then apply revisions before continuing.
