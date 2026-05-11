# Interpretation heuristics

Patterns to recognize and what they typically mean. Use as evidence weighting, not as conclusions.

## "TCA cycle uniformly upregulated"
- **Likely**: switch to a more oxidative metabolism. Common when carbon source becomes limiting or shifts from a fermentable to a non-fermentable substrate (glucose → acetate, glucose → glycerol).
- **Check**: glyoxylate shunt (aceA/aceB) — if also up, strong evidence for C2 substrate utilization.
- **Falsifier**: if respiratory chain (cyo, cyd, nuo, sdh) is NOT also up, consider fermentation-to-respiration is incomplete.

## "Ribosomal proteins uniformly downregulated"
- **Likely**: growth slowing. Stringent response (ppGpp), nutrient limitation, or stationary phase entry.
- **Check**: amino acid biosynthesis enzymes — if up, stringent response. If also down, general growth arrest.
- **Falsifier**: ribosomal proteins are abundant; even small fold changes are detectable, so not all "ribosome down" signals are biologically meaningful. Require ≥2× and ≥1/3 of subunits.

## "Sigma factor (e.g., RpoS, RpoH, SigB) up + downstream regulon up"
- **Strong evidence** for transcriptional reprogramming via that sigma factor.
- **Cite the regulon source** in report (RegulonDB for E. coli; similar resources per organism).

## "Iron uptake machinery (siderophores, transporters) up"
- **Likely**: iron limitation. Even when not the experimental variable.
- **Check**: Fur regulon coordinately changing? Then iron is genuinely limiting.
- **Implication**: may be a confound — iron limitation triggers many secondary effects.

## "Many membrane / transport proteins changing, no clear pathway"
- **Be cautious**: extraction efficiency for membrane proteins varies. Could be technical.
- **Check**: are they all from the same membrane fraction?

## "Heat shock proteins (DnaK, GroEL, ClpB, IbpA) up"
- **Likely**: protein folding stress. Could be heat, oxidative stress, ribosome stalling, or accumulation of misfolded proteins.
- **Check**: protease genes (Lon, ClpP, FtsH) co-regulated? Then proteostasis stress.

## "Cell envelope / wall biosynthesis (mur*, dac*) up"
- **Likely**: cells dividing or remodeling envelope (e.g., antibiotic exposure, osmotic challenge, biofilm transition).

## "Single very strong hit, no pathway support"
- **Treat as candidate**, not conclusion. Could be:
  - Real specific regulator
  - Off-target effect of perturbation
  - Technical artifact (peptide-level issue inflating one protein)
- **Recommend** orthogonal validation (Western blot, targeted MS, mutant phenotype).

## "Up and down enrichments orthogonal (different processes)"
- Common in real biology — resource reallocation rather than one-directional response.
- **Frame** as: "X is upregulated to enable Y, Z is downregulated to free resources."

## "No enrichment despite many DE hits"
- Possible explanations:
  - Background isn't being applied correctly
  - Annotation is sparse for this organism
  - Real biology is novel and not in annotated pathways
- **Action**: switch to gene-set lookup by COG category, or rank-based GSEA.

## Bacterial regulatory PTMs to know

- **Phospho-Asp** (response regulators, two-component systems): hard to detect by standard MS — but if you see phospho-S/T/Y on a known response regulator, suspect cross-talk.
- **Lysine acetylation on metabolic enzymes** (E. coli AceCS, GltA, GapA): regulatory, often coupled to acetyl-CoA pool.
- **Cysteine S-thiolation** under oxidative stress: glutathionylation, sulfenylation. Rarely captured in standard runs.

## Confidence calibration

When writing the report:

| Evidence level | Language |
|---|---|
| Single hit | "candidate" |
| Pathway-level coordinated change (≥3 genes, p_adj<0.05 each) | "consistent with" |
| Pathway + literature precedent | "supports the model that" |
| Pathway + structural/PTM evidence + literature | "demonstrates" / "establishes" |

Default to weaker language. Reviewers respect calibrated claims; overclaiming is a one-shot reputational cost.
