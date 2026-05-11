# Microbial proteomics playbook

Domain knowledge for interpreting microbial DIA-NN/Spectronaut/MaxQuant results. Use this as a thinking aid, not a script.

## What you're typically looking at

A microbial proteomics dataset with t-test results between two conditions usually represents one of:
- **Growth-condition contrast**: rich vs minimal media, carbon source A vs B, with vs without N source, aerobic vs anaerobic
- **Stress response**: heat shock, oxidative stress, antibiotic exposure, low pH, high osmolarity
- **Genetic perturbation**: WT vs knockout, overexpression vs control
- **Strain comparison**: WT vs evolved, lab strain vs environmental isolate
- **Time course endpoints**: lag vs exponential vs stationary

The biological story differs hugely. Always ask the user which it is in Stage 0 if the file context doesn't make it clear.

## Reference functional categories (microbial)

| Category | Typical signal in microbial proteomics |
|---|---|
| **TCA cycle** | Coordinated up = oxidative metabolism / respiration. Coordinated down = fermentation switch. Genes: gltA, sucA, sucB, sdhABCD, fumABC, mdh, icd. |
| **Glycolysis / EMP** | Up under glucose growth. Genes: pgi, pfkA, fbaA, gapA, eno, pyk. |
| **Pentose phosphate** | Up when NADPH demand high (biosynthesis, oxidative stress). Genes: zwf, gnd, tktA, talA. |
| **Ribosomal proteins** | Down = growth slowing or stringent response. Up = fast growth. Genes: rpsA-Z, rplA-Z. |
| **Stringent response** | RelA/SpoT, ppGpp synthesis. Hard to detect at protein level — look for downstream: ribosome down + amino acid biosynthesis up. |
| **Stress response (general)** | Sigma factor changes (RpoS, RpoH, RpoE), chaperones (DnaK, GroEL, GroES, IbpA, ClpB), proteases (Lon, ClpP, FtsH). |
| **Oxidative stress** | KatA/KatE catalase, SodA/SodB superoxide dismutase, AhpC/AhpF, Trx, glutathione system. |
| **Osmotic stress** | Trehalose synthesis (otsA, otsB), proline/glycine betaine transport (proU, betA). |
| **Carbon limitation** | High-affinity transporters (PTS reorganization), gluconeogenesis (pckA, ppsA), glyoxylate shunt (aceA, aceB). |
| **Nitrogen limitation** | NtrC regulon, glutamine synthetase (GlnA), GS-GOGAT, ammonium transporters (Amt). |
| **Iron limitation** | Siderophore synthesis (entABCDEF in E. coli, pvdA-Z in P. aeruginosa), Fur regulon — many proteins co-regulated. |
| **Sulfur metabolism** | Cys regulon, CysB, sulfate transporters. |
| **Quorum sensing / biofilm** | LasR, RhlR (P. aeruginosa); LuxR/LuxI; biofilm matrix proteins. |
| **Two-component systems** | Histidine kinase + response regulator pairs. Often the response regulator is more abundant. |

## Patterns that should trigger deeper inquiry

- **Coordinated up/down of an entire pathway**: strong, real signal. The gold standard for confidence.
- **Up of biosynthesis + down of import**: cell is becoming self-sufficient for that nutrient.
- **Up of import + down of biosynthesis**: cell is exploiting external supply.
- **Single high-magnitude hit, no pathway context**: could be a genuine specific regulator, or a noise spike. Check replicate variance.
- **Sigma factor up + its known regulon up**: strong evidence for transcriptional reprogramming.
- **Many hits but no enriched terms**: either dispersed novel biology, or noise. Look at COG/KOG categories instead of GO.

## Microbe-specific notes

### E. coli
- Most extensively characterized. EcoCyc + KEGG very complete.
- Key regulators to know: ArcA/B, FNR, Crp, Lrp, IHF, H-NS, RpoS, RpoH, RpoE.
- Beware: lab strains (K-12, B) accumulate mutations.

### Pseudomonas spp. (aeruginosa, putida, fluorescens)
- Diverse metabolism, strong biotech interest (P. putida).
- Two-component systems heavily used. GacS/GacA, RetS, LadS regulate virulence/biofilm.
- KEGG codes: pae (aeruginosa PAO1), ppu (putida KT2440), pfo (fluorescens).

### Bacillus spp.
- Sporulation regulon (Spo0A) dominates many stress conditions.
- Sigma factor cascade: SigA → SigH → SigF → SigE → SigG → SigK.

### Rhodopseudomonas / phototrophic bacteria
- Photosynthetic apparatus reorganization. PufBALMC, BchXYZ.
- Metabolic flexibility: photoheterotrophy, anaerobic respiration, N2 fixation.

### Cyanobacteria (Synechocystis, Synechococcus)
- Diel cycle effects huge. Always check time-of-day.
- Photosystem stoichiometry shifts under stress.

### Methanogens / archaea
- KEGG annotation often sparser. Use additional DBs (BacDive, IMG/M).
- Watch for unique cofactors: F420, methanopterin, CoM.

## Common pitfalls

1. **Membrane proteins underrepresented** in DIA-NN/MaxQuant — extraction bias. Don't over-interpret absence.
2. **Highly basic / very small proteins** often missed. Same caveat.
3. **PTM-modified peptides** routed away from base peptide — can artifact-deplete the parent protein quantification if PTM is differentially abundant.
4. **Iso-form / paralog conflicts** — `Protein.Group` lumps them; treat as a unit, not as one accession.
5. **Time vs treatment confounding** — if conditions were sampled at different times, growth-phase changes can dominate.
6. **Carry-over from previous batch** — see if any "DE" proteins are abundant ones from a prior sample.

## When the data tells no story

If enrichment is sparse and the network is fragmented:
- Confirm the perturbation actually worked (talk to user about positive controls)
- Check QC: was extraction even? CVs reasonable?
- Try GSEA (rank-based) instead of over-representation
- Look at top-N genes by raw rank, not significance — sometimes biology hides below the threshold
