# External benchmark: this project vs. AlphaMissense

*Draft working notes, not submission-ready.*

## What this is

AlphaMissense (Google DeepMind) is a genome-scale deep-learning
pathogenicity predictor, trained on ~71M variants using evolutionary
conservation and structural context. This project's structural-clustering
model is ~5 orders of magnitude smaller and was never intended to
compete with it — the point of this comparison is to use it as an
independent check: does an established, unrelated method back up (or
diverge from) what this project's simple geometric/statistical approach
finds?

Matching was done via AlphaFold DB's per-protein "aa-substitutions" CSV
(one file per UniProt accession, covering every possible amino-acid
substitution at every canonical-sequence position) — exact match on
wild-type residue + position + mutant residue, no genomic-coordinate or
transcript-mapping ambiguity. 22 of 24 genes had AlphaMissense coverage
(COL1A1 and HNF1A returned 404 — AlphaMissense's own published coverage
is ~92% of the proteome, so 2/24 missing is consistent with a known,
genuine gap in AlphaMissense itself, not a bug in this pipeline).
96% of all variants across the 22 covered genes matched successfully.

## 1. AlphaMissense's own accuracy on this project's exact labeled variants

| | This project (RF, gene-held-out) | AlphaMissense (own classification) |
|---|---|---|
| Mean balanced accuracy | ~0.70 | **0.862** |

AlphaMissense is meaningfully more accurate — expected, and consistent
with the framing maintained throughout this project: a 9-feature
structural-clustering model trained on ~8,000 variants was never going
to outperform a genome-scale deep learning model. This number is the
honest ceiling this project is working under, not a competition it wins.

Per-gene range: 0.70 (MYBPC3) to 0.92 (TP53) — full table in
`alphamissense_per_gene_accuracy.csv`. Notably, MYBPC3 was also this
project's own weakest fold (dropping the LogReg/RF held-out result well
below the mean) — a sign that gene may simply be harder to predict in
general, by any method, not a project-specific weakness.

## 2. Correlation with this project's model, across all VUS

31,459 VUS scored by both methods.

- Pearson r = **0.368** (p ≈ 0)
- Spearman r = **0.376** (p ≈ 0)

A moderate, highly significant positive correlation — real partial
overlap in signal, but far from redundant. This project's model is
capturing *something* related to what AlphaMissense captures (both are
ultimately about structural/functional disruption), but a substantial
amount of what each method sees, the other doesn't.

## 3. The key number: does AlphaMissense back up this project's own top candidates?

Of the **3,426 VUS** flagged by *both* this project's geometric rule
*and* its trained model (the project's highest-confidence output, now
using dual Cα/Cβ distance — see the pipeline-wide fix note below),
AlphaMissense independently rates **2,023 (59.0%)** as pathogenic-like
(score ≥ 0.564, AlphaMissense's own published threshold).

Among *all* VUS in these same genes, regardless of this project's
flagging, 33.8% score AlphaMissense-pathogenic-like. So the flagged
candidate list is enriched **1.75×** over the baseline rate — a real,
quantifiable signal that this project's flagging concentrates
AlphaMissense-supported candidates well above chance, and (unlike the
earlier pre-fix run) a clear majority of the top candidate list is now
independently corroborated.

## 4. A pipeline-wide fix that changed these numbers substantially

Between the previous version of this comparison and this one, two
changes were made to the core distance calculation: (1) side-chain (Cβ)
distance was added alongside backbone (Cα) distance — a candidate now
qualifies if *either* atom type is within 6 Å; (2) a masking bug was
fixed where the pipeline only ever considered the single globally-nearest
pathogenic residue, so a trivial sequence-adjacent neighbor could hide a
more distant, and more informative, qualifying candidate that was also
within 6 Å. Both changes substantially increased the number of flagged
candidates across every gene (Phase 1 flagged-VUS rate roughly tripled
project-wide) and, as a result, both this project's own cross-validated
accuracy and the AlphaMissense enrichment factor above reflect the
corrected pipeline, not the earlier Cα-only, single-nearest-neighbor
version.

## 5. The genuinely divergent cases

1,403 of the 3,426 candidates are ones this project's model rates highly
confident (often >0.9 probability) while AlphaMissense rates as
Likely-benign or Ambiguous with a low score. Two honest interpretations,
not favored over the other without further work:

- **Complementary signal**: pure 3D spatial proximity to a known
  pathogenic residue is a genuinely different kind of evidence than
  evolutionary conservation — it's plausible this project's method
  catches real cases AlphaMissense's conservation-driven approach
  under-weights.
- **False positives**: this project's simpler model, trained on far
  less data, may simply be wrong on these specific cases.

Top divergent examples (full list in
`high_confidence_candidates_vs_alphamissense.csv`):

| Gene | Variant | This project | AlphaMissense |
|---|---|---|---|
| SCN1A | p.Gln1719Glu | 1.00 | 0.38 (Ambiguous) |
| SCN8A | p.Ile1624Val | 1.00 | 0.51 (Ambiguous) |
| MSH2 | p.Ile691Ser | 0.99 | 0.49 (Ambiguous) |
| VHL | p.Arg64Cys | 0.99 | 0.28 (Likely benign) |
| ABCA4 | p.Val1774Ala | 0.99 | 0.50 (Ambiguous) |

These are exactly the kind of cases worth a closer manual look (same
process as the TSC2 lead) — not because either method is "right," but
because genuine disagreement between two independent, methodologically
unrelated predictors is itself informative.

## Outputs

- `all_genes_features_with_alphamissense.csv` — every variant, all
  buckets, with AlphaMissense score/class merged in
- `alphamissense_per_gene_accuracy.csv` — AlphaMissense's own accuracy
  per gene against ClinVar truth
- `high_confidence_candidates_vs_alphamissense.csv` — the 3,426
  candidates with AlphaMissense's verdict alongside
- `ALPHAMISSENSE_COMPARISON.md` — this document

## Limitations

- AlphaMissense's 0.564 pathogenic threshold is applied uniformly here;
  its native 3-class system (Likely benign / Ambiguous / Likely
  pathogenic) has a real "ambiguous" middle band this binary framing
  collapses.
- COL1A1 and HNF1A have no AlphaMissense coverage at all — comparisons
  above exclude them entirely.
- This comparison used AlphaMissense's static, precomputed scores; no
  aspect of this project's model was tuned to agree with it, and
  AlphaMissense had no access to this project's ClinVar labels or
  structural-clustering features.
