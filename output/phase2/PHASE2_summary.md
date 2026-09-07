# Phase 2 — cross-gene structural clustering classifier

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as Phase 1: **hypothesis-generation only**. This tests
whether the structural-clustering signal from Phase 1, pooled across
genes, carries any generalizable predictive information — it is not
competing with, or intended to replace, genome-scale tools like
AlphaMissense, and nothing here should be read as "likely pathogenic,"
"reclassified," or "confirmed."

## Project history: 3 genes → 15 genes

The first pass of Phase 2 (TP53, BRCA1, PTEN; 902 labeled variants) found
real but modest signal, resting on only one statistically trustworthy
gene-held-out fold (BRCA1) and a real confound risk (pLDDT/domain
membership possibly doing much of the work). To address both, 65
ACMG SF v3.2 candidate genes were screened against a 30-Pathogenic-AND-
30-Benign data bar plus structure/numbering validation; **12 more genes
passed** (VHL, MLH1, MSH2, MSH6, PMS2, RB1, RET, TSC2, LDLR, MYBPC3,
HNF1A, ENG), bringing the total to **15 genes and 3903 labeled
variants** — a 4.3x increase. Two genes that passed the data bar were
caught and excluded by structure validation rather than forced through:
MUTYH (2247 ClinVar/UniProt/structure numbering mismatches — smells
like an isoform/transcript issue) and FBN1 (no canonical single-chain
AlphaFold model available, 404). Everything below reflects the
15-gene dataset.

## A bug caught and fixed mid-build, worth knowing about

While building leave-one-out distance features (needed so a Pathogenic
variant doesn't trivially "match itself" at 0 Å), the exclusion of a
row's own residue position was initially applied to every row, VUS and
Benign included. That was wrong: when a VUS sits at the *same residue*
as a separately classified Pathogenic variant (a different amino-acid
substitution — this is real, established ACMG signal, the "PM5"
criterion), excluding that position threw away the single strongest
piece of evidence. Fixed: the self-position exclusion now applies only
to Pathogenic-bucket rows. Everything in this document reflects the
corrected features.

## Task 1–2: pooled data and class balance

3903 labeled variants pooled (2406 Pathogenic, 1497 Benign — 62%
Pathogenic overall). Per-gene balance varies widely (from LDLR's
739:41 Pathogenic-heavy skew to MSH2's 153:304 Benign-heavy skew), which
is exactly why gene-held-out validation matters more than a pooled
number. Class weighting (`class_weight='balanced'`) applied,
recomputed independently within each training fold. PTEN remains the
one fold where Benign-class metrics aren't statistically meaningful
(n=3) — every other gene now clears the 30/30 bar by construction, so
this is no longer a project-wide problem, just a known single-gene
limitation.

## Task 3: models

Logistic regression and a shallow (max_depth=4) random forest.

## Feature design note

"Domain/region" is encoded as a gene-agnostic `disordered_region` flag,
not one-hot of literal per-gene domain names (which would let the model
key off gene identity and undermine gene-held-out validation). pLDDT
(continuous) + `disordered_region` (binary) stands in for domain/region
across genes.

## Task 4: gene-held-out validation — 15 folds, not 1

With 15 genes there are now 15 held-out folds instead of 3, and 14 of
them (all but PTEN) have a statistically usable Benign class. This is
the real payoff of the expansion: an actual distribution of
generalization performance instead of one trustworthy data point.

**Headline (balanced accuracy = mean of Benign recall and Pathogenic
recall, so it's not skewed by per-gene class imbalance):**

| Model | Mean balanced accuracy, all 15 folds | Mean balanced accuracy, 14 statistically usable folds |
|---|---|---|
| Logistic Regression | 0.68 | 0.69 |
| Random Forest | 0.72 | 0.72 |

That's the honest number: **real signal, clearly above the 0.50 random
baseline, but well short of anything clinically usable** — consistent
with (not better or worse than) what the 3-gene BRCA1-only estimate
suggested, but now backed by 14 data points instead of 1.

**Performance varies a lot by gene** — full per-fold table in
`cv_results_logreg.csv` / `cv_results_rf.csv`. Some notable patterns:
- **Strong folds**: TP53 (Path recall 0.96–0.97, Benign recall
  0.72), LDLR (Path recall 0.92–0.96, though Benign n=41 is thin),
  VHL (Path recall 0.97–0.99).
- **Weak folds**: MYBPC3 (Path recall only 0.11–0.14 — the model
  misses almost all true Pathogenic MYBPC3 variants), PMS2 (Path
  recall 0.24–0.32), RB1 (Path recall 0.32–0.41). These are real,
  informative failures, not hidden.
- **TSC2 and ENG** land closest to a genuinely balanced, moderate
  result (~0.66–0.78 both classes) — probably the two best
  single-gene pictures of what this approach can honestly deliver.

This spread is itself the finding: the clustering signal generalizes
inconsistently across genes, which is exactly what a hypothesis-
generation tool's limitations section should say plainly rather than
averaging away.

## Task 5: feature importance — the clustering signal is now dominant

| Feature | LogReg (mean \|coef\|) | Random Forest (mean importance) |
|---|---|---|
| n_pathogenic_within_threshold | **1.08** (rank 1) | 0.23 (rank 2) |
| dist3d_A | 0.62 (rank 2) | **0.37** (rank 1) |
| plddt | 0.32 (rank 3) | 0.15 (rank 3) |
| disordered_region | 0.19 (rank 4) | 0.01 (rank 7) |
| seqdist_nearest3d | 0.05 (rank 6) | 0.10 (rank 4) |
| seq_over_3d_ratio | 0.06 (rank 5) | 0.07 (rank 6) |
| seq_minus_3d_diff | 0.02 (rank 7) | 0.07 (rank 5) |

**This is a meaningfully better result than the 3-gene version.** There,
pLDDT competed for the #1 spot in both models, raising a real confound
concern. At 15 genes, the two clustering-specific features
(`n_pathogenic_within_threshold`, `dist3d_A`) are unambiguously ranked
1st and 2nd in both models, with pLDDT dropping to 3rd. More genes
diluted whatever gene-specific structural-confidence pattern the model
could lean on in the smaller dataset — the clustering signal is the
part that held up.

## pLDDT ablation — reconfirmed at 15-gene scale

Rerunning the same 15-fold CV with `plddt` removed entirely:

| Model | Balanced accuracy, all 15 folds | With pLDDT | Without pLDDT |
|---|---|---|---|
| Logistic Regression | — | 0.683 | 0.679 |
| Random Forest | — | 0.722 | 0.711 |

Same conclusion as the 3-gene version, now on much firmer footing:
removing pLDDT costs ~0.01–0.011 balanced accuracy — negligible. **The
structural clustering signal is doing genuine, largely independent
predictive work, confirmed at 4x the original scale.**

## Task 6: applying the model to VUS, vs. Phase 1's geometric flags

All 23,800 VUS across the 15 genes were scored with a final model
trained on all 3903 pooled labeled variants (separate from the
held-out-validation models used for the numbers above).

| | Both flag | Neither flags | Geometric only | Model only |
|---|---|---|---|---|
| Count | **861** | 17,502 | 2 | 5,435 |

Same unpacking as before: of the 5,435 "model only" cases, **100% are
within 10 residues in sequence of a pathogenic residue** — the trivial
chain-adjacent cases Phase 1's rule was deliberately built to exclude.
Unlike the 3-gene run (which had 178 genuinely non-trivial "model only"
disagreements), **at 15-gene scale there are zero** — every VUS the
model rates highly with real sequence separation from a pathogenic
residue turns out to already be geometrically flagged. The two methods
have converged: with more training data, the model's non-trivial
high-confidence calls now track the geometric rule closely rather than
finding independent signal beyond it.

**The 861 "agree: both flag" cases are the strongest candidate list
this project can currently produce** — flagged independently by both
the geometric rule and a cross-validated statistical model — up from 93
at the 3-gene scale. Saved to `vus_scored_phase1_vs_phase2.csv`
(filter to `agreement == 'agree: both flag'`).

The 2 "geometric only" cases (TP53 Ser94Thr, HNF1A Gly181Glu) are both
in mid-confidence regions (pLDDT 48–76) — plausible cases where the
model's probability is more measured than Phase 1's binary flag.

## Outputs

- `training_pooled.csv` — 3903 pooled labeled variants + features
- `all_genes_features.csv` — all 33,075 rows (all buckets, 15 genes) + features
- `cv_results_logreg.csv`, `cv_results_rf.csv` — per-fold precision/recall/balanced accuracy
- `logreg_coefficients_per_fold.csv`, `rf_importances_per_fold.csv` (+ `_no_plddt` ablation versions)
- `vus_scored_phase1_vs_phase2.csv` — all 23,800 VUS, both methods' verdicts
- `model_only_nontrivial_disagreements.csv` — empty at 15-gene scale (see Task 6)
- `PHASE2_summary.md` — this document

## Limitations to carry into any future write-up

- 3903 labeled examples across 15 genes vs. ~71M genome-wide for
  AlphaMissense — still roughly 4 orders of magnitude smaller.
- Performance varies substantially by gene (0.11 to 0.99 Pathogenic
  recall across folds) — there is no single "accuracy" number that
  honestly represents this model; report the distribution, not an
  average, when this matters.
- PTEN remains the one gene where Benign-class validation isn't
  statistically meaningful (n=3).
- MUTYH and FBN1 were excluded for structure-validation reasons, not
  data-quality reasons — worth reconsidering individually if useful
  later (MUTYH's mismatch pattern in particular looks fixable with the
  right transcript).
- Same caveats as Phase 1 (single static model per gene, no
  substitution-severity feature, BRCA2 excluded project-wide for lack
  of a usable single-chain AlphaFold model).
