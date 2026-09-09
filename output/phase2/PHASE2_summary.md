# Phase 2 — cross-gene structural clustering classifier

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as Phase 1: **hypothesis-generation only**. This tests
whether the structural-clustering signal from Phase 1, pooled across
genes, carries any generalizable predictive information — it is not
competing with, or intended to replace, genome-scale tools like
AlphaMissense, and nothing here should be read as "likely pathogenic,"
"reclassified," or "confirmed."

## Project history: 3 genes → 15 genes → 24 genes, then a dual-atom fix

The first pass of Phase 2 (TP53, BRCA1, PTEN; 902 labeled variants) found
real but modest signal, resting on only one statistically trustworthy
gene-held-out fold (BRCA1) and a real confound risk (pLDDT/domain
membership possibly doing much of the work). To address both, 65
ACMG SF v3.2 candidate genes were screened against a 30-Pathogenic-AND-
30-Benign data bar plus structure/numbering validation; 12 more genes
passed, bringing the total to 15 genes. A further screening round added
9 more (ABCA4, ATP7A, COL1A1, COL2A1, COL5A1, MYO7A, SCN1A, SCN2A,
SCN8A), bringing the total to **24 genes and 8,374 labeled variants**.
MUTYH and FBN1 were screened but excluded by structure validation rather
than forced through (MUTYH: 2247 ClinVar/UniProt/structure numbering
mismatches; FBN1: no canonical single-chain AlphaFold model available).

**A further, more significant fix was applied after the 24-gene dataset
was first built.** A one-gene sensitivity check (comparing backbone Cα
distance against side-chain Cβ distance on TP53) surfaced not just an
expected atom-choice difference but a real algorithmic bug: the
nearest-pathogenic-residue search only ever considered the single
globally-closest pathogenic residue. When that residue happened to be a
trivial sequence-adjacent neighbor, it silently masked other valid,
more distant qualifying candidates that also fell within the 6 Å
threshold — undercounting real structural clustering. This was fixed by
restricting the nearest-residue search to only sequence-distance-
qualifying candidates from the start (`_nearest_qualifying` in
`clustering_lib.py`), and Cβ (side-chain) distance was integrated
alongside Cα (backbone) distance throughout the pipeline at the same
time, since both changes touched the same core distance function. A
candidate now qualifies if *either* atom type is within 6 Å, and a
candidate confirmed by *both* atom types (`confirmed_by_both_atoms`) is
treated as stronger evidence. Every number below reflects the fixed,
dual-atom pipeline — Phase 1 flagged-VUS counts increased substantially
across almost every gene as a result (project-wide flagged rate went
from a low single-digit percentage to ~11%).

## Task 1–2: pooled data and class balance

8,374 labeled variants pooled (6,054 Pathogenic, 2,320 Benign — 72%
Pathogenic overall). Per-gene balance varies widely (from LDLR's
739:41 Pathogenic-heavy skew to MSH2's 153:304 Benign-heavy skew), which
is exactly why gene-held-out validation matters more than a pooled
number. Class weighting (`class_weight='balanced'`) applied,
recomputed independently within each training fold. PTEN remains the
one fold where Benign-class metrics aren't statistically meaningful
(n=3) — every other gene clears the 30/30 bar by construction, so this
is no longer a project-wide problem, just a known single-gene
limitation.

## Task 3: models

Logistic regression and a shallow (max_depth=4) random forest.

## Feature design note

"Domain/region" is encoded as a gene-agnostic `disordered_region` flag,
not one-hot of literal per-gene domain names (which would let the model
key off gene identity and undermine gene-held-out validation). pLDDT
(continuous) + `disordered_region` (binary) stands in for domain/region
across genes.

## Task 4: gene-held-out validation — 24 folds

With 24 genes there are 24 held-out folds, and 23 of them (all but
PTEN) have a statistically usable Benign class.

**Headline (balanced accuracy = mean of Benign recall and Pathogenic
recall, so it's not skewed by per-gene class imbalance):**

| Model | Mean balanced accuracy, all 24 folds | Mean balanced accuracy, 23 statistically usable folds |
|---|---|---|
| Logistic Regression | 0.685 | 0.694 |
| Random Forest | 0.703 | 0.705 |

That's the honest number: **real signal, clearly above the 0.50 random
baseline, but well short of anything clinically usable** — essentially
unchanged from the 15-gene estimate (0.72 RF) despite 60% more genes and
the dual-atom/masking-bug fix, which is itself informative: the signal
isn't an artifact of any particular gene subset or of the pre-fix
distance calculation.

**Performance varies a lot by gene** — full per-fold table in
`cv_results_logreg.csv` / `cv_results_rf.csv` (Random Forest shown):
- **Strong folds**: SCN1A (0.847), SCN2A (0.839), TP53 (0.831), HNF1A
  (0.828), SCN8A (0.818).
- **Weak folds**: COL2A1 (0.530 — Pathogenic recall only 0.16, the
  model misses most true Pathogenic COL2A1 variants), COL1A1 (0.565),
  COL5A1 (0.578), MYBPC3 (0.606), RET (0.607). These are real,
  informative failures, not hidden — collagen genes in particular
  appear to be a systematically harder case for this approach.

This spread is itself the finding: the clustering signal generalizes
inconsistently across genes, which is exactly what a hypothesis-
generation tool's limitations section should say plainly rather than
averaging away.

## Task 5: feature importance — Cβ (side-chain) distance is now the top feature

| Feature | LogReg (mean \|coef\|) | Random Forest (mean importance) |
|---|---|---|
| dist3d_cb_A (Cβ distance) | **3.63** (rank 1) | **0.379** (rank 1) |
| dist3d_A (Cα distance) | 2.53 (rank 2) | 0.224 (rank 2) |
| n_pathogenic_within_threshold | 1.04 (rank 3) | 0.205 (rank 3) |
| disordered_region | 0.64 (rank 4) | 0.049 (rank 5) |
| seq_over_3d_ratio | 0.51 (rank 5) | 0.024 (rank 6) |
| seq_minus_3d_diff | 0.26 (rank 6) | 0.009 (rank 9) |
| seqdist_nearest3d | 0.16 (rank 7) | 0.011 (rank 8) |
| seqdist_cb | 0.11 (rank 8) | 0.014 (rank 7) |
| plddt | 0.03 (rank 9) | 0.084 (rank 4) |

Adding side-chain (Cβ) distance did more than fix the masking bug — it
turned out to be the single most informative feature in both models,
ahead of backbone (Cα) distance. This makes physical sense: side chains,
not backbones, are what actually mediate most residue-residue
interactions (binding pockets, hydrophobic packing, hydrogen bonding),
so Cβ-Cβ distance is a more direct proxy for "does this VUS sit in the
same functional neighborhood as a known pathogenic residue" than Cα-Cα
distance alone. pLDDT drops even further down the ranking than in the
15-gene version, reinforcing that the clustering signal is not a
structural-confidence proxy.

## pLDDT ablation — reconfirmed at 24-gene scale

Rerunning the same 24-fold CV with `plddt` removed entirely:

| Model | With pLDDT | Without pLDDT |
|---|---|---|
| Logistic Regression | 0.685 | 0.686 |
| Random Forest | 0.703 | 0.702 |

Removing pLDDT costs essentially nothing (≤0.001 balanced accuracy on
Random Forest) — even more negligible than at 15-gene scale. **The
structural clustering signal is doing genuine, largely independent
predictive work.**

## Task 6: applying the model to VUS, vs. Phase 1's geometric flags

All 32,387 VUS across the 24 genes were scored with a final model
trained on all 8,374 pooled labeled variants (separate from the
held-out-validation models used for the numbers above).

| | Both flag | Neither flags | Geometric only | Model only |
|---|---|---|---|---|
| Count | **3,505** | 22,849 | 144 | 5,889 |

Unlike the 15-gene, pre-fix run (where "model only" cases were 100%
trivial chain-adjacent noise and thus zero non-trivial disagreements
remained), at 24-gene, dual-atom scale **all 5,889 "model only" cases
have genuine sequence separation (>10 residues) from their nearest
pathogenic residue** — the model is now finding non-trivial candidates
Phase 1's strict either-Cα-or-Cβ-within-6Å rule doesn't flag. This is a
real behavioral change introduced by the masking-bug fix and warrants
treating "model only" as a secondary, lower-confidence lead list, not
noise to ignore, though it hasn't been externally validated the way the
"both flag" list has (see `ALPHAMISSENSE_COMPARISON.md`).

**The 3,505 "agree: both flag" cases are the strongest candidate list
this project can currently produce** — flagged independently by both
the geometric rule and a cross-validated statistical model — up from
861 at the pre-fix 15-gene scale. Saved to `vus_scored_phase1_vs_phase2.csv`
(filter to `agreement == 'agree: both flag'`). 59.0% of the subset with
AlphaMissense coverage are independently rated pathogenic-like by
AlphaMissense — see `ALPHAMISSENSE_COMPARISON.md` for the full
external-validation comparison.

The 144 "geometric only" cases are worth a closer look individually —
plausible cases where the model's more conservative probability
disagrees with Phase 1's binary flag.

## Task 7: unbiased test on genes excluded from training

`evaluate_wild_genes.py` tests the final model on 20 genes that were
deliberately excluded from the curated 24-gene training set — not
cherry-picked for good statistics, but rejected specifically because
their ClinVar data failed the 30-Pathogenic/30-Benign bar or structure
validation. Of the 20, 14 had enough post-filtering labeled data to
test (6 were skipped for structure-validation or numbering-mismatch
reasons, a data-quality problem unrelated to this project's model).

| | Random Forest | Logistic Regression |
|---|---|---|
| Mean balanced accuracy, 14 wild genes | **0.686** | 0.684 |
| Mean balanced accuracy, 24 curated genes (for comparison) | 0.703 | 0.685 |

Performance on the unbiased wild-gene set is essentially identical to the
curated 24-gene result (within ~0.02 balanced accuracy either direction)
— real evidence this project's model isn't simply overfit to hand-picked
genes. Full per-gene results in `wild_gene_evaluation.csv`.

## Outputs

- `training_pooled.csv` — 8,374 pooled labeled variants + features
- `all_genes_features.csv` — all 47,999 rows (all buckets, 24 genes) + features
- `cv_results_logreg.csv`, `cv_results_rf.csv` — per-fold precision/recall/balanced accuracy
- `logreg_coefficients_per_fold.csv`, `rf_importances_per_fold.csv` (+ `_no_plddt` ablation versions)
- `vus_scored_phase1_vs_phase2.csv` — all 32,387 VUS, both methods' verdicts
- `wild_gene_evaluation.csv` — unbiased test on 20 excluded genes (Task 7)
- `PHASE2_summary.md` — this document

## Limitations to carry into any future write-up

- 8,374 labeled examples across 24 genes vs. ~71M genome-wide for
  AlphaMissense — still roughly 4 orders of magnitude smaller.
- Performance varies substantially by gene (0.16 to 0.96 Pathogenic
  recall across folds) — there is no single "accuracy" number that
  honestly represents this model; report the distribution, not an
  average, when this matters.
- PTEN remains the one gene where Benign-class validation isn't
  statistically meaningful (n=3).
- MUTYH and FBN1 were excluded for structure-validation reasons, not
  data-quality reasons — worth reconsidering individually if useful
  later (MUTYH's mismatch pattern in particular looks fixable with the
  right transcript).
- Cβ is a reasonable static proxy for side-chain position but still
  can't capture side-chain rotamer flexibility, and (like Cα) is
  computed on a single static AlphaFold model per gene.
- Same caveats as Phase 1 (single static model per gene, BRCA2 excluded
  project-wide for lack of a usable single-chain AlphaFold model).
