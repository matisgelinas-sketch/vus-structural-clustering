# Structural Clustering of Cancer-Risk Gene Variants

A hypothesis-generation tool that uses AlphaFold-predicted protein structures
and a lightweight machine learning model to help prioritize which genetic
**variants of uncertain significance (VUS)** might be worth investigating
further — by testing whether a VUS sits close in 3D space to an
already-known disease-causing residue, even when it's far away in the raw
sequence.

> **Framing, read this first:** this project produces *leads*, not
> *verdicts*. Computational/structural evidence sits at the lowest
> evidentiary tier under ACMG/AMP guidelines and can never, on its own,
> justify a classification change. Nothing in this repository should be
> read as "likely pathogenic," "reclassified," or "confirmed" — the ceiling
> everywhere is **candidate for further investigation**.

## What it does

**Phase 1 — structural flagging.** For each gene, every VUS in ClinVar is
mapped onto its AlphaFold structure and flagged if it sits within 6 Å of a
known pathogenic residue in 3D space — by **either backbone (Cα) or
side-chain (Cβ) distance** — *and* more than 10 residues away from it in
the linear sequence. The sequence-separation condition matters, since
without it almost everything flagged turns out to be trivially adjacent
in the chain, not a meaningful structural signal. Checking both atom
types matters too: a residue's side chain, not its backbone, is often
what's actually close to a known pathogenic residue, and a candidate
confirmed by both atom types is stronger evidence than one found by only
one (`confirmed_by_both_atoms` column). The "nearest qualifying residue"
search only considers candidates that already pass the sequence-distance
filter, so a genuinely close but non-trivial pathogenic residue is never
masked by a closer, trivially-adjacent one.

**Phase 2 — cross-gene validation.** The structural-clustering signal is
pooled across genes and used to train a simple, interpretable model
(logistic regression / small random forest) on nine features per variant
(Cα and Cβ 3D distance, sequence distance, their ratio and difference,
count of nearby pathogenic residues, structural confidence, and whether
the region is disordered). Validated with **gene-held-out cross-validation**
— the model is tested only on genes it never saw during training — and
benchmarked against
[AlphaMissense](https://github.com/google-deepmind/alphamissense) as an
independent external check.

## Key results

- **24 genes** analyzed, screened from 125 candidates against a 30
  Pathogenic / 30 Benign data-bar plus structure/numbering validation
  (~48,000 classified germline missense variants total, 8,374 labeled
  Pathogenic/Benign used for training)
- **Balanced accuracy ≈ 0.70** across gene-held-out folds — real signal,
  clearly above the 0.50 random baseline, but genuinely uneven gene to
  gene. Side-chain (Cβ) distance is the single strongest feature in both
  models, ahead of backbone (Cα) distance.
- Confirmed the clustering signal is **not just a proxy for structural
  confidence** — removing pLDDT as a feature cost almost nothing
  (0.703 → 0.702 balanced accuracy)
- Benchmarked against AlphaMissense: the project's highest-confidence
  candidates are **1.75× enriched** for independent AlphaMissense agreement
  over the baseline VUS rate (59.0% vs. 33.8%) — real external validation,
  not proof of any single candidate
- Also validated on genes deliberately *excluded* from training (rejected
  for having unbalanced ClinVar data) — of 20 such genes, 14 had enough
  data to test, and mean balanced accuracy (0.686) held within ~0.02 of
  the curated 24-gene set — evidence this isn't just overfit to
  hand-picked genes

Full write-ups: [`output/phase2/PHASE2_summary.md`](output/phase2/PHASE2_summary.md),
[`output/phase2/ALPHAMISSENSE_COMPARISON.md`](output/phase2/ALPHAMISSENSE_COMPARISON.md).

An interactive results dashboard and a presentation deck (with a slide-by-slide
speaker script) are in [`presentation/`](presentation/).

## Repository layout

```
scripts/            all pipeline code
  clustering_lib.py       core parsing/distance/structure logic
  domains.py               auto-generated domain annotations (do not hand-edit)
  gene_pipeline.py         add new genes: fetch + validate + save, resumable
  run_new_genes_phase1.py  Phase 1 for every auto-sourced gene
  build_phase2_dataset.py  pools features across all genes
  train_phase2_model.py    gene-held-out CV, feature importance
  score_vus_phase2.py      scores every VUS, compares to Phase 1 flags
  compare_alphamissense.py external benchmark against AlphaMissense
  evaluate_wild_genes.py   tests the model on genes excluded from training
  phase2b/                 NCBI E-utilities ClinVar sourcing + domain generation

Genes/{GENE}/        per-gene input data (sequence always included;
                      structures and auto-fetched ClinVar data are
                      .gitignored — see "Reproducing" below)

output/{GENE}/        per-gene Phase 1 results (flagged candidates, summary)
output/phase2/         pooled Phase 2 results, model outputs, comparisons

presentation/         slide deck + speaker script
```

## Reproducing / regenerating excluded data

This repo excludes large, fully regenerable files (raw AlphaFold structures,
auto-fetched ClinVar exports, bulky intermediate feature tables — see
`.gitignore`) to keep it lean. To rebuild them:

```bash
pip install -r requirements.txt

# optional but ~3x faster: a free NCBI API key
# https://www.ncbi.nlm.nih.gov/account/settings/ ("API Key Management")
export NCBI_API_KEY=your_key_here

# (re-)fetch + validate genes listed in scripts/candidate_genes.json,
# or pass gene symbols directly
python3 scripts/gene_pipeline.py

# Phase 1 for every auto-sourced gene (TP53/BRCA1/PTEN use their own
# run_{gene}_phase1_part1.py / part2.py scripts, since they were sourced
# from manually-exported ClinVar TSVs rather than the API)
python3 scripts/run_new_genes_phase1.py

# Phase 2: pool, train, score, benchmark
python3 scripts/build_phase2_dataset.py
python3 scripts/train_phase2_model.py
python3 scripts/train_phase2_model.py --no-plddt   # ablation
python3 scripts/score_vus_phase2.py
python3 scripts/compare_alphamissense.py
python3 scripts/evaluate_wild_genes.py
```

Every script is resumable/idempotent — safe to interrupt and rerun.

## Limitations

- 8,374 labeled variants vs. ~71M for genome-scale tools like
  AlphaMissense — this project speaks to whether a structural-clustering
  signal exists and generalizes, not to real-world predictive performance
  at that scale.
- Performance varies substantially by gene; there is no single "accuracy"
  number that honestly represents the whole model.
- Single static AlphaFold structure per gene — no conformational
  flexibility, complexes, or bound states modeled. Cβ is a reasonable
  static proxy for side-chain position but still can't capture side-chain
  rotamer flexibility.
- BRCA2 is excluded project-wide: no usable single-chain AlphaFold model
  is available due to protein size.
