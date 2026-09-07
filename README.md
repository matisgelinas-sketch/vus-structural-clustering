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
known pathogenic residue in 3D space *and* more than 10 residues away from
it in the linear sequence — the second condition matters, since without it
almost everything flagged turns out to be trivially adjacent in the chain,
not a meaningful structural signal.

**Phase 2 — cross-gene validation.** The structural-clustering signal is
pooled across genes and used to train a simple, interpretable model
(logistic regression / small random forest) on seven features per variant
(3D distance, sequence distance, their ratio and difference, count of
nearby pathogenic residues, structural confidence, and whether the region
is disordered). Validated with **gene-held-out cross-validation** — the
model is tested only on genes it never saw during training — and benchmarked
against [AlphaMissense](https://github.com/google-deepmind/alphamissense)
as an independent external check.

## Key results

- **24 genes** analyzed, screened from 125 candidates against a 30
  Pathogenic / 30 Benign data-bar plus structure/numbering validation
  (~48,000 classified germline missense variants total)
- **Balanced accuracy ≈ 0.72** across gene-held-out folds — real signal,
  clearly above the 0.50 random baseline, but genuinely uneven gene to
  gene (0.52–0.90)
- Confirmed the clustering signal is **not just a proxy for structural
  confidence** — removing pLDDT as a feature cost almost nothing (0.72 → 0.71)
- Benchmarked against AlphaMissense: the project's highest-confidence
  candidates are **1.46× enriched** for independent AlphaMissense agreement
  over the baseline rate — real external validation, not proof of any
  single candidate
- Also validated on genes deliberately *excluded* from training (rejected
  for having unbalanced ClinVar data) — performance held at roughly the
  same level, evidence this isn't just overfit to hand-picked genes

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

- ~48,000 labeled variants vs. ~71M for genome-scale tools like
  AlphaMissense — this project speaks to whether a structural-clustering
  signal exists and generalizes, not to real-world predictive performance
  at that scale.
- Performance varies substantially by gene; there is no single "accuracy"
  number that honestly represents the whole model.
- Single static AlphaFold structure per gene — no conformational
  flexibility, complexes, or bound states modeled.
- BRCA2 is excluded project-wide: no usable single-chain AlphaFold model
  is available due to protein size.
