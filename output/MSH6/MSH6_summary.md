# MSH6 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to MSH6, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 4142 germline missense variants: 73
  Pathogenic/Likely pathogenic, 81 Benign/Likely benign,
  3988 VUS.
- Structure: AlphaFold model, 1360 residues
  (matches UniProt canonical
  length 1360).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

1072 of 4061
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 129 flagged candidates,
0 are low-confidence.

## Results

**129 of 3988 VUS (3.2%)**
flagged as candidates.

By domain: {'Unannotated / linker': 129}

## Outputs

- `MSH6_phase1_flagged_candidates.csv`
- `MSH6_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `MSH6_structure_viz.html`
