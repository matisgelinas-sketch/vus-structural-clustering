# LDLR — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to LDLR, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1637 germline missense variants: 739
  Pathogenic/Likely pathogenic, 41 Benign/Likely benign,
  857 VUS.
- Structure: AlphaFold model, 860 residues
  (matches UniProt canonical
  length 860).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

284 of 1596
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 27 flagged candidates,
0 are low-confidence.

## Results

**27 of 857 VUS (3.2%)**
flagged as candidates.

By domain: {'LDL-receptor class B 2': 9, 'LDL-receptor class B 5': 5, 'LDL-receptor class B 6': 4, 'LDL-receptor class B 3': 3, 'LDL-receptor class B 4': 2, 'EGF-like 3': 2, 'LDL-receptor class A 6': 1, 'LDL-receptor class A 3': 1}

## Outputs

- `LDLR_phase1_flagged_candidates.csv`
- `LDLR_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `LDLR_structure_viz.html`
