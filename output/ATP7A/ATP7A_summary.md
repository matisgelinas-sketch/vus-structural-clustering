# ATP7A — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to ATP7A, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 914 germline missense variants: 42
  Pathogenic/Likely pathogenic, 202 Benign/Likely benign,
  670 VUS.
- Structure: AlphaFold model, 1500 residues
  (matches UniProt canonical
  length 1500).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å in 3D AND > 10 residues apart in
  sequence** (same dual criterion validated on TP53, to exclude trivial
  chain-adjacent cases).

## Structural confidence

206 of 712
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 17 flagged candidates,
0 are low-confidence.

## Results

**17 of 670 VUS (2.5%)**
flagged as candidates.

By domain: {'Unannotated / linker': 17}

## Outputs

- `ATP7A_phase1_flagged_candidates.csv`
- `ATP7A_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `ATP7A_structure_viz.html`
