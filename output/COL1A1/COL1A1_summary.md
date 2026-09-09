# COL1A1 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to COL1A1, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1178 germline missense variants: 500
  Pathogenic/Likely pathogenic, 82 Benign/Likely benign,
  596 VUS.
- Structure: AlphaFold model, 1464 residues
  (matches UniProt canonical
  length 1464).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

894 of 1096
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 54 flagged candidates,
11 are low-confidence.

## Results

**54 of 596 VUS (9.1%)**
flagged as candidates.

By domain: {'Fibrillar collagen NC1': 42, 'Disordered': 11, 'VWFC': 1}

## Outputs

- `COL1A1_phase1_flagged_candidates.csv`
- `COL1A1_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `COL1A1_structure_viz.html`
