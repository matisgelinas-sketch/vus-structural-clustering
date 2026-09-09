# COL2A1 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to COL2A1, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1284 germline missense variants: 489
  Pathogenic/Likely pathogenic, 129 Benign/Likely benign,
  666 VUS.
- Structure: AlphaFold model, 1487 residues
  (matches UniProt canonical
  length 1487).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

975 of 1155
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 24 flagged candidates,
2 are low-confidence.

## Results

**24 of 666 VUS (3.6%)**
flagged as candidates.

By domain: {'Fibrillar collagen NC1': 15, 'VWFC': 7, 'Disordered': 2}

## Outputs

- `COL2A1_phase1_flagged_candidates.csv`
- `COL2A1_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `COL2A1_structure_viz.html`
