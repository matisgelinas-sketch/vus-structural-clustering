# MYBPC3 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to MYBPC3, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1689 germline missense variants: 37
  Pathogenic/Likely pathogenic, 35 Benign/Likely benign,
  1617 VUS.
- Structure: AlphaFold model, 1274 residues
  (matches UniProt canonical
  length 1274).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

306 of 1654
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 127 flagged candidates,
2 are low-confidence.

## Results

**127 of 1617 VUS (7.9%)**
flagged as candidates.

By domain: {'Fibronectin type-III 1': 35, 'Ig-like C2-type 1': 24, 'Ig-like C2-type 7': 22, 'Ig-like C2-type 3': 12, 'Ig-like C2-type 2': 11, 'Unannotated / linker': 11, 'Ig-like C2-type 5': 10, 'Ig-like C2-type 4': 2}

## Outputs

- `MYBPC3_phase1_flagged_candidates.csv`
- `MYBPC3_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `MYBPC3_structure_viz.html`
