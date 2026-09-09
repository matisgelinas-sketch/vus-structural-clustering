# TSC2 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to TSC2, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 4160 germline missense variants: 186
  Pathogenic/Likely pathogenic, 198 Benign/Likely benign,
  3776 VUS.
- Structure: AlphaFold model, 1807 residues
  (matches UniProt canonical
  length 1807).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

1375 of 3962
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 350 flagged candidates,
8 are low-confidence.

## Results

**350 of 3776 VUS (9.3%)**
flagged as candidates.

By domain: {'Rap-GAP': 168, 'Unannotated / linker': 145, 'Required for interaction with TSC1': 33, 'Disordered': 4}

## Outputs

- `TSC2_phase1_flagged_candidates.csv`
- `TSC2_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `TSC2_structure_viz.html`
