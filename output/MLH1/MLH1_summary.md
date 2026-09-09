# MLH1 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to MLH1, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 2011 germline missense variants: 200
  Pathogenic/Likely pathogenic, 37 Benign/Likely benign,
  1774 VUS.
- Structure: AlphaFold model, 756 residues
  (matches UniProt canonical
  length 756).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

468 of 1974
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 391 flagged candidates,
8 are low-confidence.

## Results

**391 of 1774 VUS (22.0%)**
flagged as candidates.

By domain: {'Unannotated / linker': 301, 'Interaction with EXO1': 88, 'Disordered': 2}

## Outputs

- `MLH1_phase1_flagged_candidates.csv`
- `MLH1_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `MLH1_structure_viz.html`
