# PMS2 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to PMS2, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 2626 germline missense variants: 37
  Pathogenic/Likely pathogenic, 50 Benign/Likely benign,
  2539 VUS.
- Structure: AlphaFold model, 862 residues
  (matches UniProt canonical
  length 862).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

826 of 2576
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 122 flagged candidates,
4 are low-confidence.

## Results

**122 of 2539 VUS (4.8%)**
flagged as candidates.

By domain: {'Unannotated / linker': 118, 'Disordered': 4}

## Outputs

- `PMS2_phase1_flagged_candidates.csv`
- `PMS2_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `PMS2_structure_viz.html`
