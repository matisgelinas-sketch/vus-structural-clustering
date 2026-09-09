# ENG — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to ENG, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 624 germline missense variants: 94
  Pathogenic/Likely pathogenic, 110 Benign/Likely benign,
  420 VUS.
- Structure: AlphaFold model, 658 residues
  (matches UniProt canonical
  length 658).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

101 of 514
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 92 flagged candidates,
2 are low-confidence.

## Results

**92 of 420 VUS (21.9%)**
flagged as candidates.

By domain: {'Required for interaction with GDF2': 60, 'ZP': 24, 'Unannotated / linker': 8}

## Outputs

- `ENG_phase1_flagged_candidates.csv`
- `ENG_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `ENG_structure_viz.html`
