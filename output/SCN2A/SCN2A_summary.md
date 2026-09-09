# SCN2A — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to SCN2A, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1591 germline missense variants: 420
  Pathogenic/Likely pathogenic, 31 Benign/Likely benign,
  1140 VUS.
- Structure: AlphaFold model, 2005 residues
  (matches UniProt canonical
  length 2005).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

584 of 1560
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 96 flagged candidates,
8 are low-confidence.

## Results

**96 of 1140 VUS (8.4%)**
flagged as candidates.

By domain: {'I': 37, 'II': 22, 'IV': 18, 'III': 13, 'Unannotated / linker': 6}

## Outputs

- `SCN2A_phase1_flagged_candidates.csv`
- `SCN2A_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `SCN2A_structure_viz.html`
