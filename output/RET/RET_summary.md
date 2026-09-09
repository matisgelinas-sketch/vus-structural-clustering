# RET — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to RET, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1989 germline missense variants: 81
  Pathogenic/Likely pathogenic, 47 Benign/Likely benign,
  1861 VUS.
- Structure: AlphaFold model, 1114 residues
  (matches UniProt canonical
  length 1114).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

519 of 1942
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 95 flagged candidates,
0 are low-confidence.

## Results

**95 of 1861 VUS (5.1%)**
flagged as candidates.

By domain: {'Unannotated / linker': 32, 'Protein kinase': 28, 'Cadherin-like region 1 (CLD1)': 13, 'Cadherin': 12, 'Cadherin-like region 3 (CLD3)': 6, 'Cadherin-like region 4 (CLD4)': 4}

## Outputs

- `RET_phase1_flagged_candidates.csv`
- `RET_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `RET_structure_viz.html`
