# MYO7A — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to MYO7A, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1568 germline missense variants: 176
  Pathogenic/Likely pathogenic, 39 Benign/Likely benign,
  1353 VUS.
- Structure: AlphaFold model, 2215 residues
  (matches UniProt canonical
  length 2215).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

300 of 1529
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 132 flagged candidates,
1 are low-confidence.

## Results

**132 of 1353 VUS (9.8%)**
flagged as candidates.

By domain: {'Myosin motor': 65, 'FERM 1': 15, 'Unannotated / linker': 14, 'FERM 2': 13, 'MyTH4 2': 12, 'MyTH4 1': 8, 'SH3': 3, 'SAH': 2}

## Outputs

- `MYO7A_phase1_flagged_candidates.csv`
- `MYO7A_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `MYO7A_structure_viz.html`
