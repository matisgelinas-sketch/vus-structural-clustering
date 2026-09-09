# LDLR — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to LDLR, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1637 germline missense variants: 739
  Pathogenic/Likely pathogenic, 41 Benign/Likely benign,
  857 VUS.
- Structure: AlphaFold model, 860 residues
  (matches UniProt canonical
  length 860).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

284 of 1596
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 322 flagged candidates,
6 are low-confidence.

## Results

**322 of 857 VUS (37.6%)**
flagged as candidates.

By domain: {'LDL-receptor class B 3': 40, 'LDL-receptor class B 6': 39, 'LDL-receptor class B 5': 37, 'LDL-receptor class B 2': 33, 'LDL-receptor class B 1': 29, 'EGF-like 3': 20, 'LDL-receptor class A 7': 20, 'LDL-receptor class B 4': 18, 'EGF-like 2; calcium-binding': 15, 'EGF-like 1': 13, 'LDL-receptor class A 2': 12, 'LDL-receptor class A 1': 12, 'LDL-receptor class A 5': 9, 'LDL-receptor class A 3': 9, 'LDL-receptor class A 4': 7, 'LDL-receptor class A 6': 6, 'Unannotated / linker': 3}

## Outputs

- `LDLR_phase1_flagged_candidates.csv`
- `LDLR_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `LDLR_structure_viz.html`
