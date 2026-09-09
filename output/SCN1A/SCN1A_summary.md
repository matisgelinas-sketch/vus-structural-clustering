# SCN1A — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to SCN1A, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 2420 germline missense variants: 1037
  Pathogenic/Likely pathogenic, 63 Benign/Likely benign,
  1320 VUS.
- Structure: AlphaFold model, 2009 residues
  (matches UniProt canonical
  length 2009).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

642 of 2357
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 270 flagged candidates,
20 are low-confidence.

## Results

**270 of 1320 VUS (20.5%)**
flagged as candidates.

By domain: {'III': 71, 'IV': 66, 'I': 60, 'II': 42, 'Unannotated / linker': 28, 'Disordered': 3}

## Outputs

- `SCN1A_phase1_flagged_candidates.csv`
- `SCN1A_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `SCN1A_structure_viz.html`
