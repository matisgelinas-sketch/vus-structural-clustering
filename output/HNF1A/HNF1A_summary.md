# HNF1A — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to HNF1A, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 510 germline missense variants: 143
  Pathogenic/Likely pathogenic, 35 Benign/Likely benign,
  332 VUS.
- Structure: AlphaFold model, 631 residues
  (matches UniProt canonical
  length 631).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

205 of 475
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 28 flagged candidates,
2 are low-confidence.

## Results

**28 of 332 VUS (8.4%)**
flagged as candidates.

By domain: {'Homeobox; HNF1-type': 19, 'POU-specific atypical': 9}

## Outputs

- `HNF1A_phase1_flagged_candidates.csv`
- `HNF1A_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `HNF1A_structure_viz.html`
