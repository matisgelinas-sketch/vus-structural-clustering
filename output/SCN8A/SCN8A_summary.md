# SCN8A — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready.*

## What this is, and isn't

Same framing as TP53/BRCA1/PTEN: **hypothesis-generation only**. Every
flagged variant is a **candidate for further investigation**, not a
verdict.

## Data and methodology

- Source: NCBI ClinVar via E-utilities (germline-only, transcript-matched
  to SCN8A, protein change parsed from the variant title) — same
  filtering methodology as the original three genes, sourced
  programmatically this round.
- 1313 germline missense variants: 284
  Pathogenic/Likely pathogenic, 39 Benign/Likely benign,
  990 VUS.
- Structure: AlphaFold model, 1980 residues
  (matches UniProt canonical
  length 1980).
- Residue numbering cross-checked across ClinVar, UniProt, and
  structure: **0 mismatches** found.
- Flagging threshold: **≤ 6 Å (backbone Cα *or* side-chain Cβ distance)
  AND > 10 residues apart in sequence** (same dual criterion validated on
  TP53, to exclude trivial chain-adjacent cases; computing both atom
  types catches cases where a residue's side chain, not its backbone,
  is what's actually close — a candidate confirmed by both is stronger
  evidence than one found by only one).

## Structural confidence

539 of 1274
Pathogenic+VUS residues fall in low-confidence regions (pLDDT<70). Kept
in every table, not excluded. Of the 69 flagged candidates,
5 are low-confidence.

## Results

**69 of 990 VUS (7.0%)**
flagged as candidates.

By domain: {'IV': 22, 'I': 18, 'II': 17, 'III': 8, 'Unannotated / linker': 4}

## Outputs

- `SCN8A_phase1_flagged_candidates.csv`
- `SCN8A_all_VUS_annotated.csv`
- `numbering_mismatches.csv`
- `low_confidence_flagged.csv`
- `SCN8A_structure_viz.html`
