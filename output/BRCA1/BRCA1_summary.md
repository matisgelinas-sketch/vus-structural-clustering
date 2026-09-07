# BRCA1 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready. Prepared as a starting point
for discussion with a researcher, not a finished report.*

## What this is, and isn't

Same framing as the TP53 write-up: this is a **hypothesis-generation
exercise**, not a classification tool. Structural proximity on a
predicted model is the lowest tier of ACMG/AMP evidence and cannot by
itself justify a classification change. Every flagged variant below is a
**candidate for further investigation**, not a verdict.

## Data and methodology notes (same pipeline as TP53)

- Source: local ClinVar TSV export, filtered to germline classifications
  only (`G:` prefix), restricted to rows whose "Variation" transcript
  resolves to BRCA1 itself.
- Protein-level changes parsed from the ClinVar "Variation" column, not
  the "Protein change" column (validated during pre-flight: the
  aggregate column mixes in numbering from unrelated transcripts —
  BRCA1's version of this was especially severe, listing up to ~76
  different residue numbers per row).
- 2330 germline missense rows kept: 120 Pathogenic/Likely pathogenic,
  360 Benign/Likely benign, 1832 VUS, 18 unrecognized label. Dropped: 59
  unparseable, 1070 non-germline/blank (BRCA1 has a much higher fraction
  of somatic-only ClinVar submissions than TP53 or PTEN), 8
  stop-gained/synonymous, 0 wrong-gene-transcript.
- Structure: canonical AlphaFold model `AF-P38398-F1-model_v6.pdb`
  (1863 residues, matches UniProt P38398 exactly). Note: the originally
  provided structure file was for UniProt isoform 2 (63 residues only,
  wrong protein for this purpose) — this was caught during pre-flight
  and replaced with the canonical full-length model from AlphaFold DB.
- Residue numbering cross-checked across ClinVar, UniProt, and structure
  for all 2330 rows: **zero mismatches**.
- Same 3D/sequence distance definitions as TP53. Flagging threshold:
  **≤ 6 Å in 3D AND > 10 residues apart in sequence** — the sequence-
  separation requirement was added after TP53's initial ≤6 Å-only pass
  showed 91% of flagged cases were just chain-adjacent (trivial), not
  genuine tertiary clustering; applied identically here for consistency.

## Structural confidence caveat — bigger factor here than in TP53

**1454 of 1952 Pathogenic+VUS residues (74.5%) fall in AlphaFold
low-confidence regions (pLDDT < 70)** — far more than TP53's 30%. This
is expected biology, not a modeling failure: BRCA1 is mostly an
intrinsically disordered ~1550-residue chain flanking two small folded
domains — an N-terminal RING-type zinc finger (residues 24–65) and
tandem BRCT domains at the C-terminus (1642–1736, 1756–1855). Most of
the protein between those domains has no stable fold to predict
confidently.

The flagged candidates don't inherit that problem proportionally: of 34
flagged VUS, only 3 (9%) are low-confidence — because both pathogenic
residues and their close-in-3D VUS neighbors concentrate almost
entirely in the folded RING/BRCT regions, where pLDDT is high.
Low-confidence residues are kept in every table, explicitly flagged,
not excluded.

## Results

- **34 of 1832 VUS (1.9%)** flagged as candidates under the dual
  criterion (≤6 Å in 3D, >10 residues apart in sequence) — a much
  smaller fraction than TP53's 3%, consistent with BRCA1 being an
  extended multi-domain protein rather than one compact folded domain.
- By domain: 15 in BRCT domain 2, 7 in BRCT domain 1, 2 in the RING
  domain, 10 in unannotated/linker regions immediately flanking those
  domains.
- Ranking is by (sequence distance / 3D distance) ratio — prioritizing
  cases where 3D proximity isn't just linear adjacency.

## Outputs

- `BRCA1_phase1_flagged_candidates.csv` — 34 flagged candidates
- `BRCA1_all_VUS_annotated.csv` — all 1832 VUS, with a `flagged_candidate`
  column, flagged or not (nothing silently dropped)
- `numbering_mismatches.csv` — empty for BRCA1
- `low_confidence_flagged.csv` — all 1454 low-pLDDT Pathogenic/VUS rows
- `BRCA1_structure_viz.html` — interactive 3D viewer (offline-capable),
  red/green/blue scheme, same as TP53
- `BRCA1_summary.md` — this document

## Limitations to carry into any future write-up

- Same general limitations as TP53 (single static model, no
  substitution-severity feature yet, distance clustering doesn't
  distinguish functional-pocket proximity from packing coincidence).
- BRCA1-specific: the very large disordered fraction means most of the
  protein is structurally unassessable by this method at all — this
  pipeline is really only informative for BRCA1 variants in or very
  near the RING and BRCT domains. Anything flagged deep in the
  disordered linker (low-confidence, if any) should be treated as
  close to uninterpretable structurally.
- 1070 dropped ClinVar rows had no germline classification at all
  (likely somatic-only submissions) — worth being explicit that this
  analysis says nothing about those variants.
