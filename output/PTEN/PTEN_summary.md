# PTEN — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready. Prepared as a starting point
for discussion with a researcher, not a finished report.*

## What this is, and isn't

Same framing as TP53/BRCA1: **hypothesis-generation only**. Structural
proximity on a predicted model is the lowest tier of ACMG/AMP evidence
and cannot by itself justify a classification change. Every flagged
variant below is a **candidate for further investigation**, not a
verdict.

## Data and methodology notes (same pipeline as TP53/BRCA1)

- Source: local ClinVar TSV export, filtered to germline classifications
  only (`G:` prefix), restricted to rows whose "Variation" transcript
  resolves to PTEN itself — this drops 4 rows that were actually **KLLN**
  variants (a gene overlapping PTEN's promoter region), caught and
  validated during pre-flight.
- Protein-level changes parsed from the "Variation" column, not the
  "Protein change" column, for the same reason as TP53/BRCA1.
- 898 germline missense rows kept: 180 Pathogenic/Likely pathogenic,
  **3 Benign/Likely benign**, 715 VUS. Dropped: 95 unparseable, 9
  non-germline/blank, 4 wrong-gene-transcript (KLLN).
- Structure: AlphaFold model `AF-P60484-F1-model_v6.pdb` (403 residues,
  matches UniProt P60484 exactly).
- Residue numbering cross-checked across ClinVar, UniProt, and structure
  for all 898 rows: **zero mismatches**.
- Same distance definitions as the other two genes. Flagging threshold:
  **≤ 6 Å in 3D AND > 10 residues apart in sequence**. A first pass with
  ≤6 Å alone flagged 408/715 VUS (57%) — even more than TP53's original
  8–10 Å-only result — because PTEN is an even more compact,
  single-domain protein, so a distance-only rule barely discriminates.
  Checking showed the same pattern as TP53: the great majority of those
  were just chain-adjacent to a pathogenic residue, not genuine tertiary
  clustering. Adding the >10-residue sequence-separation requirement
  (confirmed with you, applied identically across all three genes) fixes
  this.

## Worth flagging honestly

**Almost no germline Benign variants (3 total).** PTEN's ClinVar
submissions skew heavily toward Pathogenic/VUS, likely because PTEN
testing happens mostly in a clinical (Cowden syndrome/PHTS) context
where benign findings are less often submitted. This doesn't affect
Phase 1's distance-flagging (which only uses Pathogenic residues as
reference points), but it's a real limitation for Phase 2's pooled
training — PTEN will contribute almost no Benign examples to the
combined model, and any gene-held-out validation using PTEN as the test
gene will have essentially no Benign cases to evaluate against.

## Structural confidence

146 of 895 Pathogenic+VUS residues (16.3%) fall in low-confidence
regions (pLDDT<70) — almost entirely the disordered C-terminal tail
(residues ~352–403). None of the 47 flagged candidates are
low-confidence.

## Results

- **47 of 715 VUS (6.6%)** flagged as candidates under the dual
  criterion (≤6 Å in 3D, >10 residues apart in sequence).
- By domain: 35 in the C2 domain, 12 in the phosphatase domain.
- Ranked by (sequence distance / 3D distance) ratio, same as the other
  two genes.

## Outputs

- `PTEN_phase1_flagged_candidates.csv` — 47 flagged candidates
- `PTEN_all_VUS_annotated.csv` — all 715 VUS, with a `flagged_candidate`
  column, flagged or not (nothing silently dropped)
- `numbering_mismatches.csv` — empty for PTEN
- `low_confidence_flagged.csv` — 146 low-pLDDT Pathogenic/VUS rows
- `PTEN_structure_viz.html` — interactive 3D viewer (offline-capable),
  red/green/blue scheme, same as the other two genes
- `PTEN_summary.md` — this document

## Limitations to carry into any future write-up

- Same general limitations as TP53/BRCA1 (single static model, no
  substitution-severity feature, distance clustering doesn't distinguish
  functional relevance from packing coincidence).
- PTEN-specific: near-total absence of germline Benign labels limits how
  much this gene can teach Phase 2's classifier about what "not
  clustered" looks like from confirmed-benign examples.
- BRCA2 is excluded project-wide: no usable single-structure AlphaFold
  model is available due to protein size, which should be noted as a
  scope limitation in any final write-up.
