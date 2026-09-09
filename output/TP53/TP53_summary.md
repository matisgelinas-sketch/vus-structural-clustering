# TP53 — structural clustering of VUS against known pathogenic residues

*Draft working notes, not submission-ready. Prepared as a starting point
for discussion with a researcher, not a finished report.*

## What this is, and isn't

This is a **hypothesis-generation exercise**: it looks at whether a TP53
variant of uncertain significance (VUS) sits close in 3D space, on an
AlphaFold-predicted structure, to a residue with an existing germline
Pathogenic/Likely pathogenic ClinVar classification. Structural proximity
on a predicted model is the lowest tier of evidence under ACMG/AMP
guidelines and **cannot by itself justify a classification change**.
Nothing in this document should be read as "likely pathogenic,"
"reclassified," or "confirmed." The correct framing for every flagged
variant below is **candidate for further investigation**.

## Data and methodology notes

- Source: local ClinVar TSV export, filtered to germline classifications
  only, restricted to rows whose ClinVar "Variation" transcript resolves
  to TP53 itself (this excludes 2 rows that were actually WRAP53
  variants picked up by ClinVar's genomic-region search).
- Protein-level changes were parsed from the ClinVar "Variation" column
  (e.g. `p.His296Arg`), not the "Protein change" column, which mixes in
  numbering from unrelated transcripts and is not reliable for this use.
- 748 germline missense rows kept: 142 Pathogenic/Likely pathogenic, 97
  Benign/Likely benign, 508 VUS, 1 with an unrecognized classification
  label, 11 dropped as non-germline (somatic/oncogenicity) or blank,
  24 dropped as not parseable as a simple single missense change.
- Structure: AlphaFold model `AF-P04637-F1-model_v6.pdb` (393 residues,
  matches UniProt P04637 canonical length exactly).
- Residue numbering was cross-checked across ClinVar, the UniProt
  canonical sequence, and the AlphaFold structure for all 748 rows:
  **zero mismatches**.
- 3D distance = both **backbone (Cα–Cα)** and **side-chain (Cβ–Cβ,
  glycine falls back to Cα since it has no side chain)** Euclidean
  distance (Å) from each VUS residue to its nearest *qualifying*
  Pathogenic-bucket residue. "Qualifying" is important: the search only
  considers candidates that already satisfy the sequence-separation
  requirement below, so a genuinely close but non-trivial pathogenic
  residue is never masked by a closer, trivially-adjacent one (an
  earlier version of this pipeline had that bug — fixed).
- Flagging threshold: **≤ 6 Å by either Cα or Cβ distance AND > 10
  residues apart in sequence** (tightened from an initial 8–10 Å,
  3D-distance-only starting point — at 8–10 Å alone, 57–62% of all VUS
  were flagged, which turned out to be mostly a trivial artifact: 91% of
  those flagged cases were within 5 residues of their nearest pathogenic
  residue *in sequence*, i.e. simply chain-adjacent, not evidence of
  tertiary-structure clustering. Requiring >10 residues of sequence
  separation removes that noise and keeps only proximity that isn't
  explained by linear adjacency.) Checking both atom types matters
  because a residue's side chain, not its backbone, is often what's
  actually close to a pathogenic residue; a candidate confirmed by
  *both* atom types (`confirmed_by_both_atoms` column) is stronger
  evidence than one found by only one.

## Structural confidence caveat

198 of 650 Pathogenic+VUS residues (30%) fall in AlphaFold
low-confidence regions (pLDDT < 70) — almost entirely the N-terminal
transactivation domain (~residues 1–93) and the disordered C-terminal
tail, which is consistent with known TP53 biology (these regions are
intrinsically disordered, not just poorly predicted). The folded
DNA-binding domain (residues ~94–292), where all but one flagged
candidate below falls, is high-confidence. Low-confidence residues are
**not excluded** from the tables — they're explicitly flagged
(`low_confidence` column) so you can weight them accordingly. None of
the 16 flagged candidates fall in low-confidence regions.

## Results

- **132 of 508 VUS (26%)** flagged as candidates under the dual
  criterion (≤6 Å by Cα or Cβ, >10 residues apart in sequence); 63 of
  those are confirmed by both atom types. This is a substantially
  higher count than an earlier, Cα-only pass (16 flagged) — most of the
  increase traces to a masking-bug fix (see project-level `PHASE2_summary.md`)
  rather than the Cβ addition alone.
- By domain: 127 in the DNA-binding domain, 3 in the TAD2/proline-rich
  linker, 1 in the DBD-tetramerization linker, 1 unannotated. None are
  in low-confidence (pLDDT<70) regions.
- Ranked by (sequence distance / 3D distance) ratio — prioritizing cases
  where 3D proximity isn't explained by simply being near a pathogenic
  residue in the linear sequence.

## Outputs

- `TP53_phase1_flagged_candidates.csv` — the 132 flagged candidates,
  sorted by (sequence distance / 3D distance) ratio, most non-trivial
  clustering first.
- `TP53_all_VUS_annotated.csv` — all 508 VUS with distances, domain, and
  confidence annotations, plus a `flagged_candidate` column, flagged or
  not (nothing silently dropped).
- `numbering_mismatches.csv` — would list any ClinVar/UniProt/structure
  numbering disagreements; empty for TP53.
- `low_confidence_flagged.csv` — all Pathogenic/VUS residues in
  low-pLDDT regions.
- `TP53_structure_viz.html` — interactive 3D viewer (open in a browser;
  no internet connection required). Red = known pathogenic, green =
  flagged VUS candidate in high-confidence structure, blue = flagged VUS
  candidate in a low-confidence region.

## Limitations to carry into any future write-up

- Single static AlphaFold model — no conformational flexibility,
  ligand/DNA-bound states, or oligomeric (tetramer) context considered,
  despite TP53 functioning as a tetramer.
- Distance-based clustering will inherently flag anything in a densely
  packed domain; it does not distinguish "near a pathogenic residue
  because same functional pocket" from "near a pathogenic residue by
  packing coincidence."
- No amino-acid-substitution severity (e.g. side chain size/charge
  change) is factored in — Cβ distance captures side-chain *position*,
  not substitution severity or rotamer flexibility.
- This is one gene's output in isolation; Phase 2 will test whether this
  signal is predictive when pooled and validated across genes.
