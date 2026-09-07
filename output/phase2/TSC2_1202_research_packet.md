# TSC2 residue 1202 — research packet

*Compiled from this project's own validated data. For manual follow-up research,
not a finished analysis. Framing reminder: nothing here is a classification —
the ceiling is "candidate for further investigation."*

## Why this residue

Selected from 861 cross-validated candidates as the strongest "open question"
case: 542 residues apart from the nearest known pathogenic residue in the raw
sequence, but only 5.05Å apart in the folded 3D structure (a 107x ratio — the
most extreme in the entire candidate list). High structural confidence
(pLDDT 89). Sits in a region with no named UniProt functional domain — the
"unmapped" middle of a 1807-residue protein.

**The structural story:** the nearest known pathogenic residue is position
1744, which sits inside TSC2's Rap-GAP domain (1531–1758) — the enzyme's
catalytic core. Residue 1202 is far from that domain in sequence, but folds
back to sit right beside it in 3D.

## Identifiers
/Users/matisgelinas/Desktop/Passion Project/output/phase2/TSC2_1202_research_packet.md
- Gene: **TSC2**
- UniProt: **P49815**
- RefSeq transcript: **NM_000548.5**
- Position: **1202** (wild-type = Proline)

| Variant | ClinVar Variation ID | Status (as of this project's data pull) |
|---|---|---|
| c.3605C>T (p.Pro1202Leu) | 4192124 | VUS — 1 submitter |
| c.3604C>T (p.Pro1202Ser) | 3004908 | VUS — multiple submitters, no conflicts; linked to "Tuberous sclerosis 2" |
| c.3605C>G (p.Pro1202Arg) | 65155 | no classification provided |
| c.3605C>A (p.Pro1202His) | 65081 | no classification provided |

## Research checklist

1. **ClinVar** (check for reclassification since data pull, submitter notes):
   - https://www.ncbi.nlm.nih.gov/clinvar/variation/4192124/
   - https://www.ncbi.nlm.nih.gov/clinvar/variation/3004908/
2. **gnomAd** (population frequency — do this early):
   https://gnomad.broadinstitute.org — search `TSC2 P1202L` / `TSC2 P1202S`
3. **Ensembl VEP** (AlphaMissense + REVEL + other predictors):
   https://www.ensembl.org/Homo_sapiens/Tools/VEP
   Input: `NM_000548.5:c.3605C>T` and `NM_000548.5:c.3604C>T`
4. **AlphaFold DB** (structure + AlphaMissense pathogenicity track):
   https://alphafold.ebi.ac.uk/entry/P49815
5. **MaveDB** (existing functional assay / deep mutational scanning data):
   https://www.mavedb.org — search "TSC2"
6. **PubMed**:
   https://pubmed.ncbi.nlm.nih.gov — try `TSC2 Pro1202`, `TSC2 P1202L`,
   `TSC2 Rap-GAP domain structure`
7. **UniProt**:
   https://www.uniprot.org/uniprotkb/P49815

## This project's own data on this residue

- Domain: Unannotated / linker (no named UniProt feature covers position 1202)
- pLDDT (structural confidence): 89.0
- Nearest pathogenic residue: 1744 (in the Rap-GAP domain)
- 3D distance: 5.05 Å
- Sequence distance: 542 residues
- Phase 2 model probability (Pathogenic-like): 0.89
- Flagged by both the geometric rule and the cross-validated model
