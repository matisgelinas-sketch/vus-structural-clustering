"""
Orchestrates Steps 1-3 of the gene-expansion task: for each candidate
gene, fetch ClinVar germline missense data, apply the 30/30 data bar,
and for genes that clear it, fetch+validate structure/sequence and save
everything into Genes/{GENE}/. Processes a slice of the candidate list
(so it can be run in time-bounded chunks) and appends results to a JSON
log for later aggregation.

Usage: python3 run_gene_expansion.py <start_idx> <end_idx>
"""
import sys
import json
import traceback
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent  # project root

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
from ncbi_clinvar import (
    fetch_clinvar_germline_missense, meets_data_bar, screen_and_fetch_gene,
    fetch_domain_features,
)
from clustering_lib import cross_check_numbering

CANDIDATES = [
    ("APC", "P25054"), ("BMPR1A", "P36894"), ("CDH1", "P12830"), ("MEN1", "O00255"),
    ("MLH1", "P40692"), ("MSH2", "P43246"), ("MSH6", "P52701"), ("PMS2", "P54278"),
    ("MUTYH", "Q9UIF7"), ("NF2", "P35240"), ("PALB2", "Q86YC2"), ("RB1", "P06400"),
    ("RET", "P07949"), ("SDHAF2", "Q9NX18"), ("SDHB", "P21912"), ("SDHC", "Q99643"),
    ("SDHD", "O14521"), ("SMAD4", "Q13485"), ("STK11", "Q15831"), ("TSC1", "Q92574"),
    ("TSC2", "P49815"), ("VHL", "P40337"), ("WT1", "P19544"), ("ACTA2", "P62736"),
    ("ACTC1", "P68032"), ("CACNA1S", "Q13698"), ("CASQ2", "O14958"), ("COL3A1", "P02461"),
    ("DES", "P17661"), ("DSC2", "Q02487"), ("DSG2", "Q14126"), ("DSP", "P15924"),
    ("FBN1", "P35555"), ("FLNC", "Q14315"), ("GLA", "P06280"), ("KCNH2", "Q12809"),
    ("KCNQ1", "P51787"), ("LDLR", "P01130"), ("LMNA", "P02545"), ("MYBPC3", "Q14896"),
    ("MYH11", "P35749"), ("MYH7", "P12883"), ("MYL2", "P10916"), ("MYL3", "P08590"),
    ("PCSK9", "Q8NBP7"), ("PKP2", "Q99959"), ("PRKAG2", "Q9UGJ0"), ("SCN5A", "Q14524"),
    ("TGFBR1", "P36897"), ("TGFBR2", "P37173"), ("TMEM43", "Q9BTV4"), ("TNNC1", "P63316"),
    ("TNNI3", "P19429"), ("TNNT2", "P45379"), ("TPM1", "P09493"), ("TRDN", "Q13061"),
    ("ATP7B", "P35670"), ("BTD", "P43251"), ("CASR", "P41180"), ("GAA", "P10253"),
    ("HNF1A", "P20823"), ("OTC", "P00480"), ("TTR", "P02766"), ("ACVRL1", "P37023"),
    ("ENG", "P17813"),
]

GENES_BASE = (PROJECT_ROOT / "Genes")
LOG_PATH = (PROJECT_ROOT / "scripts/phase2b/expansion_log.json")

start = int(sys.argv[1])
end = int(sys.argv[2])
slice_ = CANDIDATES[start:end]

log = json.loads(LOG_PATH.read_text()) if LOG_PATH.exists() else {}

for gene, accession in slice_:
    print(f"\n=== {gene} ({accession}) ===", flush=True)
    try:
        clinvar_df, dropped = fetch_clinvar_germline_missense(gene, gene)
        bar_ok, n_patho, n_benign = meets_data_bar(clinvar_df, min_per_class=30)
        n_vus = int((clinvar_df['bucket'] == 'VUS').sum()) if not clinvar_df.empty else 0
        print(f"  ClinVar: Pathogenic={n_patho} Benign={n_benign} VUS={n_vus} "
              f"dropped={dropped} bar_ok={bar_ok}", flush=True)

        if not bar_ok:
            log[gene] = dict(accession=accession, n_pathogenic=n_patho, n_benign=n_benign,
                              n_vus=n_vus, bar_passed=False, structure_passed=None,
                              reason=f"Failed 30/30 bar (Pathogenic={n_patho}, Benign={n_benign})")
            LOG_PATH.write_text(json.dumps(log, indent=2))
            continue

        struct_result = screen_and_fetch_gene(gene, gene, accession, GENES_BASE)
        if not struct_result.get("passed"):
            print(f"  STRUCTURE FAILED: {struct_result.get('reason')}", flush=True)
            log[gene] = dict(accession=accession, n_pathogenic=n_patho, n_benign=n_benign,
                              n_vus=n_vus, bar_passed=True, structure_passed=False,
                              reason=struct_result.get("reason"))
            LOG_PATH.write_text(json.dumps(log, indent=2))
            continue

        seq = struct_result["seq"]
        struct_df = struct_result["struct_df"]
        af_version = struct_result["af_version"]

        mismatches = cross_check_numbering(clinvar_df, seq, struct_df)
        if len(mismatches):
            print(f"  NUMBERING MISMATCH: {len(mismatches)} rows -- excluding gene", flush=True)
            log[gene] = dict(accession=accession, n_pathogenic=n_patho, n_benign=n_benign,
                              n_vus=n_vus, bar_passed=True, structure_passed=False,
                              reason=f"{len(mismatches)} ClinVar/UniProt/structure numbering mismatches")
            LOG_PATH.write_text(json.dumps(log, indent=2))
            continue

        # --- all checks passed: save everything ---
        gene_dir = GENES_BASE / gene
        gene_dir.mkdir(parents=True, exist_ok=True)

        (gene_dir / f"{gene}_sequence.fasta").write_text(f">sp|{accession}|{gene}_HUMAN\n{seq}\n")

        pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{af_version}.pdb"
        import urllib.request
        with urllib.request.urlopen(pdb_url, timeout=60) as r:
            (gene_dir / f"AF-{accession}-F1-model_v{af_version}.pdb").write_bytes(r.read())

        clinvar_df.to_csv(gene_dir / f"{gene}_clinvar_germline_missense.csv", index=False)

        domain_feats = fetch_domain_features(accession)

        print(f"  PASSED -- saved to {gene_dir}, {len(domain_feats)} domain features", flush=True)
        log[gene] = dict(accession=accession, n_pathogenic=n_patho, n_benign=n_benign,
                          n_vus=n_vus, bar_passed=True, structure_passed=True,
                          reason=None, domain_features=domain_feats,
                          af_version=af_version, length=len(seq))
        LOG_PATH.write_text(json.dumps(log, indent=2))

    except Exception as e:
        print(f"  ERROR: {e}", flush=True)
        traceback.print_exc()
        log[gene] = dict(accession=accession, bar_passed=None, structure_passed=None,
                          reason=f"Unhandled error: {e}")
        LOG_PATH.write_text(json.dumps(log, indent=2))

print(f"\nDone with slice [{start}:{end}]. Log has {len(log)} entries.", flush=True)
