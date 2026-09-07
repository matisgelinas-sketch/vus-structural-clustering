"""
Stage 1 (fast prescreen): for each ACMG SF v3.2-derived candidate gene,
check UniProt canonical length (<=3000 for AlphaFold single-chain
modeling) and a rough ClinVar missense-variant hit count (cheap, 1
esearch call per gene, no per-variant fetching yet). Genes that clear
both go on to stage 2 (exact Pathogenic/Benign counting).
"""
import time
import json
import urllib.request
import urllib.parse

CANDIDATES = [
    # Cancer predisposition (ACMG SF v3.2)
    "APC", "BMPR1A", "CDH1", "MEN1", "MLH1", "MSH2", "MSH6", "PMS2",
    "MUTYH", "NF2", "PALB2", "RB1", "RET", "SDHAF2", "SDHB", "SDHC",
    "SDHD", "SMAD4", "STK11", "TSC1", "TSC2", "VHL", "WT1",
    # Cardiovascular (ACMG SF v3.2)
    "ACTA2", "ACTC1", "APOB", "CACNA1S", "CASQ2", "COL3A1", "DES",
    "DSC2", "DSG2", "DSP", "FBN1", "FLNC", "GLA", "KCNH2", "KCNQ1",
    "LDLR", "LMNA", "MYBPC3", "MYH11", "MYH7", "MYL2", "MYL3", "PCSK9",
    "PKP2", "PRKAG2", "RYR1", "RYR2", "SCN5A", "TGFBR1", "TGFBR2",
    "TMEM43", "TNNC1", "TNNI3", "TNNT2", "TPM1", "TRDN",
    # Other ACMG SF v3.2
    "ATP7B", "BTD", "CASR", "GAA", "HNF1A", "OTC", "TTR", "ACVRL1", "ENG",
]

def uniprot_length(gene):
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" +
           urllib.parse.quote(f"gene:{gene} AND organism_id:9606 AND reviewed:true") +
           "&format=json&fields=accession,length,gene_names")
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            d = json.load(r)
        if not d.get("results"):
            return None, None
        res = d["results"][0]
        return res["primaryAccession"], res["sequence"]["length"]
    except Exception as e:
        return None, f"ERROR: {e}"

def clinvar_missense_count(gene):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=" +
           urllib.parse.quote(f"{gene}[gene] AND missense_variant[Molecular consequence]") +
           "&retmode=json&retmax=0")
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            d = json.load(r)
        return int(d["esearchresult"]["count"])
    except Exception as e:
        return f"ERROR: {e}"

print(f"{'Gene':<10} {'Accession':<12} {'Length':<8} {'MissenseHits':<14} Verdict")
passed = []
for gene in CANDIDATES:
    acc, length = uniprot_length(gene)
    time.sleep(0.34)  # be polite to UniProt
    hits = clinvar_missense_count(gene)
    time.sleep(0.34)  # NCBI: max ~3 req/sec without API key

    if acc is None:
        print(f"{gene:<10} {'N/A':<12} {'N/A':<8} {'N/A':<14} SKIP (UniProt lookup failed)")
        continue
    if isinstance(hits, str):
        print(f"{gene:<10} {acc:<12} {length:<8} {'N/A':<14} SKIP (ClinVar lookup failed: {hits})")
        continue

    size_ok = length <= 3000
    data_ok = hits >= 80  # rough proxy floor; exact P/B counts checked in stage 2
    verdict = "PASS -> stage 2" if (size_ok and data_ok) else \
              ("FAIL: too large" if not size_ok else "FAIL: too few missense hits")
    print(f"{gene:<10} {acc:<12} {length:<8} {hits:<14} {verdict}")
    if size_ok and data_ok:
        passed.append((gene, acc, length, hits))

print(f"\n{len(passed)} / {len(CANDIDATES)} candidates passed prescreen:")
for gene, acc, length, hits in passed:
    print(f"  {gene} ({acc}, {length}aa, {hits} missense hits)")
