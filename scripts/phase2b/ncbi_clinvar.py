"""
Reusable ClinVar-via-NCBI-E-utilities loader for gene-expansion work.
Produces the same shape of data as clustering_lib.parse_clinvar() did for
the manually-exported TSVs (TP53/BRCA1/PTEN), but sourced programmatically
via esearch + esummary, so no manual website export is needed for new genes.

Validated fields (checked against a live esummary response for VHL,
2026-08-17): title (same "NM_...(GENE):c...(p.Xxx##Yyy)" format as the
TSV export's "Variation" column), germline_classification.description
(cleanly separated from somatic_clinical_impact / oncogenicity_classification
-- no G:/S:/O: prefix parsing needed here, NCBI already separates them).
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from clustering_lib import AA3TO1, PATHOGENIC_LABELS, BENIGN_LABELS, VUS_LABELS

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

# NCBI: ~3 req/sec without an API key, ~10 req/sec with one. Get a free
# key at https://www.ncbi.nlm.nih.gov/account/settings/ ("API Key
# Management") and either export NCBI_API_KEY or pass --api-key to
# gene_pipeline.py -- roughly 3x faster fetching for free.
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
SLEEP = 0.11 if NCBI_API_KEY else 0.34

PROTEIN_CHANGE_RE = re.compile(r'p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})')


def uniprot_lookup(gene_symbol: str):
    """Resolve a gene symbol to (accession, canonical_length) via UniProt
    search, so candidate genes only need to be named, not pre-researched."""
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" +
           urllib.parse.quote(f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true") +
           "&format=json&fields=accession,length")
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            d = json.load(r)
        if not d.get("results"):
            return None, None
        res = d["results"][0]
        return res["primaryAccession"], res["sequence"]["length"]
    except Exception:
        return None, None


def _with_key(params):
    if NCBI_API_KEY and "eutils.ncbi.nlm.nih.gov" in ESEARCH:
        params = dict(params)
        params["api_key"] = NCBI_API_KEY
    return params


def _get_json(url, params, retries=3):
    if "eutils.ncbi.nlm.nih.gov" in url:
        params = _with_key(params)
    full = url + "?" + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(full, timeout=30) as r:
                return json.load(r)
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(1.5 * (attempt + 1))


def _post_json(url, data, retries=3):
    if "eutils.ncbi.nlm.nih.gov" in url:
        data = _with_key(data)
    body = urllib.parse.urlencode(data).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=body, method="POST")
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(1.5 * (attempt + 1))


def fetch_clinvar_germline_missense(gene: str, gene_tag: str):
    """
    Equivalent of clustering_lib.parse_clinvar(), sourced from NCBI
    E-utilities instead of a manual TSV export. Same dropped-counts
    shape, same bucketing rules (germline-only, transcript-matched).
    Returns (DataFrame, dropped_dict).
    """
    dropped = {
        'wrong_gene_transcript': 0,
        'no_protein_change_parsed': 0,
        'non_germline_or_blank': 0,
        'stop_gained_or_synonymous': 0,
    }

    search = _get_json(ESEARCH, {
        "db": "clinvar",
        "term": f"{gene}[gene] AND missense_variant[Molecular consequence]",
        "retmode": "json",
        "retmax": 100000,
    })
    time.sleep(SLEEP)
    uids = search["esearchresult"]["idlist"]

    rows = []
    batch_size = 250
    for i in range(0, len(uids), batch_size):
        batch = uids[i:i + batch_size]
        summ = _post_json(ESUMMARY, {"db": "clinvar", "id": ",".join(batch), "retmode": "json"})
        time.sleep(SLEEP)
        result = summ.get("result", {})
        for uid in result.get("uids", batch):
            rec = result.get(uid)
            if rec is None:
                continue
            title = rec.get("title", "")
            if f'({gene_tag})' not in title:
                dropped['wrong_gene_transcript'] += 1
                continue
            m = PROTEIN_CHANGE_RE.search(title)
            if not m:
                dropped['no_protein_change_parsed'] += 1
                continue
            wt3, pos, mt3 = m.group(1), int(m.group(2)), m.group(3)
            if wt3 not in AA3TO1 or mt3 not in AA3TO1:
                dropped['no_protein_change_parsed'] += 1
                continue
            if mt3 == 'Ter' or wt3 == mt3:
                dropped['stop_gained_or_synonymous'] += 1
                continue

            germline = rec.get("germline_classification") or {}
            classification = (germline.get("description") or "").strip()
            if not classification:
                dropped['non_germline_or_blank'] += 1
                continue

            if classification in PATHOGENIC_LABELS:
                bucket = 'Pathogenic'
            elif classification in BENIGN_LABELS:
                bucket = 'Benign'
            elif classification in VUS_LABELS:
                bucket = 'VUS'
            else:
                bucket = None

            rows.append({
                'variation': title,
                'position': pos,
                'wt_aa': AA3TO1[wt3],
                'mt_aa': AA3TO1[mt3],
                'classification_raw': classification,
                'bucket': bucket,
                'variation_id': rec.get('uid', ''),
                'review_status': germline.get('review_status', ''),
                'condition': "; ".join(
                    t.get('trait_name', '') for t in germline.get('trait_set', [])
                ) if germline.get('trait_set') else '',
            })

    df = pd.DataFrame(rows)
    return df, dropped


def meets_data_bar(df: pd.DataFrame, min_per_class: int = 30):
    if df.empty:
        return False, 0, 0
    n_patho = (df['bucket'] == 'Pathogenic').sum()
    n_benign = (df['bucket'] == 'Benign').sum()
    return (n_patho >= min_per_class and n_benign >= min_per_class), int(n_patho), int(n_benign)


def screen_and_fetch_gene(gene: str, gene_tag: str, accession: str, out_base: Path):
    """
    Step 2: for a gene that already cleared the ClinVar data bar, fetch
    UniProt sequence + AlphaFold structure, cross-validate all three
    numbering systems, fetch domain annotations, and save everything
    into Genes/{gene}/ in the same layout as TP53/BRCA1/PTEN.

    Returns a dict describing the outcome (success/fail + reason +
    counts), and on success also returns the parsed clinvar_df,
    struct_df, seq, and domain feature list so the caller doesn't have
    to re-fetch/re-parse.
    """
    from Bio.PDB import PDBParser
    import io

    result = {"gene": gene, "accession": accession}

    # --- UniProt canonical sequence ---
    uni = _get_json(f"https://rest.uniprot.org/uniprotkb/{accession}.json", {})
    time.sleep(SLEEP)
    seq = uni["sequence"]["value"]
    uni_length = uni["sequence"]["length"]

    # --- AlphaFold prediction metadata ---
    try:
        af_meta_req = urllib.request.Request(
            f"https://alphafold.ebi.ac.uk/api/prediction/{accession}")
        with urllib.request.urlopen(af_meta_req, timeout=30) as r:
            af_meta = json.load(r)
        time.sleep(SLEEP)
    except Exception as e:
        result.update(passed=False, reason=f"AlphaFold API lookup failed: {e}")
        return result

    if not af_meta:
        result.update(passed=False, reason="No AlphaFold model available")
        return result

    af_seq = af_meta[0]["sequence"]
    af_version = af_meta[0]["latestVersion"]
    if len(af_seq) != uni_length:
        result.update(passed=False,
                       reason=f"AlphaFold sequence length ({len(af_seq)}) != "
                              f"UniProt canonical length ({uni_length}) -- "
                              f"likely fragment/isoform mismatch")
        return result

    # --- download + parse structure ---
    pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{af_version}.pdb"
    try:
        with urllib.request.urlopen(pdb_url, timeout=60) as r:
            pdb_text = r.read().decode()
        time.sleep(SLEEP)
    except Exception as e:
        result.update(passed=False, reason=f"PDB download failed: {e}")
        return result

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('model', io.StringIO(pdb_text))
    model = next(structure.get_models())
    chain = next(model.get_chains())
    struct_rows = []
    for residue in chain:
        if 'CA' not in residue:
            continue
        ca = residue['CA']
        resname3 = residue.get_resname().capitalize()
        struct_rows.append({
            'position': residue.id[1],
            'resname_1letter': AA3TO1.get(resname3, 'X'),
            'x': ca.coord[0], 'y': ca.coord[1], 'z': ca.coord[2],
            'plddt': ca.get_bfactor(),
        })
    struct_df = pd.DataFrame(struct_rows).sort_values('position').reset_index(drop=True)

    if len(struct_df) != uni_length:
        result.update(passed=False,
                       reason=f"Structure Ca count ({len(struct_df)}) != "
                              f"UniProt canonical length ({uni_length})")
        return result

    # --- fetch ClinVar again (caller already validated the bar, but we
    # need the df here too; caller passes it in to avoid double-fetch --
    # see fetch_and_validate_gene below which composes both steps) ---
    result.update(passed=True, seq=seq, struct_df=struct_df, af_version=af_version,
                   pdb_text=pdb_text)
    return result


def process_gene(gene: str, log: dict, genes_base: Path, log_path: Path,
                  min_per_class: int = 30, max_length: int = 3000, accession: str = None):
    """
    Full pipeline for one gene: skip if already logged, resolve its
    UniProt accession if not given, screen by size, fetch+bar-check
    ClinVar, fetch+validate structure, cross-check numbering, save to
    Genes/{gene}/, fetch domain features. Writes `log` to `log_path`
    after this gene regardless of outcome (checkpointing), so a run
    interrupted mid-list loses at most the one gene in progress.

    Returns the log entry dict for this gene (also stored in `log[gene]`).
    """
    from clustering_lib import cross_check_numbering, AA3TO1

    if gene in log:
        return log[gene]  # already processed in a prior run -- resumable for free

    entry = {"gene": gene}
    try:
        if accession is None:
            accession, length = uniprot_lookup(gene)
            if accession is None:
                entry.update(bar_passed=None, structure_passed=None,
                              reason="UniProt lookup failed (bad gene symbol?)")
                log[gene] = entry
                log_path.write_text(json.dumps(log, indent=2))
                return entry
        else:
            _, length = uniprot_lookup(gene)

        entry["accession"] = accession
        if length and length > max_length:
            entry.update(bar_passed=None, structure_passed=False, length=length,
                          reason=f"Too large for single-chain AlphaFold model ({length} > {max_length} aa)")
            log[gene] = entry
            log_path.write_text(json.dumps(log, indent=2))
            return entry

        clinvar_df, dropped = fetch_clinvar_germline_missense(gene, gene)
        bar_ok, n_patho, n_benign = meets_data_bar(clinvar_df, min_per_class=min_per_class)
        n_vus = int((clinvar_df['bucket'] == 'VUS').sum()) if not clinvar_df.empty else 0
        entry.update(n_pathogenic=n_patho, n_benign=n_benign, n_vus=n_vus)

        if not bar_ok:
            entry.update(bar_passed=False, structure_passed=None,
                          reason=f"Failed {min_per_class}/{min_per_class} bar "
                                 f"(Pathogenic={n_patho}, Benign={n_benign})")
            log[gene] = entry
            log_path.write_text(json.dumps(log, indent=2))
            return entry

        struct_result = screen_and_fetch_gene(gene, gene, accession, genes_base)
        if not struct_result.get("passed"):
            entry.update(bar_passed=True, structure_passed=False,
                          reason=struct_result.get("reason"))
            log[gene] = entry
            log_path.write_text(json.dumps(log, indent=2))
            return entry

        seq = struct_result["seq"]
        struct_df = struct_result["struct_df"]
        af_version = struct_result["af_version"]
        pdb_text = struct_result["pdb_text"]

        mismatches = cross_check_numbering(clinvar_df, seq, struct_df)
        if len(mismatches):
            entry.update(bar_passed=True, structure_passed=False,
                          reason=f"{len(mismatches)} ClinVar/UniProt/structure numbering mismatches")
            log[gene] = entry
            log_path.write_text(json.dumps(log, indent=2))
            return entry

        gene_dir = genes_base / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        (gene_dir / f"{gene}_sequence.fasta").write_text(f">sp|{accession}|{gene}_HUMAN\n{seq}\n")
        (gene_dir / f"AF-{accession}-F1-model_v{af_version}.pdb").write_text(pdb_text)
        clinvar_df.to_csv(gene_dir / f"{gene}_clinvar_germline_missense.csv", index=False)

        domain_feats = fetch_domain_features(accession)
        entry.update(bar_passed=True, structure_passed=True, reason=None,
                      domain_features=domain_feats, af_version=af_version, length=len(seq))
        log[gene] = entry
        log_path.write_text(json.dumps(log, indent=2))
        return entry

    except Exception as e:
        entry.update(bar_passed=entry.get("bar_passed"), structure_passed=None,
                      reason=f"Unhandled error: {e}")
        log[gene] = entry
        log_path.write_text(json.dumps(log, indent=2))
        return entry


def fetch_domain_features(accession: str):
    uni = _get_json(f"https://rest.uniprot.org/uniprotkb/{accession}.json", {})
    time.sleep(SLEEP)
    feats = []
    for feat in uni.get('features', []):
        if feat['type'] in ('Domain', 'Region', 'Zinc finger', 'DNA binding', 'Repeat'):
            loc = feat['location']
            start = loc['start']['value']
            end = loc['end']['value']
            desc = feat.get('description', '') or feat['type']
            feats.append((start, end, desc))
    return feats
