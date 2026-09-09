"""
Shared parsing/analysis helpers for the VUS structural-clustering project.

Framing reminder (see project spec): this pipeline produces hypothesis-
generation "candidates for further investigation," never classifications.
Do not print/label anything as "likely pathogenic," "reclassify," or
"confirmed" anywhere downstream of these helpers.
"""
import re
import subprocess
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser

AA3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
    'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
    'Tyr': 'Y', 'Val': 'V', 'Ter': '*',
}

PATHOGENIC_LABELS = {'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'}
BENIGN_LABELS = {'Benign', 'Likely benign', 'Benign/Likely benign'}
VUS_LABELS = {'Uncertain significance'}


def extract_fasta(rtf_or_fasta_path: str):
    """Return (header, sequence) from an RTF-wrapped or plain FASTA file."""
    path = Path(rtf_or_fasta_path)
    text = path.read_text(errors='replace')
    if text.startswith('{\\rtf'):
        result = subprocess.run(
            ['textutil', '-convert', 'txt', '-stdout', str(path)],
            capture_output=True, text=True, check=True,
        )
        text = result.stdout
    lines = text.strip().splitlines()
    header = lines[0]
    seq = ''.join(l.strip() for l in lines[1:] if l.strip())
    return header, seq


def parse_clinvar(tsv_path: str, gene_tag: str):
    """
    Parse a ClinVar TSV export.

    Applies the three fixes validated during pre-flight:
      1. Protein change is parsed from the 'Variation' column's parenthetical
         (p.XxxNNNYyy), NOT the aggregate 'Protein change' column, which mixes
         numbering across unrelated transcripts/paralogs.
      2. Rows are kept only if the Variation column's transcript resolves to
         this gene, i.e. contains '({gene_tag})' — this drops contamination
         from overlapping genes (e.g. WRAP53 in TP53, KLLN in PTEN).
      3. Classification is restricted to germline ('G:' prefix) rows only;
         somatic ('S:') and oncogenicity ('O:') rows use a different
         evidentiary framework and are excluded, not silently merged in.

    Returns a DataFrame with all rows that passed the transcript filter and
    had a parseable protein-level missense change, plus a `dropped` dict of
    counts for what was excluded and why (for transparent reporting).
    """
    rows = []
    dropped = {
        'wrong_gene_transcript': 0,
        'no_protein_change_parsed': 0,
        'non_germline_or_blank': 0,
        'stop_gained_or_synonymous': 0,
    }
    with open(tsv_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            variation = row['Variation']
            if f'({gene_tag})' not in variation:
                dropped['wrong_gene_transcript'] += 1
                continue
            m = re.search(r'p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})', variation)
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

            classification_raw = row['Classification'].strip()
            if not classification_raw.startswith('G:'):
                dropped['non_germline_or_blank'] += 1
                continue
            classification = classification_raw[2:].strip()

            if classification in PATHOGENIC_LABELS:
                bucket = 'Pathogenic'
            elif classification in BENIGN_LABELS:
                bucket = 'Benign'
            elif classification in VUS_LABELS:
                bucket = 'VUS'
            else:
                bucket = None  # e.g. 'not provided'

            rows.append({
                'variation': variation,
                'position': pos,
                'wt_aa': AA3TO1[wt3],
                'mt_aa': AA3TO1[mt3],
                'classification_raw': classification,
                'bucket': bucket,
                'variation_id': row.get('Variation ID', ''),
                'review_status': row.get('Review status', ''),
                'condition': row.get('Condition', ''),
            })
    df = pd.DataFrame(rows)
    return df, dropped


def cross_check_numbering(clinvar_df: pd.DataFrame, seq: str, struct_df: pd.DataFrame):
    """
    Check that ClinVar wt residue at each position matches (a) the UniProt
    canonical sequence and (b) the AlphaFold structure's residue at that
    position. Returns a DataFrame of any mismatching rows (empty if clean).
    """
    struct_by_pos = struct_df.set_index('position')['resname_1letter'].to_dict()
    mismatches = []
    for _, row in clinvar_df.iterrows():
        pos = row['position']
        seq_aa = seq[pos - 1] if pos - 1 < len(seq) else None
        struct_aa = struct_by_pos.get(pos)
        if seq_aa != row['wt_aa'] or struct_aa != row['wt_aa']:
            mismatches.append({
                'variation': row['variation'],
                'position': pos,
                'clinvar_wt': row['wt_aa'],
                'uniprot_seq_aa': seq_aa,
                'structure_aa': struct_aa,
            })
    return pd.DataFrame(mismatches)


def parse_structure(struct_path: str):
    """
    Parse an AlphaFold PDB/mmCIF model. Returns a DataFrame indexed by
    residue position with both backbone (Ca) and side-chain (Cb) x/y/z
    coordinates, pLDDT (stored in the B-factor field by AlphaFold DB), and
    1-letter residue name.

    Glycine has no side chain (just a hydrogen), so it has no Cb atom --
    falls back to using Ca for glycine's "cb_*" columns, flagged via
    is_glycine_cb_fallback. This is the standard convention wherever
    Cb-based distance is used in structural biology.
    """
    path = Path(struct_path)
    if path.suffix.lower() in ('.cif', '.mmcif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure('model', str(path))
    model = next(structure.get_models())
    chain = next(model.get_chains())

    rows = []
    for residue in chain:
        if 'CA' not in residue:
            continue
        ca = residue['CA']
        cb = residue['CB'] if 'CB' in residue else ca
        resname3 = residue.get_resname().capitalize()
        rows.append({
            'position': residue.id[1],
            'resname_1letter': AA3TO1.get(resname3, 'X'),
            'x': ca.coord[0], 'y': ca.coord[1], 'z': ca.coord[2],
            'cb_x': cb.coord[0], 'cb_y': cb.coord[1], 'cb_z': cb.coord[2],
            'is_glycine_cb_fallback': resname3 == 'Gly',
            'plddt': ca.get_bfactor(),
        })
    return pd.DataFrame(rows).sort_values('position').reset_index(drop=True)


def _nearest_qualifying(pos, ref_positions, ref_coords, v, min_seqdist, count_threshold_A):
    """
    Shared core: among ref_positions/ref_coords, find the CLOSEST one that
    is more than min_seqdist residues away in sequence -- not simply the
    globally closest one, which is almost always a trivial sequence
    neighbor (adjacent residues' backbones/side-chains are always close by
    construction) and would silently mask a real, more distant structural
    contact that also falls within range and would otherwise qualify.

    (Found via a Cb sensitivity check on TP53: several VUS had a valid,
    sequence-distant pathogenic residue within 6A, but the globally-nearest
    match was always a same-neighborhood residue that failed the sequence
    filter, so the real candidate was never even considered.)

    Returns None if no reference position qualifies (nothing to report),
    else a dict with dist3d, seqdist, nearest_pos, n_within_threshold
    (count of ALL ref_positions within count_threshold_A, qualifying or
    not -- local density is informative on its own, independent of the
    trivial-neighbor question).
    """
    d3 = np.linalg.norm(ref_coords - v, axis=1)
    seqd = np.abs(np.array(ref_positions) - pos)
    qualifies = seqd > min_seqdist
    n_within = int((d3 <= count_threshold_A).sum())
    if not qualifies.any():
        return None
    d3q = np.where(qualifies, d3, np.inf)
    best = int(np.argmin(d3q))
    return dict(dist3d=float(d3[best]), seqdist=float(seqd[best]),
                nearest_pos=ref_positions[best], n_within=n_within)


def compute_distances(clinvar_df: pd.DataFrame, struct_df: pd.DataFrame,
                       min_seqdist: int = 10, count_threshold_A: float = 6.0):
    """
    For every VUS residue, compute 3D distance (Angstrom) to the nearest
    QUALIFYING Pathogenic-bucket residue -- one that is also more than
    min_seqdist residues away in sequence -- and the linear sequence
    distance to that same residue. Computed independently for backbone
    (Ca) and side-chain (Cb) coordinates, since they can identify a
    different "nearest" residue (a residue's backbone can be moderately
    close while its side chain points toward or away from a given site).
    Adds columns to a copy of clinvar_df (VUS rows only get values;
    non-VUS rows get NaN).
    """
    coords_ca = struct_df.set_index('position')[['x', 'y', 'z']]
    coords_cb = struct_df.set_index('position')[['cb_x', 'cb_y', 'cb_z']]
    patho_positions = sorted(
        clinvar_df.loc[clinvar_df['bucket'] == 'Pathogenic', 'position'].unique()
    )
    patho_positions = [p for p in patho_positions if p in coords_ca.index]
    if not patho_positions:
        raise ValueError('No pathogenic residues with structure coordinates found.')
    patho_coords_ca = coords_ca.loc[patho_positions].to_numpy()
    patho_coords_cb = coords_cb.loc[patho_positions].to_numpy()

    out = clinvar_df.copy()
    for col in ['dist3d_A', 'seqdist_nearest3d', 'nearest_pathogenic_pos',
                'dist3d_cb_A', 'seqdist_cb', 'nearest_pathogenic_pos_cb']:
        out[col] = np.nan

    vus_mask = out['bucket'] == 'VUS'
    for idx, row in out[vus_mask].iterrows():
        pos = row['position']
        if pos not in coords_ca.index:
            continue
        ca_res = _nearest_qualifying(pos, patho_positions, patho_coords_ca,
                                      coords_ca.loc[pos].to_numpy(), min_seqdist, count_threshold_A)
        if ca_res:
            out.loc[idx, 'dist3d_A'] = ca_res['dist3d']
            out.loc[idx, 'seqdist_nearest3d'] = ca_res['seqdist']
            out.loc[idx, 'nearest_pathogenic_pos'] = ca_res['nearest_pos']

        cb_res = _nearest_qualifying(pos, patho_positions, patho_coords_cb,
                                      coords_cb.loc[pos].to_numpy(), min_seqdist, count_threshold_A)
        if cb_res:
            out.loc[idx, 'dist3d_cb_A'] = cb_res['dist3d']
            out.loc[idx, 'seqdist_cb'] = cb_res['seqdist']
            out.loc[idx, 'nearest_pathogenic_pos_cb'] = cb_res['nearest_pos']

    return out


def compute_distances_leave_one_out(clinvar_df: pd.DataFrame, struct_df: pd.DataFrame,
                                     count_threshold_A: float = 6.0, min_seqdist: int = 10):
    """
    Phase 2 feature builder. For EVERY row (Pathogenic, Benign, or VUS),
    compute distance to the nearest QUALIFYING Pathogenic-bucket residue
    (more than min_seqdist residues away in sequence, at a DIFFERENT
    sequence position than the row's own position), independently for
    backbone (Ca) and side-chain (Cb) coordinates.

    The "different position" exclusion matters for Pathogenic-labeled
    rows: without it, every Pathogenic row would trivially match itself
    at 0 Angstrom (or match another Pathogenic variant at the exact same
    residue), which would leak the label into the distance feature and
    make the classifier trivially, meaninglessly perfect. Excluding the
    row's own position (not just the row itself) closes that leak
    correctly even when multiple Pathogenic variants share a position.

    This exclusion is applied ONLY when the row's own bucket is
    'Pathogenic'. For VUS/Benign rows there is no self-leakage risk, and
    a separate Pathogenic-bucket variant at that same residue (a
    different amino-acid substitution) is genuine, informative signal
    (ACMG's PM5 "different pathogenic missense at the same residue"
    logic) -- it must NOT be excluded, or real evidence gets thrown away
    in favor of a much more distant, weaker match.

    The min_seqdist requirement (fixed alongside the Ca/Cb addition) stops
    a trivial sequence-adjacent pathogenic residue from masking a real,
    more distant structural contact that also falls within range -- see
    _nearest_qualifying's docstring.

    Adds: dist3d_A, seqdist_nearest3d, nearest_pathogenic_pos,
    n_pathogenic_within_threshold (Ca-based), and the Cb-based
    equivalents dist3d_cb_A, seqdist_cb, nearest_pathogenic_pos_cb,
    n_pathogenic_within_threshold_cb.
    """
    coords_ca = struct_df.set_index('position')[['x', 'y', 'z']]
    coords_cb = struct_df.set_index('position')[['cb_x', 'cb_y', 'cb_z']]
    all_patho_positions = sorted(
        clinvar_df.loc[clinvar_df['bucket'] == 'Pathogenic', 'position'].unique()
    )
    all_patho_positions = [p for p in all_patho_positions if p in coords_ca.index]
    if not all_patho_positions:
        raise ValueError('No pathogenic residues with structure coordinates found.')

    out = clinvar_df.copy()
    for col in ['dist3d_A', 'seqdist_nearest3d', 'nearest_pathogenic_pos', 'n_pathogenic_within_threshold',
                'dist3d_cb_A', 'seqdist_cb', 'nearest_pathogenic_pos_cb', 'n_pathogenic_within_threshold_cb']:
        out[col] = np.nan

    for idx, row in out.iterrows():
        pos = row['position']
        if pos not in coords_ca.index:
            continue
        if row['bucket'] == 'Pathogenic':
            ref_positions = [p for p in all_patho_positions if p != pos]
        else:
            ref_positions = all_patho_positions
        if not ref_positions:
            continue

        ref_coords_ca = coords_ca.loc[ref_positions].to_numpy()
        ca_res = _nearest_qualifying(pos, ref_positions, ref_coords_ca,
                                      coords_ca.loc[pos].to_numpy(), min_seqdist, count_threshold_A)
        if ca_res:
            out.loc[idx, 'dist3d_A'] = ca_res['dist3d']
            out.loc[idx, 'seqdist_nearest3d'] = ca_res['seqdist']
            out.loc[idx, 'nearest_pathogenic_pos'] = ca_res['nearest_pos']
            out.loc[idx, 'n_pathogenic_within_threshold'] = ca_res['n_within']

        ref_coords_cb = coords_cb.loc[ref_positions].to_numpy()
        cb_res = _nearest_qualifying(pos, ref_positions, ref_coords_cb,
                                      coords_cb.loc[pos].to_numpy(), min_seqdist, count_threshold_A)
        if cb_res:
            out.loc[idx, 'dist3d_cb_A'] = cb_res['dist3d']
            out.loc[idx, 'seqdist_cb'] = cb_res['seqdist']
            out.loc[idx, 'nearest_pathogenic_pos_cb'] = cb_res['nearest_pos']
            out.loc[idx, 'n_pathogenic_within_threshold_cb'] = cb_res['n_within']

    return out
