"""
Regenerates scripts/domains.py from scratch: the original 3 genes'
hand-curated maps stay hardcoded here (they predate the automated
pipeline), and every other gene's domain map is generated fresh from
scripts/phase2b/expansion_log.json's `domain_features` field. Run this
any time expansion_log.json changes -- domains.py never needs manual
editing again.
"""
import json
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent  # project root

LOG_PATH = (PROJECT_ROOT / "scripts/phase2b/expansion_log.json")
DOMAINS_PATH = (PROJECT_ROOT / "scripts/domains.py")

# Hand-curated for the original 3 genes (predates the automated pipeline).
MANUAL_DOMAINS = {
    "TP53": [
        (102, 292, "DNA-binding domain"),
        (325, 356, "Oligomerization (tetramerization) domain"),
        (1, 44, "Transactivation domain (acidic)"),
        (50, 96, "Disordered (TAD2/proline-rich linker)"),
        (282, 325, "Disordered (DBD-tetramerization linker)"),
        (351, 393, "Disordered (C-terminal regulatory)"),
    ],
    "BRCA1": [
        (24, 65, "RING-type zinc finger"),
        (1642, 1736, "BRCT domain 1"),
        (1756, 1855, "BRCT domain 2"),
        (230, 270, "Disordered"),
        (306, 338, "Disordered"),
        (534, 570, "Disordered"),
        (654, 709, "Disordered"),
        (1181, 1216, "Disordered"),
        (1322, 1387, "Disordered"),
        (1440, 1505, "Disordered"),
        (1565, 1596, "Disordered"),
    ],
    "PTEN": [
        (14, 185, "Phosphatase (tensin-type) domain"),
        (190, 350, "C2 (tensin-type) domain"),
        (352, 403, "Disordered (C-terminal tail)"),
    ],
}

HEADER = '''"""
Curated domain/region annotations. TP53/BRCA1/PTEN are hand-curated
(predate the automated pipeline); every other gene is generated
automatically from UniProt features by scripts/phase2b/generate_domains.py
-- do not hand-edit those blocks, they get overwritten on the next run.
Ranges are UniProt canonical numbering (1-based, inclusive). Order
matters: first matching range wins (most specific/structurally
meaningful first).
"""
'''

def emit_gene_block(gene, ranges):
    lines = [f"{gene}_DOMAINS = ["]
    for start, end, label in ranges:
        label_escaped = label.replace('"', '\\"')
        lines.append(f'    ({start}, {end}, "{label_escaped}"),')
    lines.append("]")
    return "\n".join(lines)

log = json.loads(LOG_PATH.read_text()) if LOG_PATH.exists() else {}
auto_genes = {g: e["domain_features"] for g, e in log.items()
              if e.get("structure_passed") and e.get("domain_features")}

all_genes = list(MANUAL_DOMAINS.keys()) + sorted(auto_genes.keys())

blocks = []
for gene in MANUAL_DOMAINS:
    blocks.append(emit_gene_block(gene, MANUAL_DOMAINS[gene]))
for gene in sorted(auto_genes):
    blocks.append(emit_gene_block(gene, auto_genes[gene]))

map_entries = "\n".join(f'    "{g}": {g}_DOMAINS,' for g in all_genes)

out = HEADER + "\n" + "\n\n".join(blocks) + f"""

DOMAIN_MAPS = {{
{map_entries}
}}


def annotate_domain(gene: str, position: int) -> str:
    for start, end, label in DOMAIN_MAPS[gene]:
        if start <= position <= end:
            return label
    return "Unannotated / linker"
"""

DOMAINS_PATH.write_text(out)
print(f"Regenerated {DOMAINS_PATH} with {len(all_genes)} genes "
      f"({len(MANUAL_DOMAINS)} manual + {len(auto_genes)} from log).")
