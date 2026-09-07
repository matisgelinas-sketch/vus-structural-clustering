"""
Single entry point for adding genes to the project. Consolidates what
used to be a two-stage, two-fork manual process (prescreen script, then
a separate fetch+validate+save script, with domains.py hand-edited
afterward) into one resumable pass.

Usage:
  python3 scripts/gene_pipeline.py                          # process scripts/candidate_genes.json
  python3 scripts/gene_pipeline.py GENE1 GENE2 ...           # process specific genes
  NCBI_API_KEY=xxxx python3 scripts/gene_pipeline.py         # ~3x faster (10 req/s vs 3 req/s)

Safe to re-run or interrupt: genes already in expansion_log.json are
skipped instantly, and the log is written after every single gene, not
just at the end.

After this finishes, domains.py is regenerated automatically. The
downstream scripts (run_new_genes_phase1.py, build_phase2_dataset.py)
auto-discover genes from expansion_log.json -- no manual edits needed
anywhere to start using newly-added genes.
"""
import json
import subprocess
import sys
import time
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # project root

sys.path.insert(0, str(Path(__file__).parent / "phase2b"))
sys.path.insert(0, str(Path(__file__).parent))
from ncbi_clinvar import process_gene

GENES_BASE = (PROJECT_ROOT / "Genes")
LOG_PATH = (PROJECT_ROOT / "scripts/phase2b/expansion_log.json")
CANDIDATES_PATH = (PROJECT_ROOT / "scripts/candidate_genes.json")

if len(sys.argv) > 1:
    candidates = sys.argv[1:]
else:
    if not CANDIDATES_PATH.exists():
        print(f"No genes given and {CANDIDATES_PATH} doesn't exist. "
              f"Pass gene symbols as arguments, or create that file with a JSON list.")
        sys.exit(1)
    candidates = json.loads(CANDIDATES_PATH.read_text())

log = json.loads(LOG_PATH.read_text()) if LOG_PATH.exists() else {}

already_done = [g for g in candidates if g in log]
todo = [g for g in candidates if g not in log]
if already_done:
    print(f"Skipping {len(already_done)} already-processed genes: {', '.join(already_done)}")
print(f"Processing {len(todo)} genes...\n")

t0 = time.time()
for i, gene in enumerate(todo, 1):
    print(f"[{i}/{len(todo)}] {gene} ...", end=" ", flush=True)
    entry = process_gene(gene, log, GENES_BASE, LOG_PATH)
    if entry.get("structure_passed"):
        print(f"PASSED (P={entry['n_pathogenic']} B={entry['n_benign']} V={entry['n_vus']})")
    elif entry.get("bar_passed") is False:
        print(f"failed data bar (P={entry.get('n_pathogenic','?')} B={entry.get('n_benign','?')})")
    else:
        print(f"excluded -- {entry.get('reason')}")

elapsed = time.time() - t0
print(f"\nDone in {elapsed/60:.1f} min.")

passed = [g for g, e in log.items() if e.get("structure_passed")]
print(f"\n{len(passed)} genes fully validated and ready to train on: {', '.join(sorted(passed))}")

print("\nRegenerating domains.py from the log...")
subprocess.run([sys.executable, str(Path(__file__).parent / "phase2b" / "generate_domains.py")], check=True)

print("\nNext: run_new_genes_phase1.py, then build_phase2_dataset.py / train_phase2_model.py / "
      "score_vus_phase2.py -- all auto-discover genes from expansion_log.json now, no edits needed.")
