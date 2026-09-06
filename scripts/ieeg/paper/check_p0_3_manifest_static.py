"""Static sanity check of p0_3_ablation/run_manifest.json structure.

Pre-flight for the validator's addendum-5.9 cross-check: verifies manifest
filenames sit under the canonical root, per-cell counts/seed sets are complete,
and entry cell_id matches the path cell. Read-only; not part of the pipeline.
"""
import collections
import json
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
MANIFEST = os.path.join(ROOT, "processed", "subject", "group",
                        "ieeg_methods_paper", "p0_3_ablation", "run_manifest.json")

m = json.load(open(MANIFEST))
prefix = os.path.join(os.path.abspath(ROOT), "processed", "subject", "group",
                      "ieeg_methods_paper", "p0_3_ablation") + os.sep
fns = [e["filename"] for e in m["checkpoint_sha256"]]
print("entries:", len(fns), "unique:", len(set(fns)))
print("all under canonical root:", all(f.startswith(prefix) for f in fns))
counts = collections.Counter()
seeds = collections.defaultdict(list)
okcell = True
for e in m["checkpoint_sha256"]:
    rel = e["filename"][len(prefix):]
    arm, checkpoints, cell, fn = rel.split(os.sep)
    seed = int(fn.replace("world_", "").replace(".mat", ""))
    counts[cell] += 1
    seeds[cell].append(seed)
    okcell = okcell and e["cell_id"] == cell
print("per-cell counts:", sorted(counts.items()))
print("entry cell_id == path cell:", okcell)
print("seed sets 1..64 per cell:", all(sorted(v) == list(range(1, 65)) for v in seeds.values()))
print("worlds:", m["worlds"], "status:", m["status"], "workers:", m["workers"],
      "contract:", m["contract_sha256"][:16])
missing = [f for f in fns if not os.path.isfile(f)]
print("files missing on disk:", len(missing))
ok = (len(fns) == 768 and len(set(fns)) == 768
      and all(f.startswith(prefix) for f in fns)
      and sorted(counts.values()) == [64] * 12
      and all(sorted(v) == list(range(1, 65)) for v in seeds.values())
      and okcell and not missing
      and m["worlds"] == 768 and m["status"] == "complete")
sys.exit(0 if ok else 1)
