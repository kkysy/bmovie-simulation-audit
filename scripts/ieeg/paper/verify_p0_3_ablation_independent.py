#!/usr/bin/env python3
"""Independent inference-layer re-verification of all 768 P0-3 checkpoints.

Second implementation, deliberately sharing no inference code with the
production runners or the MATLAB validator:

  * Identity layer (all 768, pure h5py): numeric identity only -- seed_index
    against filename, world_seed against the contract seed formulas, A_uV
    against the contract (Mpower). MATLAB string fields (status/arm/cell_id/
    SHA) are opaque heap objects to h5py; string identity is covered by the
    MATLAB validator (5830 assertions) and by the bridge's per-file
    cell/seed assertion, cross-checked here against the path.
  * Mphase (128 worlds, pure h5py): subject_mphase_vector is a plain
    float64 dataset; the full exact sign-flip family (two-sided p and
    studentized max-T p) is recomputed in numpy and compared against the
    stored mphase_shape fields.
  * Mpower (640 worlds): subject vectors arrive via the zero-logic bridge
    dump_p0_3_subject_vectors.m (MATLAB table serialization is not
    h5py-readable); the exact sign-flip family is recomputed in numpy from
    those vectors and compared against the stored signflip reduced fields,
    which are read directly from the checkpoint via h5py (not through the
    bridge).

p values are integer counts over 2^n; the reference comparison is the
count difference. count diff 0 = exact match. diff 1 is only accepted when
the boundary statistic lies within 1e-9 (relative) of |observed| (a genuine
cross-language ulp tie); anything larger fails. Verdict rows and a summary
are written next to the aggregated tables.
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("E:/HighEV_sampling/Bmovie")
BASE = ROOT / "processed/subject/group/ieeg_methods_paper/p0_3_ablation"
CONTRACT = ROOT / "scripts/ieeg/paper/acc00_sim_p0_3_ablation_contract.json"
BRIDGE_TSV = BASE / "tables/p0_3_subject_vectors_bridge.tsv"
REPORT_TSV = BASE / "tables/p0_3_independent_verification.tsv"
ULP_TIE_TOL = 1e-9  # relative to max(1, |observed|)
# Production studentization uses MATLAB std(stats,0,1): w=0 normalizes by
# N-1 (sample SD). numpy therefore uses ddof=1 below. (2026-09-03: an
# earlier revision of this verifier used ddof=0 and reported a spurious
# "MATLAB builtin std() numerical deviation" of exactly sqrt(N/(N-1))-1;
# GPT root-caused it to this definitional mismatch -- there was no std()
# arithmetic defect.)
SD_RTOL = 1e-9

failures: list[str] = []
rows: list[dict] = []
sd_max_rel = 0.0


def check(cond: bool, msg: str) -> None:
    if not cond:
        failures.append(msg)
        print("FAIL:", msg)


def signflip_family(x: np.ndarray):
    """Exact exhaustive sign-flip family, second implementation (numpy)."""
    n = x.shape[0]
    bits = (np.arange(2 ** n, dtype=np.uint64)[:, None] >> np.arange(n - 1, -1, -1, dtype=np.uint64)) & 1
    signs = 1 - 2 * bits.astype(np.float64)  # row 0 = all-plus (observed)
    stats = signs @ x / n
    obs = stats[0]
    sd = stats.std(axis=0, ddof=1)  # sample SD over the enumeration (MATLAB w=0)
    scale = np.where(sd == 0, 1.0, sd)
    stud = stats / scale
    stud_obs = stud[0]
    maxstats = np.abs(stud).max(axis=1)
    maxobs = float(np.abs(stud_obs).max())
    return obs, stats, sd, stud_obs, maxstats, maxobs


def count_diff(tag: str, stored_p: float, null_stats: np.ndarray, obs_abs: float) -> None:
    """Exact-count comparison; a +-1 shift is accepted only with a genuine
    near-boundary statistic (cross-language ulp tie), excluding the observed
    row's exact self-tie at index 0."""
    n_enum = null_stats.shape[0]
    count_exact = int(np.count_nonzero(null_stats >= obs_abs))
    diff = count_exact - round(stored_p * n_enum)
    if diff != 0:
        dist = np.abs(null_stats - obs_abs) / max(1.0, obs_abs)
        # exact self-ties (the observed row and its sign-mirror) sit at dist 0
        # in both implementations; the diagnostic needs the nearest strictly
        # positive distance.
        positive = dist[dist > 0]
        near = float(np.min(positive)) if positive.size else np.inf
        check(abs(diff) == 1 and near < ULP_TIE_TOL,
              f"{tag}: p count diff {diff} (stored {stored_p}, recomputed {count_exact / n_enum}, nearest boundary reldist {near:.2e})")
    rows[-1]["metrics"] += f" [{tag.split('/')[-1]} p={stored_p:.6g} dcount={diff}]"


def sd_family_check(tag: str, name: str, stored: np.ndarray, recomputed: np.ndarray) -> None:
    global sd_max_rel
    s = np.asarray(stored, dtype=np.float64).ravel()
    r = np.asarray(recomputed, dtype=np.float64).ravel()
    denom = np.maximum(np.abs(r), 1e-300)
    rel = float(np.max(np.where(s == 0, np.abs(r) / denom, np.abs(s - r) / denom)))
    sd_max_rel = max(sd_max_rel, rel)
    check(rel <= SD_RTOL, f"{tag}: {name} mismatch stored={s} recomputed={r} (rel {rel:.2e})")


def verify_signflip(tag: str, x: np.ndarray, stored: dict) -> None:
    obs, stats, sd, stud_obs, maxstats, maxobs = signflip_family(x)
    check(np.allclose(obs, stored["observed"], rtol=0, atol=1e-12 * max(1.0, float(np.max(np.abs(obs))))),
          f"{tag}: observed mismatch stored={stored['observed']} recomputed={obs}")
    sd_family_check(tag, "studentization_sd", stored["studentization_sd"], sd)
    sd_family_check(tag, "studentized_observed", stored["studentized_observed"], stud_obs)
    sd_family_check(tag, "max_observed", stored["max_observed"], maxobs)
    for j, name in enumerate(stored["metric_names"]):
        col = x[:, j]
        if not np.all(np.isfinite(col)):
            check(np.isnan(stored["p_two_sided"][j]), f"{tag}/{name}: non-finite column but stored p not NaN")
            continue
        count_diff(f"{tag}/{name}", float(stored["p_two_sided"][j]), np.abs(stats[:, j]), abs(float(obs[j])))
    # studentized max-T p (joint over metrics)
    count_diff(f"{tag}/maxT", float(stored["p_maxstat"]), maxstats, maxobs)


def main() -> int:
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    a_expected = float(c["additive"]["A_uV"])
    seed_p = c["seed"]["P01_to_P10"]
    seed_m = c["seed"]["M01_to_M02"]
    n_subj = dict(zip(c["sign_flip"]["runs"], c["sign_flip"]["subject_counts"]))

    # bridge CSV: subject vectors for all 640 Mpower worlds
    bridge: dict[tuple[str, int, str], list[tuple[str, float, float, float]]] = {}
    with BRIDGE_TSV.open(encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            key = (r["cell_id"], int(r["seed_index"]), r["run"])
            bridge.setdefault(key, []).append((r["subject"], float(r["Mpower"]), float(r["Mphase"]), float(r["n_pairs"])))
    check(len(bridge) == 640 * 2, f"bridge cardinality {len(bridge)} != 1280 world-runs")

    cells = [(cid, "mpower", "Mpower", int(w)) for cid, w in
             ((x["id"], x["worlds"]) for x in c["cells"]["Mpower"])]
    cells += [(x["id"], "mphase", "Mphase", int(x["worlds"])) for x in c["cells"]["Mphase"]]

    n_files = 0
    for cell_id, arm_dir, arm, n_worlds in cells:
        for si in range(1, n_worlds + 1):
            path = BASE / arm_dir / "checkpoints" / cell_id / f"world_{si:03d}.mat"
            if not path.is_file():
                check(False, f"missing checkpoint {path}")
                continue
            n_files += 1
            tag = f"{cell_id}/{si}"
            rows.append({"checkpoint": tag, "arm": arm, "metrics": ""})
            with h5py.File(path, "r") as h:
                cp = h["checkpoint"]
                check(int(cp["seed_index"][()][0, 0]) == si, f"{tag}: seed_index != filename")
                ws = float(cp["world_seed"][()][0, 0])
                if arm == "Mpower":
                    exp = float(seed_p["paper_base_seed"]) + 100000 * seed_p["scenario_index"] + 1000 * seed_p["effect_index"] + si
                    check(ws == exp, f"{tag}: world_seed {ws} != {exp}")
                    check(float(cp["injection_provenance/A_uV"][()][0, 0]) == a_expected,
                          f"{tag}: A_uV != contract {a_expected}")
                    for run in c["sign_flip"]["runs"]:
                        b = cp[f"signflip/{run}"]
                        stored = {k: b[k][()].ravel().astype(np.float64) for k in
                                  ("observed", "p_two_sided", "studentization_sd", "studentized_observed")}
                        stored["max_observed"] = float(b["max_observed"][()][0, 0])
                        stored["p_maxstat"] = float(b["p_maxstat"][()][0, 0])
                        stored["metric_names"] = c["sign_flip"]["metric_order"]
                        vec = bridge.get((cell_id, si, run))
                        if vec is None:
                            check(False, f"{tag}/{run}: missing bridge rows")
                            continue
                        vec = sorted(vec, key=lambda t: t[0])  # contract: lexicographic subject order
                        check(len(vec) == n_subj[run], f"{tag}/{run}: {len(vec)} subjects != contract {n_subj[run]}")
                        x = np.array([[m, p] for _, m, p, _ in vec], dtype=np.float64)
                        verify_signflip(f"{tag}/{run}", x, stored)
                else:
                    exp = float(seed_m["base_seed"]) + 100000 * seed_m["scenario_index"] + si
                    check(ws == exp, f"{tag}: world_seed {ws} != {exp}")
                    for run in c["sign_flip"]["runs"]:
                        b = cp[f"mphase_shape/{run}"]
                        x = b["subject_mphase_vector"][()].ravel().astype(np.float64)
                        keep = np.isfinite(x)
                        x = x[keep]
                        if x.size == 0:
                            check(np.isnan(float(b["p_two_sided"][()][0, 0])), f"{tag}/{run}: empty run but p not NaN")
                            continue
                        check(x.size == n_subj[run], f"{tag}/{run}: {x.size} subjects != contract {n_subj[run]}")
                        stored = {
                            "observed": np.array([float(b["observed"][()][0, 0])]),
                            "p_two_sided": np.array([float(b["p_two_sided"][()][0, 0])]),
                            "studentization_sd": np.array([float(b["studentization_sd"][()][0, 0])]),
                            "studentized_observed": np.array([float(b["observed_standardized_score"][()][0, 0])]),
                            "max_observed": abs(float(b["observed_standardized_score"][()][0, 0])),
                            "p_maxstat": float(b["p_maxstat"][()][0, 0]),
                            "metric_names": ["Mphase"],
                        }
                        verify_signflip(f"{tag}/{run}", x.reshape(-1, 1), stored)

    check(n_files == 768, f"checkpoint count {n_files} != 768")
    with REPORT_TSV.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, ["checkpoint", "arm", "metrics"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"verified {n_files} checkpoints; failures: {len(failures)}; "
          f"max sd-family rel deviation: {sd_max_rel:.3e}")
    print("INDEPENDENT_VERIFICATION", "FAIL" if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
