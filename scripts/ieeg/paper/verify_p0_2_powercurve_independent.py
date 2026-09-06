#!/usr/bin/env python
"""Independent post-run verification of the P0-2 986-world power curve.

ZCode-written second implementation. Does NOT import or call
summarize_p0_2_powercurve.py. Reads the frozen contract and the raw
checkpoints with h5py only (no MATLAB), then independently derives:

  1. Per-world exact Mpower two-sided p from the stored sign-flip
     enumeration, and checks it equals the stored p_two_sided (metric 0).
  2. Internal consistency of every stored flip statistic with the stored
     sign matrix: solves x from stats = (signs' x)/n by least squares and
     requires an exact-residual solution whose mean equals the stored
     observed slope (this proves the enumeration is a genuine sign-flip
     orbit of some subject vector with all-plus = observed).
  3. Per-cell power = #{p <= 0.05}/n with Clopper-Pearson 95% and MCSE.
  4. Theta-inversion coverage, recomputed with a separately written
     implementation of the frozen translation identity
        T_j(theta) = S_j - theta*mean_i(sign_j,i),  T_obs(theta) = obs - theta
     on all R2 worlds and a fixed R1 subsample (R1_COVERAGE_PER_CELL per
     cell), then compared against the aggregator's per-world table if
     present.

Exit 0 only if every hard check passes. Writes
  tables/independent_power_summary.tsv
  tables/independent_coverage_check.tsv
under the p0_2_sim_powercurve directory.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.special import betaincinv

ROOT = Path(__file__).resolve().parents[3]
PAPER = ROOT / "scripts/ieeg/paper/acc00_sim_powercurve_contract.json"
MAIN = ROOT / "processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/profiling/main"
TABLES = ROOT / "processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables"
AGG_WORLDS = TABLES / "p0_2_powercurve_worlds.tsv"

ALPHA = 0.05
R1_COVERAGE_PER_CELL = 8  # deterministic: seed_index 1..8 of each cell

failures: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        failures.append(msg)
        print(f"FAIL: {msg}")


def cp95(k: int, n: int) -> tuple[float, float]:
    lo = 0.0 if k == 0 else float(betaincinv(k, n - k + 1, 0.025))
    hi = 1.0 if k == n else float(betaincinv(k + 1, n - k, 0.975))
    return lo, hi


def main() -> int:
    contract = json.loads(PAPER.read_text(encoding="utf-8"))
    grid_ids = contract["additive"]["grid_id"]
    grid_worlds = contract["additive"]["grid_worlds"]
    a_grid = contract["additive"]["a_erp_grid_uV"]
    gain = float(contract["bridge"]["expected_slope_gain_G_log10_per_uV2_per_z"])
    runs = contract["sign_flip"]["runs"]
    n_subj = dict(zip(runs, contract["sign_flip"]["subject_counts"], strict=True))
    inv = contract["theta_inversion"]
    theta = float(inv["lower"]) + np.arange(round((inv["upper"] - inv["lower"]) / inv["step"]) + 1) * float(inv["step"])
    hit_tol = float(inv["hit_tolerance"])

    rows: list[dict] = []
    coverage_rows: list[dict] = []
    n_files = 0
    for gi, (gid, n_worlds) in enumerate(zip(grid_ids, grid_worlds, strict=True), start=1):
        expected_slope = float(a_grid[gi]) ** 2 * gain
        for si in range(1, n_worlds + 1):
            path = MAIN / f"world_g{gi:02d}_n{si:03d}.mat"
            if not path.is_file():
                check(False, f"missing checkpoint {path.name}")
                continue
            n_files += 1
            with h5py.File(path, "r") as h:
                check(abs(float(h["checkpoint/truth_bridge_slope"][()][0, 0]) - expected_slope) < 1e-18,
                      f"{path.name}: truth_bridge_slope != A^2*gain from contract")
                check(float(h["checkpoint/scenario_index"][()][0, 0]) == 2.0,
                      f"{path.name}: scenario_index != 2")
                for run in runs:
                    base = f"checkpoint/inference/{run}"
                    signs = np.asarray(h[f"{base}/signs"][()], dtype=np.float64)
                    stats = np.asarray(h[f"{base}/statistics"][()], dtype=np.float64)
                    observed = np.asarray(h[f"{base}/observed"][()], dtype=np.float64).reshape(-1, order="F")
                    p_two = np.asarray(h[f"{base}/p_two_sided"][()], dtype=np.float64).reshape(-1, order="F")
                    ns = n_subj[run]
                    n_enum = 2 ** ns
                    check(signs.shape == (ns, n_enum), f"{path.name}/{run}: signs shape {signs.shape}")
                    check(np.all(signs[:, 0] == 1.0), f"{path.name}/{run}: first flip not all-plus")
                    check(np.array_equal(stats[:, 0], observed),
                          f"{path.name}/{run}: all-plus statistic != observed")
                    # 1. exact p recompute (integer count, no tolerance)
                    obs_m = float(observed[0])
                    count = int(np.count_nonzero(np.abs(stats[0]) >= abs(obs_m)))
                    check(count / n_enum == float(p_two[0]),
                          f"{path.name}/{run}: stored Mpower p {p_two[0]} != recomputed {count}/{n_enum}")
                    # 2. sign-orbit consistency: stats = signs' x / n exactly
                    x_hat, *_ = np.linalg.lstsq(signs.T / ns, stats[0], rcond=None)
                    resid = np.max(np.abs(signs.T @ x_hat / ns - stats[0]))
                    scale = max(1.0, float(np.max(np.abs(stats[0]))))
                    check(resid / scale < 1e-12,
                          f"{path.name}/{run}: flip statistics not an exact sign orbit (rel resid {resid/scale:.2e})")
                    check(abs(float(np.mean(x_hat)) - obs_m) <= 1e-12 * scale,
                          f"{path.name}/{run}: recovered subject mean != observed slope")
                    rows.append({"grid_id": gid, "run": run, "seed_index": si,
                                 "true_bridge_slope": expected_slope, "observed_slope": obs_m,
                                 "p_mpower": count / n_enum, "reject": int(count / n_enum <= ALPHA)})
                    # 4. coverage recompute on the designated subset
                    do_cov = (run == "R2") or (run == "R1" and si <= R1_COVERAGE_PER_CELL)
                    if do_cov:
                        means = signs.mean(axis=0)
                        p_curve = np.empty(theta.size)
                        for b in range(0, theta.size, 128):
                            th = theta[b:b + 128]
                            translated = stats[0][:, None] - means[:, None] * th[None, :]
                            p_curve[b:b + th.size] = np.count_nonzero(
                                np.abs(translated) >= np.abs(obs_m - th)[None, :], axis=0) / n_enum
                        accept = p_curve > ALPHA
                        hits = np.flatnonzero(accept)
                        nearest = int(np.argmin(np.abs(theta - expected_slope)))
                        coverage_rows.append({
                            "grid_id": gid, "run": run, "seed_index": si,
                            "ci_low": float(theta[hits[0]]) if hits.size else math.nan,
                            "ci_high": float(theta[hits[-1]]) if hits.size else math.nan,
                            "truth_in_ci": bool(abs(float(theta[nearest]) - expected_slope) <= hit_tol
                                                and accept[nearest]),
                        })
            if n_files % 100 == 0:
                print(f"  ... {n_files}/986 checkpoints verified", flush=True)

    print(f"checkpoints verified: {n_files}/986")

    # 3. per-cell power summary
    power_rows: list[dict] = []
    for gid in grid_ids:
        for run in runs:
            sub = [r for r in rows if r["grid_id"] == gid and r["run"] == run]
            n = len(sub)
            k = sum(r["reject"] for r in sub)
            lo, hi = cp95(k, n)
            truth = sub[0]["true_bridge_slope"]
            obs = np.array([r["observed_slope"] for r in sub])
            power_rows.append({
                "grid_id": gid, "run": run, "n_worlds": n, "rejections": k,
                "power_hat": k / n, "cp95_low": lo, "cp95_high": hi,
                "mcse": math.sqrt((k / n) * (1 - k / n) / n),
                "true_bridge_slope": truth,
                "observed_slope_mean": float(obs.mean()),
                "bias": float(obs.mean()) - truth,
                "rmse": float(np.sqrt(np.mean((obs - truth) ** 2))),
            })
            print(f"{gid} {run}: power {k}/{n} = {k/n:.4f}  CP95 [{lo:.4f}, {hi:.4f}]"
                  f"  truth {truth:.3e}  bias {obs.mean()-truth:+.3e}")

    # coverage summary + optional comparison against aggregator per-world table
    cov_summary: dict[tuple[str, str], list[dict]] = {}
    for r in coverage_rows:
        cov_summary.setdefault((r["grid_id"], r["run"]), []).append(r)
    for (gid, run), sub in sorted(cov_summary.items()):
        k = sum(1 for r in sub if r["truth_in_ci"])
        lo, hi = cp95(k, len(sub))
        print(f"coverage {gid} {run}: {k}/{len(sub)} = {k/len(sub):.4f}  CP95 [{lo:.4f}, {hi:.4f}]")

    if AGG_WORLDS.is_file():
        agg = {(r["grid_id"], r["run"], int(r["seed_index"])): r
               for r in csv.DictReader(AGG_WORLDS.open(encoding="utf-8"), delimiter="\t")}
        n_cmp = 0
        for r in rows:
            key = (r["grid_id"], r["run"], r["seed_index"])
            a = agg.get(key)
            if a is None:
                check(False, f"aggregator world table missing {key}")
                continue
            n_cmp += 1
            check(int(a["reject_alpha_0_05"]) == r["reject"],
                  f"{key}: aggregator reject {a['reject_alpha_0_05']} != independent {r['reject']}")
            check(abs(float(a["p_two_sided_mpower"]) - r["p_mpower"]) == 0.0,
                  f"{key}: aggregator p {a['p_two_sided_mpower']} != independent {r['p_mpower']}")
        for r in coverage_rows:
            a = agg.get((r["grid_id"], r["run"], r["seed_index"]))
            if a is None:
                continue
            check(bool(int(a["truth_in_ci"])) == r["truth_in_ci"],
                  f"{r['grid_id']}/{r['run']}/n{r['seed_index']:03d}: aggregator truth_in_ci "
                  f"{a['truth_in_ci']} != independent {r['truth_in_ci']}")
            if not math.isnan(r["ci_low"]):
                check(abs(float(a["ci_theta_low"]) - r["ci_low"]) < 1e-9
                      and abs(float(a["ci_theta_high"]) - r["ci_high"]) < 1e-9,
                      f"{r['grid_id']}/{r['run']}/n{r['seed_index']:03d}: CI bounds differ from aggregator")
        print(f"aggregator cross-check: {n_cmp} world-run rows compared")
    else:
        print(f"NOTE: aggregator table {AGG_WORLDS.name} not present yet; skipped cross-check")

    TABLES.mkdir(parents=True, exist_ok=True)
    with (TABLES / "independent_power_summary.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(power_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(power_rows)
    with (TABLES / "independent_coverage_check.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(coverage_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(coverage_rows)

    if failures:
        print(f"\nINDEPENDENT VERIFICATION FAILED: {len(failures)} problem(s)")
        return 1
    print("\nINDEPENDENT VERIFICATION PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
