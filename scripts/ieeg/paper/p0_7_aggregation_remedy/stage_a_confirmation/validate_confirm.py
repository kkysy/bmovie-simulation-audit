# Release port of the analysis-workspace tool of the same name; see README in this folder.
"""Independent validator for Stage A bridge exports and main summaries."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, skew

def _repo_root():
    for _q in Path(__file__).resolve().parents:
        if (_q / "MANIFEST.md").is_file():
            return _q
    raise RuntimeError("release repository root not found")
ROOT = _repo_root()
VERIFY_DIR = ROOT / "scripts" / "ieeg" / "paper"
sys.path.insert(0, str(VERIFY_DIR))
from verify_p0_6_paired_contrast import newcombe_method10  # noqa: E402


def cp_interval(k: int, n: int) -> tuple[float, float]:
    return (0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1)),
            1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k)))


def exact_signflip(x: np.ndarray) -> tuple[float, float, float]:
    n = len(x)
    if n == 0:
        return math.nan, math.nan, math.nan
    rows = np.arange(2 ** n, dtype=np.uint32)[:, None]
    bits = (rows >> np.arange(n - 1, -1, -1, dtype=np.uint32)) & 1
    signs = 1.0 - 2.0 * bits
    signs[0, :] = 1.0
    stats = signs @ x / n
    observed = float(stats[0])
    p = float(np.mean(np.abs(stats) >= abs(observed)))
    sd = float(np.std(stats, ddof=1))
    return observed, p, sd


def close(a: float, b: float, tol: float = 1e-12) -> bool:
    return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= tol


def validate(output: Path, mode: str) -> dict:
    pairs = pd.read_csv(output / "validator_pair_rows.tsv", sep="\t")
    subjects = pd.read_csv(output / "validator_subject_rows.tsv", sep="\t")
    worlds = pd.read_csv(output / "validator_world_rows.tsv", sep="\t")
    errors: list[str] = []
    # Per (reducer, run) one-flip discrepancy stats: accumulated |recomputed p - checkpoint p|
    # and the number of worlds whose two p values straddle each reporting threshold.
    flip_stats: dict[tuple[str, str], dict] = {}
    expected = 120 if mode == "formal" else None
    seeds = sorted(pairs.seed_index.unique())
    if expected is not None and seeds != list(range(1, expected + 1)):
        errors.append("formal seed_index set is not exactly 1:120")
    if pairs.world_seed.min() < 20860828 or pairs.world_seed.max() > 20860947:
        errors.append("world seed outside frozen independent range")
    if not (pairs.mphase_surrogate_reducer == "mean").all():
        errors.append("production surrogate reducer is not uniformly mean")
    if not (pairs.refit_surrogate_count == 32).all():
        errors.append("refit surrogate count differs from 32")

    recomputed_worlds = []
    for seed in seeds:
        ps = pairs[pairs.seed_index == seed]
        for run in ("R1", "R2"):
            pr = ps[ps.run == run]
            for reducer in ("median", "mean"):
                values = []
                for subject, g in pr.groupby("subject", sort=True):
                    v = g.corrected_mphase.to_numpy(dtype=float)
                    v = v[np.isfinite(v)]
                    if not len(v):
                        continue
                    stat = float(np.median(v) if reducer == "median" else np.mean(v))
                    pair_skew = float(skew(v, bias=False)) if len(v) >= 3 and not np.all(v == v[0]) else math.nan
                    sr = subjects[(subjects.seed_index == seed) & (subjects.run == run) &
                                  (subjects.reducer == reducer) & (subjects.subject == subject)]
                    if len(sr) != 1:
                        errors.append(f"seed {seed} {run} {reducer} {subject}: subject row cardinality")
                    else:
                        row = sr.iloc[0]
                        if not close(stat, float(row.Mphase)):
                            errors.append(f"seed {seed} {run} {reducer} {subject}: reducer mismatch")
                        if int(row.n_pairs_finite) != len(v) or not close(pair_skew, float(row.pair_corrected_skewness)):
                            errors.append(f"seed {seed} {run} {reducer} {subject}: pair diagnostic mismatch")
                    values.append(stat)
                x = np.asarray(values, dtype=float)
                observed, p, sd = exact_signflip(x)
                wr = worlds[(worlds.seed_index == seed) & (worlds.run == run) & (worlds.reducer == reducer)]
                if len(wr) != 1:
                    errors.append(f"seed {seed} {run} {reducer}: world row cardinality")
                else:
                    row = wr.iloc[0]
                    # Cross-language exact sign-flip can differ by one flip (1/2^n) when a
                    # sign vector ties the observed statistic in floating point.
                    p_tol = 1.0 / float(2 ** len(x)) + 1e-12
                    pc = float(row.p_two_sided)
                    fs = flip_stats.setdefault((reducer, run), {"budget": 0.0, "s005": 0, "s050": 0, "s075": 0})
                    fs["budget"] += abs(p - pc)
                    fs["s005"] += int((p <= 0.05) != (pc <= 0.05))
                    fs["s050"] += int((p < 0.5) != (pc < 0.5))
                    fs["s075"] += int((p > 0.75) != (pc > 0.75))
                    if (not close(observed, float(row.observed))
                            or abs(p - pc) > p_tol
                            or not close(sd, float(row.studentization_sd))):
                        errors.append(f"seed {seed} {run} {reducer}: exact sign-flip mismatch")
                recomputed_worlds.append({"seed_index": seed, "run": run, "reducer": reducer,
                                          "observed": observed, "p_two_sided": p, "n_subjects": len(x)})

    rw = pd.DataFrame(recomputed_worlds)
    fpr_rows, dist_rows = [], []
    for reducer in ("median", "mean"):
        for run in ("R1", "R2"):
            g = rw[(rw.reducer == reducer) & (rw.run == run)]
            reject = g.p_two_sided <= 0.05
            k, n = int(reject.sum()), len(g)
            lo, hi = cp_interval(k, n) if n else (math.nan, math.nan)
            fpr_rows.append({"reducer": reducer, "run": run, "n_worlds": n, "rejections": k,
                             "fpr": k / n if n else math.nan, "cp_low": lo, "cp_high": hi})
            below, above = int((g.p_two_sided < 0.5).sum()), int((g.p_two_sided > 0.75).sum())
            bci, aci = cp_interval(below, n) if n else (math.nan, math.nan), cp_interval(above, n) if n else (math.nan, math.nan)
            sd = float(g.observed.std(ddof=1)) if n > 1 else math.nan
            dist_rows.append({"reducer": reducer, "run": run, "n_worlds": n, "mean_p": float(g.p_two_sided.mean()) if n else math.nan,
                              "p_lt_0p5_count": below, "p_lt_0p5_fraction": below / n if n else math.nan,
                              "p_lt_0p5_cp_low": bci[0], "p_lt_0p5_cp_high": bci[1], "p_gt_0p75_count": above,
                              "p_gt_0p75_fraction": above / n if n else math.nan, "p_gt_0p75_cp_low": aci[0],
                              "p_gt_0p75_cp_high": aci[1], "observed_mean": float(g.observed.mean()) if n else math.nan,
                              "observed_sample_sd": sd, "null_center_sd": float(g.observed.mean() / sd) if sd > 0 else math.nan})
    vfpr, vdist, paired_rows = pd.DataFrame(fpr_rows), pd.DataFrame(dist_rows), []
    for run in ("R1", "R2"):
        m = rw[(rw.reducer == "mean") & (rw.run == run)].set_index("seed_index")
        d = rw[(rw.reducer == "median") & (rw.run == run)].set_index("seed_index")
        z = m.join(d, lsuffix="_mean", rsuffix="_median", how="inner")
        xm, xd = z.p_two_sided_mean <= 0.05, z.p_two_sided_median <= 0.05
        both, mean_only = int((xm & xd).sum()), int((xm & ~xd).sum())
        median_only, neither = int((~xm & xd).sum()), int((~xm & ~xd).sum())
        delta, lo, hi = newcombe_method10(mean_only, median_only, both, len(z)) if len(z) else (math.nan,) * 3
        paired_rows.append({"run": run, "n": len(z), "both_reject": both, "mean_only_reject": mean_only,
                            "median_only_reject": median_only, "neither_reject": neither,
                            "fpr_mean": (both + mean_only) / len(z) if len(z) else math.nan,
                            "fpr_median": (both + median_only) / len(z) if len(z) else math.nan,
                            "delta_mean_minus_median": delta, "ci_low": lo, "ci_high": hi})
    vpaired = pd.DataFrame(paired_rows)
    for name, got in (("main_fpr.tsv", vfpr), ("main_p_distribution.tsv", vdist), ("main_paired_risk_difference.tsv", vpaired)):
        main = pd.read_csv(output / name, sep="\t")
        if list(main.columns) != list(got.columns) or main.shape != got.shape:
            errors.append(f"{name}: schema or shape mismatch")
            continue
        # Tolerances are derived per (reducer, run) group from the discrepancies that
        # actually occurred: threshold-count columns may move only by the number of worlds
        # whose recomputed and checkpoint p values straddle that threshold; aggregate-fraction
        # and mean_p columns inherit the same group budget. Interval columns (CP, Newcombe)
        # are checked for internal consistency against the counts printed in the same row,
        # so any tolerated count shift uses one counting logic throughout.
        keyed = "reducer" in got.columns
        for _, grow in got.iterrows():
            sel = main[(main.reducer == grow["reducer"]) & (main.run == grow["run"])] if keyed \
                else main[main.run == grow["run"]]
            if len(sel) != 1:
                errors.append(f"{name}: row key mismatch for {grow['run']}"
                              + (f"/{grow['reducer']}" if keyed else ""))
                continue
            mrow = sel.iloc[0]
            if keyed:
                fs = flip_stats[(grow["reducer"], grow["run"])]
                n_g = max(int(grow["n_worlds"]), 1)
                atol_map = {"mean_p": fs["budget"] / n_g + 1e-12,
                            "rejections": float(fs["s005"]),
                            "p_lt_0p5_count": float(fs["s050"]), "p_gt_0p75_count": float(fs["s075"]),
                            "fpr": fs["s005"] / n_g + 1e-12,
                            "p_lt_0p5_fraction": fs["s050"] / n_g + 1e-12,
                            "p_gt_0p75_fraction": fs["s075"] / n_g + 1e-12}
                cp_src = {"cp_low": ("rejections", 0), "cp_high": ("rejections", 1),
                          "p_lt_0p5_cp_low": ("p_lt_0p5_count", 0), "p_lt_0p5_cp_high": ("p_lt_0p5_count", 1),
                          "p_gt_0p75_cp_low": ("p_gt_0p75_count", 0), "p_gt_0p75_cp_high": ("p_gt_0p75_count", 1)}
            else:
                fm, fd = flip_stats[("mean", grow["run"])], flip_stats[("median", grow["run"])]
                ssum = fm["s005"] + fd["s005"]
                n_g = max(int(grow["n"]), 1)
                atol_map = {"both_reject": float(ssum), "mean_only_reject": float(ssum),
                            "median_only_reject": float(ssum), "neither_reject": float(ssum),
                            "fpr_mean": fm["s005"] / n_g + 1e-12, "fpr_median": fd["s005"] / n_g + 1e-12}
                cp_src = {}
            for col in got.columns:
                if col in cp_src:
                    ccol, cidx = cp_src[col]
                    bound = cp_interval(int(mrow[ccol]), int(mrow["n_worlds"]))[cidx]
                    okc = close(float(mrow[col]), float(bound))
                elif col in ("delta_mean_minus_median", "ci_low", "ci_high"):
                    d2, l2, h2 = newcombe_method10(int(mrow["mean_only_reject"]), int(mrow["median_only_reject"]),
                                                  int(mrow["both_reject"]), int(mrow["n"]))
                    okc = close(float(mrow[col]), {"delta_mean_minus_median": d2, "ci_low": l2, "ci_high": h2}[col])
                elif pd.api.types.is_numeric_dtype(got[col]):
                    okc = close(float(mrow[col]), float(grow[col]), atol_map.get(col, 1e-12))
                else:
                    okc = str(mrow[col]) == str(grow[col])
                if not okc:
                    errors.append(f"{name}: mismatch in {col} for {grow['run']}"
                                  + (f"/{grow['reducer']}" if keyed else ""))
    decision = json.loads((output / "main_decision.json").read_text(encoding="utf-8"))
    if mode == "formal":
        complete = expected is not None and len(seeds) == expected and all(
            int(rw[rw.run == run].n_subjects.min()) == (16 if run == "R1" else 13) for run in ("R1", "R2"))
        rule_a = all(r.cp_low <= 0.05 <= r.cp_high for r in vfpr[vfpr.reducer == "mean"].itertuples())
        rule_b = all(r.ci_high < 0 for r in vpaired.itertuples())
        rule_c = all(r.cp_low > 0.05 for r in vfpr[vfpr.reducer == "median"].itertuples())
        md = vdist[vdist.reducer == "mean"]
        rule_d = all(0.40 <= r.mean_p <= 0.60 and r.p_lt_0p5_cp_low <= 0.5 <= r.p_lt_0p5_cp_high and r.p_gt_0p75_cp_low <= 0.25 <= r.p_gt_0p75_cp_high for r in md.itertuples())
        expected_status = "INCOMPLETE" if not complete else ("CONFIRMED" if all((rule_a, rule_b, rule_c, rule_d)) else "NOT_CONFIRMED")
        if decision.get("status") != expected_status:
            errors.append("main decision mismatch")
    elif decision.get("status") != "INCOMPLETE_SMOKE" or decision.get("scientific_rules_evaluated") is not False:
        errors.append("smoke decision improperly evaluates scientific rules")
    report = {"status": "FAIL" if errors else "PASS", "mode": mode, "checkpoints": len(seeds), "errors": errors,
              "checks": ["pair reducers", "pair diagnostics", "exact sign-flip", "FPR/CP", "Newcombe method 10", "p distribution", "main-summary agreement (per-group flip budgets)", "decision"]}
    (output / "validator_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=("smoke", "formal"), required=True)
    args = parser.parse_args()
    report = validate(args.output, args.mode)
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
