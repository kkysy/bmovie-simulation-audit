# Release port of the analysis-workspace tool of the same name; see README in this folder.
"""Primary Stage A summary from checkpoint-recorded world statistics."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import beta

def _repo_root():
    for _q in Path(__file__).resolve().parents:
        if (_q / "MANIFEST.md").is_file():
            return _q
    raise RuntimeError("release repository root not found")
ROOT = _repo_root()
VERIFY_DIR = ROOT / "scripts" / "ieeg" / "paper"
sys.path.insert(0, str(VERIFY_DIR))
from verify_p0_6_paired_contrast import newcombe_method10  # noqa: E402


def cp_interval(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def paired_cells(x: pd.Series, y: pd.Series) -> tuple[int, int, int, int]:
    x = x.astype(bool).to_numpy()
    y = y.astype(bool).to_numpy()
    return int((x & y).sum()), int((x & ~y).sum()), int((~x & y).sum()), int((~x & ~y).sum())


def summarize(output: Path, mode: str) -> dict:
    worlds = pd.read_csv(output / "validator_world_rows.tsv", sep="\t")
    expected = 120 if mode == "formal" else None
    fpr_rows, dist_rows, complete = [], [], True
    for reducer in ("median", "mean"):
        for run in ("R1", "R2"):
            g = worlds[(worlds.reducer == reducer) & (worlds.run == run)].sort_values("seed_index")
            if expected is not None:
                complete &= len(g) == expected and g.seed_index.nunique() == expected
                complete &= int(g.n_subjects.min()) == (16 if run == "R1" else 13)
            reject = g.p_two_sided <= 0.05
            k, n = int(reject.sum()), len(g)
            lo, hi = cp_interval(k, n) if n else (float("nan"), float("nan"))
            fpr_rows.append({"reducer": reducer, "run": run, "n_worlds": n, "rejections": k,
                             "fpr": k / n if n else float("nan"), "cp_low": lo, "cp_high": hi})
            below, above = int((g.p_two_sided < 0.5).sum()), int((g.p_two_sided > 0.75).sum())
            bci = cp_interval(below, n) if n else (float("nan"), float("nan"))
            aci = cp_interval(above, n) if n else (float("nan"), float("nan"))
            sd = float(g.observed.std(ddof=1)) if n > 1 else float("nan")
            dist_rows.append({"reducer": reducer, "run": run, "n_worlds": n,
                              "mean_p": float(g.p_two_sided.mean()) if n else float("nan"),
                              "p_lt_0p5_count": below, "p_lt_0p5_fraction": below / n if n else float("nan"),
                              "p_lt_0p5_cp_low": bci[0], "p_lt_0p5_cp_high": bci[1],
                              "p_gt_0p75_count": above, "p_gt_0p75_fraction": above / n if n else float("nan"),
                              "p_gt_0p75_cp_low": aci[0], "p_gt_0p75_cp_high": aci[1],
                              "observed_mean": float(g.observed.mean()) if n else float("nan"),
                              "observed_sample_sd": sd,
                              "null_center_sd": float(g.observed.mean() / sd) if sd > 0 else float("nan")})

    fpr, dist, paired_rows = pd.DataFrame(fpr_rows), pd.DataFrame(dist_rows), []
    for run in ("R1", "R2"):
        m = worlds[(worlds.reducer == "mean") & (worlds.run == run)][["seed_index", "p_two_sided"]]
        d = worlds[(worlds.reducer == "median") & (worlds.run == run)][["seed_index", "p_two_sided"]]
        z = m.merge(d, on="seed_index", suffixes=("_mean", "_median"), validate="one_to_one")
        both, mean_only, median_only, neither = paired_cells(z.p_two_sided_mean <= 0.05, z.p_two_sided_median <= 0.05)
        delta, lo, hi = newcombe_method10(mean_only, median_only, both, len(z)) if len(z) else (float("nan"),) * 3
        paired_rows.append({"run": run, "n": len(z), "both_reject": both, "mean_only_reject": mean_only,
                            "median_only_reject": median_only, "neither_reject": neither,
                            "fpr_mean": (both + mean_only) / len(z) if len(z) else float("nan"),
                            "fpr_median": (both + median_only) / len(z) if len(z) else float("nan"),
                            "delta_mean_minus_median": delta, "ci_low": lo, "ci_high": hi})
    paired = pd.DataFrame(paired_rows)
    fpr.to_csv(output / "main_fpr.tsv", sep="\t", index=False)
    dist.to_csv(output / "main_p_distribution.tsv", sep="\t", index=False)
    paired.to_csv(output / "main_paired_risk_difference.tsv", sep="\t", index=False)
    if mode != "formal":
        decision = {"analysis_id": "ACC00-MPHASE-SUBJECT-REDUCER-CONFIRM-STAGE-A-v1.0.0",
                    "status": "INCOMPLETE_SMOKE", "scientific_rules_evaluated": False,
                    "reason": "Smoke validates implementation only; it is not a formal Stage A result."}
    else:
        rule_a = all(r.cp_low <= 0.05 <= r.cp_high for r in fpr[fpr.reducer == "mean"].itertuples())
        rule_b = all(r.ci_high < 0 for r in paired.itertuples())
        rule_c = all(r.cp_low > 0.05 for r in fpr[fpr.reducer == "median"].itertuples())
        md = dist[dist.reducer == "mean"]
        rule_d = all(0.40 <= r.mean_p <= 0.60 and r.p_lt_0p5_cp_low <= 0.5 <= r.p_lt_0p5_cp_high
                     and r.p_gt_0p75_cp_low <= 0.25 <= r.p_gt_0p75_cp_high for r in md.itertuples())
        rules = {"A_mean_calibration": rule_a, "B_paired_improvement": rule_b,
                 "C_median_miscalibration_sanity": rule_c, "D_mean_p_uniformity": rule_d}
        decision = {"analysis_id": "ACC00-MPHASE-SUBJECT-REDUCER-CONFIRM-STAGE-A-v1.0.0",
                    "status": "INCOMPLETE" if not complete else ("CONFIRMED" if all(rules.values()) else "NOT_CONFIRMED"),
                    "complete": complete, "rules": rules,
                    "prohibition": "Do not tune or add worlds after this decision; stop before any power-arm design."}
    (output / "main_decision.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")
    return decision


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=("smoke", "formal"), required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(args.output, args.mode), indent=2))


if __name__ == "__main__":
    main()
