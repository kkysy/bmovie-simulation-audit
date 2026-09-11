# Release port of the analysis-workspace tool of the same name; see README in this folder.
"""Primary preregistered summary for the matched J32/J256 robustness stage."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, t

def _repo_root():
    for _q in Path(__file__).resolve().parents:
        if (_q / "MANIFEST.md").is_file():
            return _q
    raise RuntimeError("release repository root not found")
ROOT = _repo_root()
VERIFY_DIR = ROOT / "scripts" / "ieeg" / "paper"
sys.path.insert(0, str(VERIFY_DIR))
from verify_p0_6_paired_contrast import newcombe_method10  # noqa: E402

ARMS = ("J32_reference", "J256_high")
REDUCERS = ("median", "mean")
RUNS = ("R1", "R2")


def cp_interval(k: int, n: int) -> tuple[float, float]:
    return (
        0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1)),
        1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k)),
    )


def paired_cells(x: pd.Series, y: pd.Series) -> tuple[int, int, int, int]:
    xa = x.astype(bool).to_numpy()
    ya = y.astype(bool).to_numpy()
    return (
        int((xa & ya).sum()),
        int((xa & ~ya).sum()),
        int((~xa & ya).sum()),
        int((~xa & ~ya).sum()),
    )


def fpr_table(worlds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        for reducer in REDUCERS:
            for run in RUNS:
                group = worlds[(worlds.arm == arm) & (worlds.reducer == reducer) & (worlds.run == run)]
                reject = group.p_two_sided <= 0.05
                k, n = int(reject.sum()), len(group)
                low, high = cp_interval(k, n) if n else (math.nan, math.nan)
                observed_sd = float(group.observed.std(ddof=1)) if n > 1 else math.nan
                rows.append(
                    {
                        "arm": arm,
                        "J": 32 if arm == "J32_reference" else 256,
                        "reducer": reducer,
                        "run": run,
                        "n_worlds": n,
                        "rejections": k,
                        "fpr": k / n if n else math.nan,
                        "cp_low": low,
                        "cp_high": high,
                        "mean_p": float(group.p_two_sided.mean()) if n else math.nan,
                        "observed_mean": float(group.observed.mean()) if n else math.nan,
                        "observed_sample_sd": observed_sd,
                        "null_center_sd": float(group.observed.mean() / observed_sd) if observed_sd > 0 else math.nan,
                    }
                )
    return pd.DataFrame(rows)


def p_distribution_table(worlds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        for reducer in REDUCERS:
            for run in RUNS:
                group = worlds[(worlds.arm == arm) & (worlds.reducer == reducer) & (worlds.run == run)]
                n = len(group)
                below = int((group.p_two_sided < 0.5).sum())
                above = int((group.p_two_sided > 0.75).sum())
                below_ci = cp_interval(below, n) if n else (math.nan, math.nan)
                above_ci = cp_interval(above, n) if n else (math.nan, math.nan)
                rows.append(
                    {
                        "arm": arm,
                        "reducer": reducer,
                        "run": run,
                        "n_worlds": n,
                        "mean_p": float(group.p_two_sided.mean()) if n else math.nan,
                        "p_lt_0p5_count": below,
                        "p_lt_0p5_fraction": below / n if n else math.nan,
                        "p_lt_0p5_cp_low": below_ci[0],
                        "p_lt_0p5_cp_high": below_ci[1],
                        "p_gt_0p75_count": above,
                        "p_gt_0p75_fraction": above / n if n else math.nan,
                        "p_gt_0p75_cp_low": above_ci[0],
                        "p_gt_0p75_cp_high": above_ci[1],
                    }
                )
    return pd.DataFrame(rows)


def stability_table(worlds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for run in RUNS:
        j32 = worlds[(worlds.arm == "J32_reference") & (worlds.reducer == "mean") & (worlds.run == run)][
            ["seed_index", "observed"]
        ].rename(columns={"observed": "observed_j32"})
        j256 = worlds[(worlds.arm == "J256_high") & (worlds.reducer == "mean") & (worlds.run == run)][
            ["seed_index", "observed"]
        ].rename(columns={"observed": "observed_j256"})
        paired = j32.merge(j256, on="seed_index", validate="one_to_one").sort_values("seed_index")
        difference = paired.observed_j256.to_numpy(float) - paired.observed_j32.to_numpy(float)
        n = len(difference)
        s32 = float(paired.observed_j32.std(ddof=1)) if n > 1 else math.nan
        mean_difference = float(np.mean(difference)) if n else math.nan
        sd_difference = float(np.std(difference, ddof=1)) if n > 1 else math.nan
        sem = sd_difference / math.sqrt(n) if n > 1 else math.nan
        critical = float(t.ppf(0.975, n - 1)) if n > 1 else math.nan
        ci_low = mean_difference - critical * sem if n > 1 else math.nan
        ci_high = mean_difference + critical * sem if n > 1 else math.nan
        rmse = float(np.sqrt(np.mean(difference**2))) if n else math.nan
        correlation = float(np.corrcoef(paired.observed_j32, paired.observed_j256)[0, 1]) if n > 1 else math.nan
        max_position = int(np.argmax(np.abs(difference))) if n else 0
        rows.append(
            {
                "run": run,
                "n": n,
                "s32": s32,
                "mean_difference": mean_difference,
                "sd_difference": sd_difference,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "bias_ratio": mean_difference / s32 if s32 > 0 else math.nan,
                "ci_low_ratio": ci_low / s32 if s32 > 0 else math.nan,
                "ci_high_ratio": ci_high / s32 if s32 > 0 else math.nan,
                "rmse": rmse,
                "rmse_ratio": rmse / s32 if s32 > 0 else math.nan,
                "sd_difference_ratio": sd_difference / s32 if s32 > 0 else math.nan,
                "median_difference": float(np.median(difference)) if n else math.nan,
                "q025_difference": float(np.quantile(difference, 0.025)) if n else math.nan,
                "q25_difference": float(np.quantile(difference, 0.25)) if n else math.nan,
                "q75_difference": float(np.quantile(difference, 0.75)) if n else math.nan,
                "q975_difference": float(np.quantile(difference, 0.975)) if n else math.nan,
                "pearson_r": correlation,
                "max_abs_difference": float(abs(difference[max_position])) if n else math.nan,
                "max_abs_seed_index": int(paired.iloc[max_position].seed_index) if n else -1,
            }
        )
    return pd.DataFrame(rows)


def reducer_paired_table(worlds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for run in RUNS:
        mean_rows = worlds[(worlds.arm == "J256_high") & (worlds.reducer == "mean") & (worlds.run == run)][
            ["seed_index", "p_two_sided"]
        ].rename(columns={"p_two_sided": "p_mean"})
        median_rows = worlds[(worlds.arm == "J256_high") & (worlds.reducer == "median") & (worlds.run == run)][
            ["seed_index", "p_two_sided"]
        ].rename(columns={"p_two_sided": "p_median"})
        paired = mean_rows.merge(median_rows, on="seed_index", validate="one_to_one")
        both, mean_only, median_only, neither = paired_cells(paired.p_mean <= 0.05, paired.p_median <= 0.05)
        delta, low, high = newcombe_method10(mean_only, median_only, both, len(paired)) if len(paired) else (math.nan,) * 3
        rows.append(
            {
                "run": run,
                "n": len(paired),
                "both_reject": both,
                "mean_only_reject": mean_only,
                "median_only_reject": median_only,
                "neither_reject": neither,
                "fpr_mean": (both + mean_only) / len(paired) if len(paired) else math.nan,
                "fpr_median": (both + median_only) / len(paired) if len(paired) else math.nan,
                "delta_mean_minus_median": delta,
                "ci_low": low,
                "ci_high": high,
            }
        )
    return pd.DataFrame(rows)


def j_paired_table(worlds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for reducer in REDUCERS:
        for run in RUNS:
            j32 = worlds[(worlds.arm == "J32_reference") & (worlds.reducer == reducer) & (worlds.run == run)][
                ["seed_index", "p_two_sided"]
            ].rename(columns={"p_two_sided": "p_j32"})
            j256 = worlds[(worlds.arm == "J256_high") & (worlds.reducer == reducer) & (worlds.run == run)][
                ["seed_index", "p_two_sided"]
            ].rename(columns={"p_two_sided": "p_j256"})
            paired = j256.merge(j32, on="seed_index", validate="one_to_one")
            both, j256_only, j32_only, neither = paired_cells(paired.p_j256 <= 0.05, paired.p_j32 <= 0.05)
            delta, low, high = newcombe_method10(j256_only, j32_only, both, len(paired)) if len(paired) else (math.nan,) * 3
            p_difference = paired.p_j256.to_numpy(float) - paired.p_j32.to_numpy(float)
            rows.append(
                {
                    "reducer": reducer,
                    "run": run,
                    "n": len(paired),
                    "both_reject": both,
                    "j256_only_reject": j256_only,
                    "j32_only_reject": j32_only,
                    "neither_reject": neither,
                    "fpr_j256": (both + j256_only) / len(paired) if len(paired) else math.nan,
                    "fpr_j32": (both + j32_only) / len(paired) if len(paired) else math.nan,
                    "delta_j256_minus_j32": delta,
                    "ci_low": low,
                    "ci_high": high,
                    "median_p_difference": float(np.median(p_difference)) if len(paired) else math.nan,
                    "q25_p_difference": float(np.quantile(p_difference, 0.25)) if len(paired) else math.nan,
                    "q75_p_difference": float(np.quantile(p_difference, 0.75)) if len(paired) else math.nan,
                    "max_abs_p_difference": float(np.max(np.abs(p_difference))) if len(paired) else math.nan,
                }
            )
    return pd.DataFrame(rows)


def formal_complete(worlds: pd.DataFrame) -> bool:
    expected_seeds = list(range(1, 121))
    if sorted(worlds.seed_index.unique().tolist()) != expected_seeds:
        return False
    for arm in ARMS:
        for reducer in REDUCERS:
            for run in RUNS:
                group = worlds[(worlds.arm == arm) & (worlds.reducer == reducer) & (worlds.run == run)]
                expected_subjects = 16 if run == "R1" else 13
                if (
                    len(group) != 120
                    or group.seed_index.nunique() != 120
                    or not (group.n_subjects == expected_subjects).all()
                ):
                    return False
    return True


def _json_default(obj):
    # numpy 2.x scalars (bool_, int64) are not JSON serializable; .item() gives the Python scalar
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")


def decision_from_tables(
    complete: bool, fpr: pd.DataFrame, stability: pd.DataFrame, reducer_paired: pd.DataFrame
) -> dict:
    rule_a_by_run = {
        row.run: bool(row.ci_low_ratio >= -0.10 and row.ci_high_ratio <= 0.10 and row.rmse_ratio <= 0.25)
        for row in stability.itertuples()
    }
    rule_b_by_run = {}
    for run in RUNS:
        mean_row = fpr[(fpr.arm == "J256_high") & (fpr.reducer == "mean") & (fpr.run == run)].iloc[0]
        median_row = fpr[(fpr.arm == "J256_high") & (fpr.reducer == "median") & (fpr.run == run)].iloc[0]
        rule_b_by_run[run] = bool(mean_row.cp_low <= 0.05 <= mean_row.cp_high and median_row.cp_low > 0.05)
    rule_c_by_run = {row.run: bool(row.ci_high < 0.0) for row in reducer_paired.itertuples()}
    rules = {
        "A_mean_statistic_stability": rule_a_by_run,
        "B_J256_calibration": rule_b_by_run,
        "C_J256_mean_minus_median": rule_c_by_run,
    }
    all_pass = all(all(values.values()) for values in rules.values())
    hard_a = any(
        row.ci_low_ratio > 0.10 or row.ci_high_ratio < -0.10 or row.rmse_ratio >= 0.50
        for row in stability.itertuples()
    )
    hard_mean = False
    hard_median = False
    for run in RUNS:
        mean_row = fpr[(fpr.arm == "J256_high") & (fpr.reducer == "mean") & (fpr.run == run)].iloc[0]
        median_row = fpr[(fpr.arm == "J256_high") & (fpr.reducer == "median") & (fpr.run == run)].iloc[0]
        hard_mean |= not (mean_row.cp_low <= 0.05 <= mean_row.cp_high)
        hard_median |= median_row.cp_high <= 0.05
    hard_contrast = any(row.ci_low >= 0.0 for row in reducer_paired.itertuples())
    hard_not_robust = bool(hard_a or hard_mean or hard_median or hard_contrast)
    if not complete:
        status = "INCOMPLETE_INVALID"
    elif all_pass:
        status = "ROBUST_TO_SURROGATE_COUNT"
    elif hard_not_robust:
        status = "NOT_ROBUST_TO_SURROGATE_COUNT"
    else:
        status = "INCONCLUSIVE_J_ROBUSTNESS"
    return {
        "analysis_id": "ACC00-MPHASE-JROBUST-v1.0.0",
        "status": status,
        "complete": complete,
        "scientific_rules_evaluated": complete,
        "rules": rules,
        "hard_not_robust": {
            "A_continuous": hard_a,
            "J256_mean_calibration": hard_mean,
            "J256_median_no_elevation": hard_median,
            "J256_contrast_reversal": hard_contrast,
        },
        "prohibition": "Do not change J, world count, seeds, thresholds, intervals, or classification after observing results.",
    }


def summarize(output: Path, mode: str) -> dict:
    worlds = pd.read_csv(output / "validator_world_rows.tsv", sep="\t")
    fpr = fpr_table(worlds)
    distribution = p_distribution_table(worlds)
    stability = stability_table(worlds)
    reducer_paired = reducer_paired_table(worlds)
    j_paired = j_paired_table(worlds)
    fpr.to_csv(output / "main_fpr.tsv", sep="\t", index=False)
    distribution.to_csv(output / "main_p_distribution.tsv", sep="\t", index=False)
    stability.to_csv(output / "main_mean_statistic_stability.tsv", sep="\t", index=False)
    reducer_paired.to_csv(output / "main_reducer_paired_difference.tsv", sep="\t", index=False)
    j_paired.to_csv(output / "main_j_paired_difference.tsv", sep="\t", index=False)
    if mode == "smoke":
        decision = {
            "analysis_id": "ACC00-MPHASE-JROBUST-v1.0.0",
            "status": "INCOMPLETE_SMOKE",
            "complete": False,
            "scientific_rules_evaluated": False,
            "reason": "Two-world smoke validates implementation and timing only.",
        }
    else:
        decision = decision_from_tables(formal_complete(worlds), fpr, stability, reducer_paired)
    (output / "main_decision.json").write_text(json.dumps(decision, indent=2, default=_json_default), encoding="utf-8")
    return decision


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=("smoke", "formal"), required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(args.output, args.mode), indent=2, default=_json_default))


if __name__ == "__main__":
    main()
