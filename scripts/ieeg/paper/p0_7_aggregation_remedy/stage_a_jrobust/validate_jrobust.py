# Release port of the analysis-workspace tool of the same name; see README in this folder.
"""Independent validator for the matched Mphase J32/J256 robustness stage."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, skew, t

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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def cp_interval(k: int, n: int) -> tuple[float, float]:
    return (
        0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1)),
        1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k)),
    )


def exact_signflip(values: np.ndarray) -> tuple[float, float, float]:
    n = len(values)
    rows = np.arange(2**n, dtype=np.uint32)[:, None]
    bits = (rows >> np.arange(n - 1, -1, -1, dtype=np.uint32)) & 1
    signs = 1.0 - 2.0 * bits
    signs[0, :] = 1.0
    statistics = signs @ values / n
    observed = float(statistics[0])
    p_value = float(np.mean(np.abs(statistics) >= abs(observed)))
    sample_sd = float(np.std(statistics, ddof=1))
    return observed, p_value, sample_sd


def close(left: float, right: float, tolerance: float = 1e-12) -> bool:
    return (math.isnan(left) and math.isnan(right)) or abs(left - right) <= tolerance


def close_series(left: pd.Series, right: pd.Series, tolerance: float) -> bool:
    a = left.to_numpy(float)
    b = right.to_numpy(float)
    return bool(np.all(np.isclose(a, b, rtol=0.0, atol=tolerance, equal_nan=True)))


def check_hash(path: Path, expected: str, errors: list[str], label: str) -> None:
    if not path.is_file():
        errors.append(f"missing frozen file: {label}: {path}")
    elif sha256_file(path) != expected.upper():
        errors.append(f"frozen file hash mismatch: {label}: {path}")


def validate_contract_files(contract: dict, errors: list[str]) -> None:
    parent = contract["parent_contract"]
    check_hash(ROOT / parent["path"], parent["sha256"], errors, "parent contract")
    reference = contract["stage_a_reference"]
    for name in ("contract", "runner", "validator", "pair_rows", "subject_rows", "world_rows"):
        item = reference[name]
        check_hash(ROOT / item["path"], item["sha256"], errors, f"Stage A {name}")
    for item in contract["production_code"]:
        check_hash(ROOT / item["path"], item["sha256"], errors, "production code")
    for item in contract["local_code"]:
        check_hash(ROOT / item["path"], item["sha256"], errors, "local code")


def validate_schema(pairs: pd.DataFrame, subjects: pd.DataFrame, worlds: pd.DataFrame, errors: list[str]) -> None:
    pair_columns = {
        "seed_index", "world_seed", "run_kind", "generator_mode", "background_world_sha256",
        "arm", "J_contract", "subject", "run", "session_id", "pair_id", "corrected_mphase",
        "observed_ppc", "mean_refit_surrogate_ppc", "n_phase_events", "quantized_feasible_count",
        "domain_fraction", "surrogate_count", "surrogate_family_mode", "refit_selection_mode",
        "refit_surrogate_count", "mphase_surrogate_reducer", "candidate_delta_sha256",
        "refit_delta_sha256", "pair_stream_initial_state_sha256", "pair_stream_final_state_sha256",
        "pair_signal_sha256", "pair_elapsed_s", "surrogate_ppc_s",
    }
    subject_columns = {
        "seed_index", "world_seed", "run_kind", "arm", "J_contract", "subject", "run", "reducer",
        "Mphase", "n_pairs_total", "n_pairs_finite", "pair_corrected_skewness",
        "pair_positive_fraction", "pair_zero_count",
    }
    world_columns = {
        "seed_index", "world_seed", "run_kind", "arm", "J_contract", "reducer", "run", "n_subjects",
        "observed", "p_two_sided", "studentization_sd", "subject_statistic_skewness", "n_positive",
        "n_negative", "n_zero", "positive_fraction", "effective_contract_sha256", "background_world_sha256",
    }
    for name, frame, required in (("pair", pairs, pair_columns), ("subject", subjects, subject_columns), ("world", worlds, world_columns)):
        missing = sorted(required - set(frame.columns))
        if missing:
            errors.append(f"{name} export missing columns: {missing}")


def validate_identity(
    contract: dict, pairs: pd.DataFrame, subjects: pd.DataFrame, worlds: pd.DataFrame, mode: str, errors: list[str]
) -> list[int]:
    expected_seeds = [1, 2] if mode == "smoke" else list(range(1, 121))
    seeds = sorted(pairs.seed_index.unique().tolist())
    if seeds != expected_seeds:
        errors.append(f"{mode} seed_index set mismatch: {seeds[:5]}...{seeds[-5:] if seeds else []}")
    for frame_name, frame in (("pairs", pairs), ("subjects", subjects), ("worlds", worlds)):
        if sorted(frame.seed_index.unique().tolist()) != seeds:
            errors.append(f"{frame_name} seed set differs from pair export")
        if not (frame.run_kind == mode).all():
            errors.append(f"{frame_name} run_kind mismatch")
        expected_world_seed = contract["seed_namespace"]["base"] + contract["seed_namespace"]["offset"] + frame.seed_index
        if not (frame.world_seed == expected_world_seed).all():
            errors.append(f"{frame_name} world_seed formula mismatch")
    if not (pairs.generator_mode == contract["world"]["generator_mode"]).all():
        errors.append("generator mode is not uniformly benchmark")
    if set(pairs.arm.unique()) != set(ARMS) or set(subjects.arm.unique()) != set(ARMS) or set(worlds.arm.unique()) != set(ARMS):
        errors.append("arm identity mismatch")
    if not (pairs.mphase_surrogate_reducer == "mean").all():
        errors.append("production surrogate reducer is not uniformly mean")
    expected_j = {"J32_reference": 32, "J256_high": 256}
    for arm, value in expected_j.items():
        for frame_name, frame in (("pairs", pairs), ("subjects", subjects), ("worlds", worlds)):
            selected = frame[frame.arm == arm]
            if selected.empty or not (selected.J_contract == value).all():
                errors.append(f"{frame_name} {arm} J identity mismatch")
    for run, expected_subjects in (("R1", 16), ("R2", 13)):
        selected = worlds[worlds.run == run]
        if selected.empty or not (selected.n_subjects == expected_subjects).all():
            errors.append(f"world export {run} subject count mismatch")
    for seed in seeds:
        selected = worlds[worlds.seed_index == seed]
        if selected.background_world_sha256.nunique() != 1:
            errors.append(f"seed {seed}: background world hash differs across arms/reducers")
        hashes = selected.groupby("arm").effective_contract_sha256.nunique()
        if any(hashes.get(arm, 0) != 1 for arm in ARMS):
            errors.append(f"seed {seed}: effective contract hash is not arm-stable")
        arm_hashes = selected.groupby("arm").effective_contract_sha256.first()
        if len(arm_hashes) == 2 and arm_hashes["J32_reference"] == arm_hashes["J256_high"]:
            errors.append(f"seed {seed}: effective J32 and J256 contract hashes are identical")
    return seeds


def validate_pair_crn(pairs: pd.DataFrame, seeds: list[int], errors: list[str]) -> None:
    keys = ["seed_index", "subject", "run", "session_id", "pair_id"]
    comparison_columns = [
        "n_phase_events", "quantized_feasible_count", "domain_fraction", "surrogate_count",
        "surrogate_family_mode", "candidate_delta_sha256", "pair_stream_initial_state_sha256",
        "pair_signal_sha256", "background_world_sha256",
    ]
    for seed in seeds:
        world = pairs[pairs.seed_index == seed]
        j32 = world[world.arm == "J32_reference"].sort_values(keys).reset_index(drop=True)
        j256 = world[world.arm == "J256_high"].sort_values(keys).reset_index(drop=True)
        if len(j32) != 351 or len(j256) != 351:
            errors.append(f"seed {seed}: expected 351 pairs per arm")
            continue
        if not j32[keys].equals(j256[keys]):
            errors.append(f"seed {seed}: pair identities differ across arms")
            continue
        if len(j32) > 1 and not np.array_equal(
            j32.pair_stream_final_state_sha256.astype(str).to_numpy()[:-1],
            j32.pair_stream_initial_state_sha256.astype(str).to_numpy()[1:],
        ):
            errors.append(f"seed {seed}: J32 sequential stream-state chain is broken")
        for column in comparison_columns:
            if pd.api.types.is_numeric_dtype(j32[column]):
                if not close_series(j32[column], j256[column], 1e-12):
                    errors.append(f"seed {seed}: CRN numeric mismatch in {column}")
            elif not (j32[column].astype(str) == j256[column].astype(str)).all():
                errors.append(f"seed {seed}: CRN identity mismatch in {column}")
        for arm, expected_j in (("J32_reference", 32), ("J256_high", 256)):
            arm_rows = world[world.arm == arm]
            expected_count = np.minimum(expected_j, arm_rows.surrogate_count.to_numpy(int))
            if not np.array_equal(arm_rows.refit_surrogate_count.to_numpy(int), expected_count):
                errors.append(f"seed {seed}: {arm} refit count rule failed")
            expected_selection = np.where(arm_rows.surrogate_count.to_numpy(int) <= expected_j, "all_candidates", "uniform_without_replacement")
            if not np.array_equal(arm_rows.refit_selection_mode.astype(str).to_numpy(), expected_selection):
                errors.append(f"seed {seed}: {arm} refit selection mode failed")
        feasible = world.quantized_feasible_count.to_numpy(int)
        candidate_count = world.surrogate_count.to_numpy(int)
        expected_candidate = np.minimum(feasible, 1024)
        if not np.array_equal(candidate_count, expected_candidate):
            errors.append(f"seed {seed}: candidate family count rule failed")
        expected_family_mode = np.where(feasible == 0, "not_estimable", np.where(feasible <= 1024, "enumerated", "sampled"))
        if not np.array_equal(world.surrogate_family_mode.astype(str).to_numpy(), expected_family_mode):
            errors.append(f"seed {seed}: candidate family mode failed")


def compare_stage_a_reference(
    contract: dict, pairs: pd.DataFrame, subjects: pd.DataFrame, worlds: pd.DataFrame, errors: list[str]
) -> None:
    tolerance = float(contract["tolerances"]["stage_a_continuous_abs"])
    new_pairs = pairs[pairs.arm == "J32_reference"].copy()
    new_subjects = subjects[subjects.arm == "J32_reference"].copy()
    new_worlds = worlds[worlds.arm == "J32_reference"].copy()
    seeds = sorted(new_pairs.seed_index.unique())
    old_pairs = pd.read_csv(ROOT / contract["stage_a_reference"]["pair_rows"]["path"], sep="\t")
    old_subjects = pd.read_csv(ROOT / contract["stage_a_reference"]["subject_rows"]["path"], sep="\t")
    old_worlds = pd.read_csv(ROOT / contract["stage_a_reference"]["world_rows"]["path"], sep="\t")
    old_pairs = old_pairs[old_pairs.seed_index.isin(seeds)]
    old_subjects = old_subjects[old_subjects.seed_index.isin(seeds)]
    old_worlds = old_worlds[old_worlds.seed_index.isin(seeds)]
    pair_keys = ["seed_index", "world_seed", "subject", "run", "session_id", "pair_id"]
    pair_text = ["mphase_surrogate_reducer"]
    pair_numeric = [
        "corrected_mphase", "observed_ppc", "mean_refit_surrogate_ppc", "n_phase_events",
        "domain_fraction", "surrogate_count", "refit_surrogate_count",
    ]
    compare_frames(old_pairs, new_pairs, pair_keys, pair_text, pair_numeric, tolerance, "Stage A pair", errors)
    subject_keys = ["seed_index", "world_seed", "subject", "run", "reducer"]
    subject_numeric = [
        "Mphase", "n_pairs_total", "n_pairs_finite", "pair_corrected_skewness",
        "pair_positive_fraction", "pair_zero_count",
    ]
    compare_frames(old_subjects, new_subjects, subject_keys, [], subject_numeric, tolerance, "Stage A subject", errors)
    world_keys = ["seed_index", "world_seed", "reducer", "run"]
    world_numeric = [
        "n_subjects", "observed", "p_two_sided", "studentization_sd", "subject_statistic_skewness",
        "n_positive", "n_negative", "n_zero", "positive_fraction",
    ]
    compare_frames(old_worlds, new_worlds, world_keys, [], world_numeric, tolerance, "Stage A world", errors)


def compare_frames(
    old: pd.DataFrame,
    new: pd.DataFrame,
    keys: list[str],
    text_columns: list[str],
    numeric_columns: list[str],
    tolerance: float,
    label: str,
    errors: list[str],
) -> None:
    old_sorted = old.sort_values(keys).reset_index(drop=True)
    new_sorted = new.sort_values(keys).reset_index(drop=True)
    if len(old_sorted) != len(new_sorted) or not old_sorted[keys].equals(new_sorted[keys]):
        errors.append(f"{label} identity/cardinality mismatch")
        return
    for column in text_columns:
        if not (old_sorted[column].astype(str) == new_sorted[column].astype(str)).all():
            errors.append(f"{label} text mismatch in {column}")
    for column in numeric_columns:
        if not close_series(old_sorted[column], new_sorted[column], tolerance):
            errors.append(f"{label} numeric mismatch in {column}")


def recompute_subjects_and_worlds(
    contract: dict,
    pairs: pd.DataFrame,
    subjects: pd.DataFrame,
    worlds: pd.DataFrame,
    seeds: list[int],
    errors: list[str],
) -> pd.DataFrame:
    tolerance = float(contract["tolerances"]["continuous_abs"])
    recomputed = []
    for seed in seeds:
        for arm in ARMS:
            world_pairs = pairs[(pairs.seed_index == seed) & (pairs.arm == arm)]
            for run in RUNS:
                run_pairs = world_pairs[world_pairs.run == run]
                for reducer in REDUCERS:
                    values = []
                    for subject, group in run_pairs.groupby("subject", sort=True):
                        pair_values = group.corrected_mphase.to_numpy(float)
                        pair_values = pair_values[np.isfinite(pair_values)]
                        statistic = float(np.median(pair_values) if reducer == "median" else np.mean(pair_values))
                        pair_skew = float(skew(pair_values, bias=False)) if len(pair_values) >= 3 and not np.all(pair_values == pair_values[0]) else math.nan
                        selected = subjects[
                            (subjects.seed_index == seed) & (subjects.arm == arm) & (subjects.run == run)
                            & (subjects.reducer == reducer) & (subjects.subject == subject)
                        ]
                        if len(selected) != 1:
                            errors.append(f"{seed}/{arm}/{run}/{reducer}/{subject}: subject row cardinality")
                        else:
                            row = selected.iloc[0]
                            if not close(statistic, float(row.Mphase), tolerance):
                                errors.append(f"{seed}/{arm}/{run}/{reducer}/{subject}: subject reducer mismatch")
                            if int(row.n_pairs_finite) != len(pair_values) or not close(pair_skew, float(row.pair_corrected_skewness), tolerance):
                                errors.append(f"{seed}/{arm}/{run}/{reducer}/{subject}: subject diagnostics mismatch")
                        values.append(statistic)
                    observed, p_value, sample_sd = exact_signflip(np.asarray(values, float))
                    selected_world = worlds[
                        (worlds.seed_index == seed) & (worlds.arm == arm) & (worlds.run == run) & (worlds.reducer == reducer)
                    ]
                    if len(selected_world) != 1:
                        errors.append(f"{seed}/{arm}/{run}/{reducer}: world row cardinality")
                    else:
                        row = selected_world.iloc[0]
                        p_tolerance = 1 / (2 ** len(values)) + float(contract["tolerances"]["signflip_one_flip_plus_abs"])
                        if not close(observed, float(row.observed), tolerance) or abs(p_value - float(row.p_two_sided)) > p_tolerance or not close(sample_sd, float(row.studentization_sd), tolerance):
                            errors.append(f"{seed}/{arm}/{run}/{reducer}: exact sign-flip mismatch")
                    recomputed.append(
                        {"seed_index": seed, "arm": arm, "reducer": reducer, "run": run,
                         "n_subjects": len(values), "observed": observed, "p_two_sided": p_value,
                         "studentization_sd": sample_sd}
                    )
    return pd.DataFrame(recomputed)


def independent_tables(worlds: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fpr_rows = []
    distribution_rows = []
    for arm in ARMS:
        for reducer in REDUCERS:
            for run in RUNS:
                group = worlds[(worlds.arm == arm) & (worlds.reducer == reducer) & (worlds.run == run)]
                n = len(group)
                k = int((group.p_two_sided <= 0.05).sum())
                low, high = cp_interval(k, n) if n else (math.nan, math.nan)
                observed_sd = float(group.observed.std(ddof=1)) if n > 1 else math.nan
                fpr_rows.append(
                    {"arm": arm, "J": 32 if arm == "J32_reference" else 256, "reducer": reducer, "run": run,
                     "n_worlds": n, "rejections": k, "fpr": k / n if n else math.nan, "cp_low": low,
                     "cp_high": high, "mean_p": float(group.p_two_sided.mean()) if n else math.nan,
                     "observed_mean": float(group.observed.mean()) if n else math.nan,
                     "observed_sample_sd": observed_sd,
                     "null_center_sd": float(group.observed.mean() / observed_sd) if observed_sd > 0 else math.nan}
                )
                below = int((group.p_two_sided < 0.5).sum())
                above = int((group.p_two_sided > 0.75).sum())
                below_ci = cp_interval(below, n) if n else (math.nan, math.nan)
                above_ci = cp_interval(above, n) if n else (math.nan, math.nan)
                distribution_rows.append(
                    {"arm": arm, "reducer": reducer, "run": run, "n_worlds": n,
                     "mean_p": float(group.p_two_sided.mean()) if n else math.nan,
                     "p_lt_0p5_count": below, "p_lt_0p5_fraction": below / n if n else math.nan,
                     "p_lt_0p5_cp_low": below_ci[0], "p_lt_0p5_cp_high": below_ci[1],
                     "p_gt_0p75_count": above, "p_gt_0p75_fraction": above / n if n else math.nan,
                     "p_gt_0p75_cp_low": above_ci[0], "p_gt_0p75_cp_high": above_ci[1]}
                )
    stability_rows = []
    reducer_paired_rows = []
    j_paired_rows = []
    for run in RUNS:
        j32 = worlds[(worlds.arm == "J32_reference") & (worlds.reducer == "mean") & (worlds.run == run)][["seed_index", "observed"]].rename(columns={"observed": "j32"})
        j256 = worlds[(worlds.arm == "J256_high") & (worlds.reducer == "mean") & (worlds.run == run)][["seed_index", "observed"]].rename(columns={"observed": "j256"})
        paired = j32.merge(j256, on="seed_index", validate="one_to_one").sort_values("seed_index")
        difference = paired.j256.to_numpy(float) - paired.j32.to_numpy(float)
        n = len(difference)
        s32 = float(paired.j32.std(ddof=1)) if n > 1 else math.nan
        mean_difference = float(np.mean(difference)) if n else math.nan
        sd_difference = float(np.std(difference, ddof=1)) if n > 1 else math.nan
        critical = float(t.ppf(0.975, n - 1)) if n > 1 else math.nan
        half_width = critical * sd_difference / math.sqrt(n) if n > 1 else math.nan
        ci_low, ci_high = mean_difference - half_width, mean_difference + half_width
        rmse = float(np.sqrt(np.mean(difference**2))) if n else math.nan
        max_position = int(np.argmax(np.abs(difference))) if n else 0
        stability_rows.append(
            {"run": run, "n": n, "s32": s32, "mean_difference": mean_difference, "sd_difference": sd_difference,
             "ci_low": ci_low, "ci_high": ci_high, "bias_ratio": mean_difference / s32 if s32 > 0 else math.nan,
             "ci_low_ratio": ci_low / s32 if s32 > 0 else math.nan, "ci_high_ratio": ci_high / s32 if s32 > 0 else math.nan,
             "rmse": rmse, "rmse_ratio": rmse / s32 if s32 > 0 else math.nan,
             "sd_difference_ratio": sd_difference / s32 if s32 > 0 else math.nan,
             "median_difference": float(np.median(difference)) if n else math.nan,
             "q025_difference": float(np.quantile(difference, 0.025)) if n else math.nan,
             "q25_difference": float(np.quantile(difference, 0.25)) if n else math.nan,
             "q75_difference": float(np.quantile(difference, 0.75)) if n else math.nan,
             "q975_difference": float(np.quantile(difference, 0.975)) if n else math.nan,
             "pearson_r": float(np.corrcoef(paired.j32, paired.j256)[0, 1]) if n > 1 else math.nan,
             "max_abs_difference": float(abs(difference[max_position])) if n else math.nan,
             "max_abs_seed_index": int(paired.iloc[max_position].seed_index) if n else -1}
        )
        mean_rows = worlds[(worlds.arm == "J256_high") & (worlds.reducer == "mean") & (worlds.run == run)][["seed_index", "p_two_sided"]].rename(columns={"p_two_sided": "mean"})
        median_rows = worlds[(worlds.arm == "J256_high") & (worlds.reducer == "median") & (worlds.run == run)][["seed_index", "p_two_sided"]].rename(columns={"p_two_sided": "median"})
        paired_reducers = mean_rows.merge(median_rows, on="seed_index", validate="one_to_one")
        both, mean_only, median_only, neither = paired_binary_cells(paired_reducers["mean"] <= 0.05, paired_reducers["median"] <= 0.05)
        delta, low, high = newcombe_method10(mean_only, median_only, both, len(paired_reducers)) if len(paired_reducers) else (math.nan,) * 3
        reducer_paired_rows.append(
            {"run": run, "n": len(paired_reducers), "both_reject": both, "mean_only_reject": mean_only,
             "median_only_reject": median_only, "neither_reject": neither,
             "fpr_mean": (both + mean_only) / len(paired_reducers) if len(paired_reducers) else math.nan,
             "fpr_median": (both + median_only) / len(paired_reducers) if len(paired_reducers) else math.nan,
             "delta_mean_minus_median": delta, "ci_low": low, "ci_high": high}
        )
    for reducer in REDUCERS:
        for run in RUNS:
            j32 = worlds[(worlds.arm == "J32_reference") & (worlds.reducer == reducer) & (worlds.run == run)][["seed_index", "p_two_sided"]].rename(columns={"p_two_sided": "j32"})
            j256 = worlds[(worlds.arm == "J256_high") & (worlds.reducer == reducer) & (worlds.run == run)][["seed_index", "p_two_sided"]].rename(columns={"p_two_sided": "j256"})
            paired = j256.merge(j32, on="seed_index", validate="one_to_one")
            both, j256_only, j32_only, neither = paired_binary_cells(paired.j256 <= 0.05, paired.j32 <= 0.05)
            delta, low, high = newcombe_method10(j256_only, j32_only, both, len(paired)) if len(paired) else (math.nan,) * 3
            p_difference = paired.j256.to_numpy(float) - paired.j32.to_numpy(float)
            j_paired_rows.append(
                {"reducer": reducer, "run": run, "n": len(paired), "both_reject": both,
                 "j256_only_reject": j256_only, "j32_only_reject": j32_only, "neither_reject": neither,
                 "fpr_j256": (both + j256_only) / len(paired) if len(paired) else math.nan,
                 "fpr_j32": (both + j32_only) / len(paired) if len(paired) else math.nan,
                 "delta_j256_minus_j32": delta, "ci_low": low, "ci_high": high,
                 "median_p_difference": float(np.median(p_difference)) if len(paired) else math.nan,
                 "q25_p_difference": float(np.quantile(p_difference, 0.25)) if len(paired) else math.nan,
                 "q75_p_difference": float(np.quantile(p_difference, 0.75)) if len(paired) else math.nan,
                 "max_abs_p_difference": float(np.max(np.abs(p_difference))) if len(paired) else math.nan}
            )
    return (
        pd.DataFrame(fpr_rows), pd.DataFrame(distribution_rows), pd.DataFrame(stability_rows),
        pd.DataFrame(reducer_paired_rows), pd.DataFrame(j_paired_rows),
    )


def paired_binary_cells(left: pd.Series, right: pd.Series) -> tuple[int, int, int, int]:
    a = left.astype(bool).to_numpy()
    b = right.astype(bool).to_numpy()
    return int((a & b).sum()), int((a & ~b).sum()), int((~a & b).sum()), int((~a & ~b).sum())


def compare_summary_file(path: Path, expected: pd.DataFrame, keys: list[str], errors: list[str]) -> None:
    if not path.is_file():
        errors.append(f"missing main summary: {path.name}")
        return
    actual = pd.read_csv(path, sep="\t")
    actual = actual.sort_values(keys).reset_index(drop=True)
    expected = expected.sort_values(keys).reset_index(drop=True)
    if list(actual.columns) != list(expected.columns) or len(actual) != len(expected):
        errors.append(f"summary schema/cardinality mismatch: {path.name}")
        return
    for column in actual.columns:
        if pd.api.types.is_numeric_dtype(actual[column]):
            if not close_series(actual[column], expected[column], 1e-12):
                errors.append(f"summary numeric mismatch: {path.name}:{column}")
        elif not (actual[column].astype(str) == expected[column].astype(str)).all():
            errors.append(f"summary text mismatch: {path.name}:{column}")


def decision_from_tables(complete: bool, fpr: pd.DataFrame, stability: pd.DataFrame, reducer_paired: pd.DataFrame) -> dict:
    rule_a = {row.run: bool(row.ci_low_ratio >= -0.10 and row.ci_high_ratio <= 0.10 and row.rmse_ratio <= 0.25) for row in stability.itertuples()}
    rule_b = {}
    for run in RUNS:
        mean_row = fpr[(fpr.arm == "J256_high") & (fpr.reducer == "mean") & (fpr.run == run)].iloc[0]
        median_row = fpr[(fpr.arm == "J256_high") & (fpr.reducer == "median") & (fpr.run == run)].iloc[0]
        rule_b[run] = bool(mean_row.cp_low <= 0.05 <= mean_row.cp_high and median_row.cp_low > 0.05)
    rule_c = {row.run: bool(row.ci_high < 0.0) for row in reducer_paired.itertuples()}
    rules = {"A_mean_statistic_stability": rule_a, "B_J256_calibration": rule_b, "C_J256_mean_minus_median": rule_c}
    all_pass = all(all(values.values()) for values in rules.values())
    hard_a = any(row.ci_low_ratio > 0.10 or row.ci_high_ratio < -0.10 or row.rmse_ratio >= 0.50 for row in stability.itertuples())
    hard_mean = False
    hard_median = False
    for run in RUNS:
        mean_row = fpr[(fpr.arm == "J256_high") & (fpr.reducer == "mean") & (fpr.run == run)].iloc[0]
        median_row = fpr[(fpr.arm == "J256_high") & (fpr.reducer == "median") & (fpr.run == run)].iloc[0]
        hard_mean |= not (mean_row.cp_low <= 0.05 <= mean_row.cp_high)
        hard_median |= median_row.cp_high <= 0.05
    hard_contrast = any(row.ci_low >= 0.0 for row in reducer_paired.itertuples())
    if not complete:
        status = "INCOMPLETE_INVALID"
    elif all_pass:
        status = "ROBUST_TO_SURROGATE_COUNT"
    elif hard_a or hard_mean or hard_median or hard_contrast:
        status = "NOT_ROBUST_TO_SURROGATE_COUNT"
    else:
        status = "INCONCLUSIVE_J_ROBUSTNESS"
    return {
        "status": status,
        "complete": complete,
        "rules": rules,
        "hard_not_robust": {"A_continuous": hard_a, "J256_mean_calibration": hard_mean,
                            "J256_median_no_elevation": hard_median, "J256_contrast_reversal": hard_contrast},
    }


def validate(output: Path, mode: str) -> dict:
    contract_path = Path(__file__).with_name("confirm_mphase_jrobust_contract.json")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    errors: list[str] = []
    validate_contract_files(contract, errors)
    pairs = pd.read_csv(output / "validator_pair_rows.tsv", sep="\t")
    subjects = pd.read_csv(output / "validator_subject_rows.tsv", sep="\t")
    worlds = pd.read_csv(output / "validator_world_rows.tsv", sep="\t")
    validate_schema(pairs, subjects, worlds, errors)
    seeds = validate_identity(contract, pairs, subjects, worlds, mode, errors)
    validate_pair_crn(pairs, seeds, errors)
    compare_stage_a_reference(contract, pairs, subjects, worlds, errors)
    recomputed_worlds = recompute_subjects_and_worlds(contract, pairs, subjects, worlds, seeds, errors)
    fpr, distribution, stability, reducer_paired, j_paired = independent_tables(worlds)
    compare_summary_file(output / "main_fpr.tsv", fpr, ["arm", "reducer", "run"], errors)
    compare_summary_file(output / "main_p_distribution.tsv", distribution, ["arm", "reducer", "run"], errors)
    compare_summary_file(output / "main_mean_statistic_stability.tsv", stability, ["run"], errors)
    compare_summary_file(output / "main_reducer_paired_difference.tsv", reducer_paired, ["run"], errors)
    compare_summary_file(output / "main_j_paired_difference.tsv", j_paired, ["reducer", "run"], errors)
    complete = mode == "formal" and seeds == list(range(1, 121)) and len(recomputed_worlds) == 960
    decision = json.loads((output / "main_decision.json").read_text(encoding="utf-8"))
    if mode == "smoke":
        if decision.get("status") != "INCOMPLETE_SMOKE" or decision.get("scientific_rules_evaluated") is not False:
            errors.append("smoke decision improperly evaluates scientific rules")
    else:
        expected = decision_from_tables(complete, fpr, stability, reducer_paired)
        if decision.get("status") != expected["status"] or decision.get("complete") != expected["complete"]:
            errors.append("formal decision status/completeness mismatch")
        if decision.get("rules") != expected["rules"] or decision.get("hard_not_robust") != expected["hard_not_robust"]:
            errors.append("formal decision rule mismatch")
    validated_status = "INCOMPLETE_INVALID" if errors else decision.get("status")
    validated_decision = {
        "analysis_id": contract["analysis_id"],
        "status": validated_status,
        "main_decision_status": decision.get("status"),
        "validator_status": "FAIL" if errors else "PASS",
        "scientific_interpretation_allowed": not errors and mode == "formal",
        "errors": errors,
    }
    (output / "validated_decision.json").write_text(json.dumps(validated_decision, indent=2), encoding="utf-8")
    report = {
        "status": "FAIL" if errors else "PASS",
        "mode": mode,
        "matched_checkpoints": len(seeds),
        "validated_decision_status": validated_status,
        "errors": errors,
        "checks": [
            "frozen parent/production/local/Stage A hashes",
            "seed/arm/J identities",
            "pair-level CRN and candidate-family hashes",
            "B/J enumeration and sampling counts",
            "Stage A J32 pair/subject/world reproduction",
            "pair-to-subject median and mean recomputation",
            "exact sign-flip column 2",
            "FPR and Clopper-Pearson intervals",
            "Newcombe method 10 paired contrasts",
            "Rule A Student-t equivalence and RMSE",
            "main-summary agreement",
            "three-branch decision",
        ],
    }
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
