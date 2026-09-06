#!/usr/bin/env python3
"""Read-only aggregation for the frozen P0-2 Mpower simulation power curve.

The main mode refuses incomplete or identity-inconsistent runs before creating
any successful paper table.  MATLAB is used only to decode its documented
``load`` representation of compact string/table objects; h5py reads all large
numeric sign-flip arrays directly.  No #refs# names are interpreted.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np
import scipy
from scipy.special import betaincinv


ROOT = Path(__file__).resolve().parents[3]
PAPER = ROOT / "scripts/ieeg/paper/acc00_sim_powercurve_contract.json"
PARENT = ROOT / "scripts/ieeg/acc00_sim_contract.json"
MAIN = ROOT / "processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/profiling/main"
V2 = ROOT / "processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/profiling/v2"
NULL = ROOT / "processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/nullpool_edgeguard_20260829"
OUT = ROOT / "processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables"
RUN_MANIFEST = ROOT / "processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/run_manifest.json"
ALPHA = 0.05
HISTORICAL_V2_PROFILING_CONTRACT_SHA = "60329E5E93664CEF1B5FCFBD9F353DDC28CCA59E857E0A5D51C75815F70BFCD2"
METADATA_CHUNK_SIZE = 100
METADATA_CHUNK_TIMEOUT_S = 300


class AggregationError(RuntimeError):
    """An unusable checkpoint collection; not a scientific decision."""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def now_utc() -> str:
    return datetime.now(UTC).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def scalar(h: h5py.File, key: str) -> float:
    if key not in h:
        raise AggregationError(f"missing HDF5 numeric field: {key}")
    x = np.asarray(h[key][()])
    if x.size != 1:
        raise AggregationError(f"non-scalar HDF5 field: {key}, shape={x.shape}")
    value = float(x.reshape(-1)[0])
    if not math.isfinite(value):
        raise AggregationError(f"non-finite HDF5 field: {key}")
    return value


def matlab_char(h: h5py.File, key: str) -> str:
    """Decode plain MATLAB char only; MATLAB string/table are deliberately not guessed."""
    if key not in h:
        raise AggregationError(f"missing HDF5 char field: {key}")
    d = h[key]
    if d.attrs.get("MATLAB_class", b"") != b"char":
        raise AggregationError(f"expected MATLAB char, not compact object: {key}")
    values = np.asarray(d[()]).reshape(-1, order="F")
    try:
        return "".join(chr(int(v)) for v in values)
    except (TypeError, ValueError) as exc:
        raise AggregationError(f"cannot decode MATLAB char: {key}") from exc


def cp_interval(k: int, n: int) -> tuple[float, float]:
    if not 0 <= k <= n or n <= 0:
        raise ValueError(f"invalid binomial count k={k}, n={n}")
    lo = 0.0 if k == 0 else float(betaincinv(k, n - k + 1, 0.025))
    hi = 1.0 if k == n else float(betaincinv(k + 1, n - k, 0.975))
    return lo, hi


def mcse(k: int, n: int) -> float:
    p = k / n
    return math.sqrt(p * (1.0 - p) / n)


def theta_grid(contract: dict[str, Any]) -> np.ndarray:
    inv = contract["theta_inversion"]
    lower, upper, step = (float(inv[k]) for k in ("lower", "upper", "step"))
    count_float = (upper - lower) / step
    count = round(count_float)
    if not math.isclose(count_float, count, rel_tol=0.0, abs_tol=1e-10):
        raise AggregationError("theta inversion upper/lower/step do not form an exact grid")
    grid = lower + np.arange(count + 1, dtype=np.float64) * step
    if not math.isclose(float(grid[-1]), upper, rel_tol=0.0, abs_tol=2e-15):
        raise AggregationError("theta inversion grid endpoint mismatch")
    return grid


def p_curve_from_arrays(signs: np.ndarray, statistics: np.ndarray, observed: float,
                        grid: np.ndarray, block: int = 64) -> np.ndarray:
    """Exact Mpower sign-flip inversion using the frozen translation identity."""
    if signs.ndim != 2 or statistics.ndim != 1 or signs.shape[1] != statistics.size:
        raise AggregationError("sign/statistic matrix dimensions are inconsistent")
    if not np.isfinite(signs).all() or not np.isfinite(statistics).all() or not math.isfinite(observed):
        raise AggregationError("non-finite exact sign-flip inputs")
    means = signs.mean(axis=0)
    p = np.empty(grid.size, dtype=np.float64)
    for begin in range(0, grid.size, block):
        theta = grid[begin:begin + block]
        translated = statistics[:, None] - means[:, None] * theta[None, :]
        observed_shifted = observed - theta
        # The frozen rule has no floating-point tolerance and includes all-plus.
        p[begin:begin + theta.size] = np.count_nonzero(
            np.abs(translated) >= np.abs(observed_shifted)[None, :], axis=0
        ) / statistics.size
    return p


def coverage_from_arrays(signs: np.ndarray, statistics: np.ndarray, observed: float,
                         truth: float, contract: dict[str, Any]) -> dict[str, Any]:
    grid = theta_grid(contract)
    p = p_curve_from_arrays(signs, statistics, observed, grid)
    accept = p > ALPHA
    hits = np.flatnonzero(accept)
    nearest = int(np.argmin(np.abs(grid - truth)))
    distance = abs(float(grid[nearest]) - truth)
    hit_tolerance = float(contract["theta_inversion"]["hit_tolerance"])
    return {
        "ci_theta_low": float(grid[hits[0]]) if hits.size else math.nan,
        "ci_theta_high": float(grid[hits[-1]]) if hits.size else math.nan,
        "ci_noncontiguous": bool(hits.size > 1 and np.any(np.diff(hits) != 1)),
        "theta_at_truth_grid": float(grid[nearest]),
        "p_at_truth_grid": float(p[nearest]),
        "truth_grid_distance": distance,
        "truth_in_ci": bool(distance <= hit_tolerance and accept[nearest]),
    }


def validate_signflip(h: h5py.File, run: str, n_subjects: int) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Confirm stored orientation from all-plus before extracting Mpower (metric index 0)."""
    base = f"checkpoint/inference/{run}"
    needed = ("signs", "statistics", "observed", "p_two_sided")
    for field in needed:
        if f"{base}/{field}" not in h:
            raise AggregationError(f"{run}: missing exact sign-flip field {field}")
    signs = np.asarray(h[f"{base}/signs"][()], dtype=np.float64)
    stats = np.asarray(h[f"{base}/statistics"][()], dtype=np.float64)
    observed = np.asarray(h[f"{base}/observed"][()], dtype=np.float64).reshape(-1, order="F")
    p_two = np.asarray(h[f"{base}/p_two_sided"][()], dtype=np.float64).reshape(-1, order="F")
    expected_enum = 2 ** n_subjects
    if signs.shape != (n_subjects, expected_enum):
        raise AggregationError(f"{run}: signs orientation/shape {signs.shape}, expected {(n_subjects, expected_enum)}")
    if stats.shape != (2, expected_enum) or observed.shape != (2,) or p_two.shape != (2,):
        raise AggregationError(f"{run}: statistic orientation/shape is not metrics x enumerations")
    if not all(np.isfinite(x).all() for x in (signs, stats, observed, p_two)):
        raise AggregationError(f"{run}: non-finite exact sign-flip field")
    if not np.all(np.isin(signs, (-1.0, 1.0))):
        raise AggregationError(f"{run}: sign matrix contains values other than -1/+1")
    # This checks orientation from data instead of assuming MATLAB's HDF5 transpose.
    if not np.all(signs[:, 0] == 1.0):
        raise AggregationError(f"{run}: first enumerated sign vector is not all-plus")
    if not np.array_equal(stats[:, 0], observed):
        raise AggregationError(f"{run}: all-plus statistic does not equal observed vector")
    for metric in range(2):
        expected_p = float(np.count_nonzero(np.abs(stats[metric]) >= abs(observed[metric])) / expected_enum)
        if p_two[metric] != expected_p:
            raise AggregationError(f"{run}: stored p_two_sided metric {metric} disagrees with raw exact enumeration")
    return float(observed[0]), float(p_two[0]), signs, stats[0]


def list_main_files(contract: dict[str, Any]) -> tuple[list[tuple[int, int, Path]], list[str]]:
    files: list[tuple[int, int, Path]] = []
    errors: list[str] = []
    counts = contract["additive"]["grid_worlds"]
    if not MAIN.is_dir():
        return files, [f"MISSING main directory: {MAIN}"] + [
            f"MISSING CELL G{i:02d}: expected {n}, found 0" for i, n in enumerate(counts, 1)
        ]
    for gi, expected in enumerate(counts, 1):
        actual = sorted(MAIN.glob(f"world_g{gi:02d}_n*.mat"))
        known = {f"world_g{gi:02d}_n{si:03d}.mat" for si in range(1, expected + 1)}
        actual_names = {p.name for p in actual}
        missing = [si for si in range(1, expected + 1) if f"world_g{gi:02d}_n{si:03d}.mat" not in actual_names]
        extras = sorted(actual_names - known)
        if len(actual) != expected or missing or extras:
            detail = f"MISSING CELL G{gi:02d}: expected {expected}, found {len(actual)}"
            if missing:
                preview = ",".join(str(x) for x in missing[:12])
                detail += f"; missing seed_index={preview}{'...' if len(missing) > 12 else ''}"
            if extras:
                detail += f"; unexpected={','.join(extras[:5])}"
            errors.append(detail)
        for si in range(1, expected + 1):
            p = MAIN / f"world_g{gi:02d}_n{si:03d}.mat"
            if p.is_file():
                files.append((gi, si, p))
    tmp = sorted(MAIN.glob("*.tmp*"))
    if tmp:
        errors.append(f"TEMPORARY checkpoint files present: {', '.join(p.name for p in tmp[:5])}")
    return files, errors


def main_manifest_errors(contract: dict[str, Any]) -> list[str]:
    """Completion marker is required in addition to the 986 independently checked files."""
    if contract.get("execution", {}).get("main_run_authorized") is not True:
        return ["current paper contract does not authorize the main-run lineage"]
    if not RUN_MANIFEST.is_file():
        return [f"MISSING final run manifest: {RUN_MANIFEST}"]
    try:
        manifest = json.loads(RUN_MANIFEST.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"invalid final run manifest: {exc}"]
    required = {"status": "complete", "workers": 12, "worlds": 986, "contract_sha256": sha256(PAPER)}
    errors = []
    for key, expected in required.items():
        if manifest.get(key) != expected:
            errors.append(f"run_manifest {key}={manifest.get(key)!r}, expected {expected!r}")
    return errors


def find_matlab() -> str | None:
    return shutil.which("matlab") or ("D:/MATLAB/bin/matlab.exe" if Path("D:/MATLAB/bin/matlab.exe").is_file() else None)


def matlab_metadata(files: Iterable[Path]) -> dict[str, dict[str, Any]]:
    """Decode MATLAB compact string/table metadata with MATLAB itself, never #refs# guesses."""
    matlab = find_matlab()
    if not matlab:
        raise AggregationError("MATLAB is required to verify compact MATLAB string/table metadata, but was not found")
    paths = list(files)
    rows: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="p0_2_metadata_") as td:
        temp = Path(td)
        # A MATLAB process is limited to 100 files per chunk.  Thus the 300-s
        # timeout is per bounded chunk rather than a false 300-s limit for all
        # 986 checkpoints; this keeps decoding recoverable on slow storage.
        for start in range(0, len(paths), METADATA_CHUNK_SIZE):
            chunk = paths[start:start + METADATA_CHUNK_SIZE]
            input_file = temp / f"paths_{start:04d}.txt"
            output_file = temp / f"metadata_{start:04d}.jsonl"
            helper = temp / f"extract_{start:04d}.m"
            input_file.write_text("\n".join(str(p.resolve()) for p in chunk) + "\n", encoding="utf-8")
            # Use MATLAB's public load/string/istable/jsonencode operations.  The helper is system-temporary.
            m_quote = lambda p: str(p).replace("\\", "/").replace("'", "''")
            helper.write_text(
                "infile='" + m_quote(input_file) + "'; outfile='" + m_quote(output_file) + "';\n"
                "paths=splitlines(string(fileread(infile))); paths=paths(strlength(paths)>0); fid=fopen(outfile,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));\n"
                "for ii=1:numel(paths), q=load(char(paths(ii)),'checkpoint'); c=q.checkpoint; assert(istable(c.subject_matrix),'subject_matrix must be table');\n"
                "m=struct('path',char(paths(ii)),'status',char(string(c.status)),'analysis_id',char(string(c.analysis_id)),"
                "'scenario',char(string(c.scenario)),'scenario_index',double(c.scenario_index),'grid_index',double(c.grid_index),"
                "'grid_id',char(string(c.grid_id)),'seed_index',double(c.seed_index),'seed',double(c.seed),'A_uV',double(c.A_uV),"
                "'beta_gen',double(c.beta_gen),'expected_slope',double(c.expected_slope),'truth_bridge_slope',double(c.truth_bridge_slope),"
                "'contract_sha256',char(string(c.contract_sha256)),'parent_historical_contract_sha256',char(string(c.parent_historical_contract_sha256)),"
                "'input_hashes',c.input_hashes,'subject_matrix_is_table',istable(c.subject_matrix),'subject_matrix_height',height(c.subject_matrix),"
                "'subject_matrix_variable_names',{cellstr(c.subject_matrix.Properties.VariableNames)},'truth_scenario',char(string(c.truth.scenario)),"
                "'started_at_utc',char(string(c.started_at_utc)),'completed_at_utc',char(string(c.completed_at_utc)));\n"
                "fprintf(fid,'%s\\n',jsonencode(m)); end; disp('P0_2_MATLAB_METADATA_OK');\n",
                encoding="utf-8",
            )
            command = [matlab, "-batch", f"run('{m_quote(helper)}')"]
            result = subprocess.run(command, text=True, encoding="utf-8", errors="replace", capture_output=True,
                                    timeout=METADATA_CHUNK_TIMEOUT_S)
            if result.returncode != 0 or not output_file.is_file():
                tail = (result.stdout + "\n" + result.stderr)[-3000:]
                raise AggregationError(f"MATLAB metadata chunk {start // METADATA_CHUNK_SIZE + 1} failed (exit {result.returncode}):\n{tail}")
            rows.extend(json.loads(line) for line in output_file.read_text(encoding="utf-8").splitlines() if line.strip())
    if len(rows) != len(paths):
        raise AggregationError(f"MATLAB metadata row count {len(rows)} != checkpoint count {len(paths)}")
    return {str(Path(row["path"]).resolve()): row for row in rows}


def validate_metadata(meta: dict[str, Any], contract: dict[str, Any], gi: int, si: int, path: Path,
                      parent_sha: str, required_contract_sha: str) -> None:
    expected_seed = int(contract["seed"]["paper_base_seed"] + 100000 * contract["seed"]["scenario_index"] + 1000 * gi + si)
    expected_slope = float(contract["additive"]["a_erp_grid_uV"][gi] ** 2 * contract["bridge"]["expected_slope_gain_G_log10_per_uV2_per_z"])
    expected = {
        "status": "complete", "analysis_id": contract["analysis_id"], "scenario": "additive",
        "scenario_index": 2, "grid_index": gi, "grid_id": contract["additive"]["grid_id"][gi - 1],
        "seed_index": si, "seed": expected_seed, "A_uV": float(contract["additive"]["a_erp_grid_uV"][gi]),
        "beta_gen": float(contract["additive"]["beta_gen_uV2_per_z"][gi - 1]), "expected_slope": expected_slope,
        "truth_bridge_slope": expected_slope, "contract_sha256": required_contract_sha,
        "parent_historical_contract_sha256": parent_sha, "truth_scenario": "additive",
    }
    for key, value in expected.items():
        if meta.get(key) != value:
            raise AggregationError(f"{path.name}: metadata {key}={meta.get(key)!r}, expected {value!r}")
    if not meta.get("subject_matrix_is_table") or int(meta.get("subject_matrix_height", 0)) <= 0:
        raise AggregationError(f"{path.name}: subject_matrix was not a nonempty MATLAB table")
    actual_hashes = meta.get("input_hashes", {})
    wanted_hashes = {key: value["sha256"] for key, value in contract["inputs"].items()}
    if actual_hashes != wanted_hashes:
        raise AggregationError(f"{path.name}: input hash manifest mismatch")
    for time_key in ("started_at_utc", "completed_at_utc"):
        try:
            datetime.fromisoformat(str(meta[time_key]).replace("Z", "+00:00"))
        except (KeyError, TypeError, ValueError) as exc:
            raise AggregationError(f"{path.name}: invalid checkpoint {time_key}") from exc


def read_world(gi: int, si: int, path: Path, meta: dict[str, Any], contract: dict[str, Any], parent_sha: str,
               required_contract_sha: str) -> list[dict[str, Any]]:
    validate_metadata(meta, contract, gi, si, path, parent_sha, required_contract_sha)
    rows: list[dict[str, Any]] = []
    subject_counts = dict(zip(contract["sign_flip"]["runs"], contract["sign_flip"]["subject_counts"], strict=True))
    with h5py.File(path, "r") as h:
        # HDF5 direct checks duplicate identity-bearing numeric fields from MATLAB load.
        if matlab_char(h, "checkpoint/grid_id") != meta["grid_id"]:
            raise AggregationError(f"{path.name}: HDF5 char grid_id disagrees with MATLAB load")
        for key in ("grid_index", "seed_index", "seed", "A_uV", "beta_gen", "expected_slope", "truth_bridge_slope"):
            if scalar(h, f"checkpoint/{key}") != float(meta[key]):
                raise AggregationError(f"{path.name}: HDF5 numeric {key} disagrees with MATLAB load")
        for run in contract["sign_flip"]["runs"]:
            observed, p_two, signs, stat = validate_signflip(h, run, subject_counts[run])
            coverage = coverage_from_arrays(signs, stat, observed, float(meta["truth_bridge_slope"]), contract)
            rows.append({
                "grid_id": meta["grid_id"], "grid_index": gi, "seed_index": si, "seed": int(meta["seed"]),
                "run": run, "A_uV": meta["A_uV"], "beta_gen": meta["beta_gen"],
                "true_bridge_slope": meta["truth_bridge_slope"], "observed_slope": observed,
                "p_two_sided_mpower": p_two, "reject_alpha_0_05": int(p_two <= ALPHA),
                **coverage,
                "robustfit_fit_count": scalar(h, "checkpoint/diagnostics/robustfit_fit_count"),
                "robustfit_iteration_limit_count": scalar(h, "checkpoint/diagnostics/robustfit_iteration_limit_count"),
                "elapsed_s": scalar(h, "checkpoint/timing/elapsed_s"), "peak_memory_bytes": scalar(h, "checkpoint/timing/peak_memory_bytes"),
                "checkpoint_filename": path.name, "checkpoint_sha256": sha256(path),
                "paper_contract_sha256": meta["contract_sha256"],
                "parent_historical_contract_sha256": meta["parent_historical_contract_sha256"],
            })
    return rows


def mean_sd(values: list[float]) -> tuple[float, float]:
    mean = sum(values) / len(values)
    sd = math.nan if len(values) < 2 else math.sqrt(sum((x - mean) ** 2 for x in values) / (len(values) - 1))
    return mean, sd


def summarize(rows: list[dict[str, Any]], contract: dict[str, Any]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(row["grid_id"], row["run"])].append(row)
    expected = dict(zip(contract["additive"]["grid_id"], contract["additive"]["grid_worlds"], strict=True))
    out: list[dict[str, Any]] = []
    for (grid_id, run), subset in sorted(grouped.items()):
        n, k = len(subset), sum(int(r["reject_alpha_0_05"]) for r in subset)
        coverage_k = sum(bool(r["truth_in_ci"]) for r in subset)
        slopes = [float(r["observed_slope"]) for r in subset]
        truth = float(subset[0]["true_bridge_slope"])
        sm, sd = mean_sd(slopes)
        lo, hi = cp_interval(k, n)
        clo, chi = cp_interval(coverage_k, n)
        total_fits = sum(float(r["robustfit_fit_count"]) for r in subset)
        iter_limits = sum(float(r["robustfit_iteration_limit_count"]) for r in subset)
        elapsed = [float(r["elapsed_s"]) for r in subset]
        memory = [float(r["peak_memory_bytes"]) for r in subset]
        out.append({
            "grid_id": grid_id, "run": run, "expected_worlds": expected[grid_id], "found_worlds": n,
            "A_uV": subset[0]["A_uV"], "beta_gen": subset[0]["beta_gen"], "true_bridge_slope": truth,
            "rejection_count": k, "power_hat": k / n, "power_cp95_low": lo, "power_cp95_high": hi, "power_mcse": mcse(k, n),
            "observed_slope_mean": sm, "observed_slope_sample_sd": sd, "bias": sm - truth,
            "rmse": math.sqrt(sum((x - truth) ** 2 for x in slopes) / n),
            "coverage_count": coverage_k, "coverage_rate": coverage_k / n, "coverage_cp95_low": clo,
            "coverage_cp95_high": chi, "coverage_mcse": mcse(coverage_k, n),
            "robustfit_iteration_limit_total": iter_limits, "robustfit_fit_total": total_fits,
            "robustfit_iteration_limit_rate": iter_limits / total_fits if total_fits else math.nan,
            "elapsed_s_mean": sum(elapsed) / n, "elapsed_s_max": max(elapsed),
            "peak_memory_bytes_mean": sum(memory) / n, "peak_memory_bytes_max": max(memory),
            "paper_contract_sha256": subset[0]["paper_contract_sha256"],
        })
    return out


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def null_reference(contract: dict[str, Any]) -> list[dict[str, Any]]:
    complete = json.loads((NULL / "nullpool_complete.json").read_text(encoding="utf-8"))
    if complete.get("status") != "complete" or complete.get("n_worlds") != contract["null_reference"]["worlds"]:
        raise AggregationError("historical A=0 null pool is incomplete")
    if complete.get("contract_sha256") != contract["null_reference"]["contract_sha256"]:
        raise AggregationError("historical A=0 null pool contract SHA mismatch")
    fpr = {(r["run"], r["test"]): r for r in read_tsv(NULL / "nullpool_fpr.tsv")}
    dist = {(r["run"], r["metric"]): r for r in read_tsv(NULL / "nullpool_distribution_summary.tsv")}
    rows: list[dict[str, Any]] = []
    for run in contract["sign_flip"]["runs"]:
        old, distribution = fpr.get((run, "Mpower")), dist.get((run, "Mpower"))
        if old is None or distribution is None:
            raise AggregationError(f"historical A=0 Mpower reference missing for {run}")
        n, k = int(old["n_worlds"]), int(old["trigger_count_alpha_0_05"])
        lo, hi = cp_interval(k, n)
        rows.append({
            "reference_type": "historical_A0_reference_not_grid_cell", "run": run, "metric": "Mpower",
            "n_worlds": n, "rejections_alpha_0_05": k, "fpr": k / n, "cp95_low_recomputed": lo,
            "cp95_high_recomputed": hi, "mcse": mcse(k, n), "null_mean": float(distribution["mean"]),
            "null_sd": float(distribution["std_sigma_null"]), "historical_contract_sha256": complete["contract_sha256"],
        })
    return rows


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise AggregationError(f"refusing to write empty table {path.name}")
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        writer.writeheader(); writer.writerows(rows)


def markdown(summary_rows: list[dict[str, Any]], null_rows: list[dict[str, Any]], contract_sha: str) -> str:
    lines = ["# P0-2 SIM power-curve aggregation", "", "Read-only aggregation of the completed frozen grid. Mpower is the sole P0-2 outcome; Mphase and p_maxstat are excluded from power/coverage conclusions.", "",
             "## Mpower grid summary", "", "|grid|run|n|rejections|power|CP95%|MCSE|truth slope|mean slope|bias|RMSE|coverage|coverage CP95%|", "|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|"]
    for r in summary_rows:
        lines.append(f"|{r['grid_id']}|{r['run']}|{r['found_worlds']}|{r['rejection_count']}|{r['power_hat']:.8f}|[{r['power_cp95_low']:.8f}, {r['power_cp95_high']:.8f}]|{r['power_mcse']:.8f}|{r['true_bridge_slope']:.10g}|{r['observed_slope_mean']:.10g}|{r['bias']:.10g}|{r['rmse']:.10g}|{r['coverage_rate']:.8f}|[{r['coverage_cp95_low']:.8f}, {r['coverage_cp95_high']:.8f}]|")
    lines += ["", "## Historical A=0 reference (not a grid cell)", "", "|run|n|rejections|FPR|CP95%|MCSE|null mean|null SD|historical contract SHA|", "|---|---:|---:|---:|---|---:|---:|---:|---|"]
    for r in null_rows:
        lines.append(f"|{r['run']}|{r['n_worlds']}|{r['rejections_alpha_0_05']}|{r['fpr']:.8f}|[{r['cp95_low_recomputed']:.8f}, {r['cp95_high_recomputed']:.8f}]|{r['mcse']:.8f}|{r['null_mean']:.10g}|{r['null_sd']:.10g}|{r['historical_contract_sha256']}|")
    lines += ["", "## Fixed calculations", "", "- Clopper--Pearson endpoints are exact beta inversions with the k=0/k=n endpoint rules; MCSE is `sqrt(p_hat*(1-p_hat)/n)`.", "- Each coverage interval uses the frozen exact Mpower sign-flip inversion on theta = -0.02:1e-5:0.02, accepting only `p(theta) > 0.05`. No bootstrap, normal, or percentile interval is substituted.", f"- Paper contract SHA-256: `{contract_sha}`.", ""]
    return "\n".join(lines)


def checkpoint_time_range(metadata: dict[str, dict[str, Any]]) -> tuple[str, str]:
    def parse(value: str) -> datetime:
        return datetime.fromisoformat(value.replace("Z", "+00:00")).astimezone(UTC)
    started = min(parse(str(row["started_at_utc"])) for row in metadata.values())
    completed = max(parse(str(row["completed_at_utc"])) for row in metadata.values())
    return (started.isoformat(timespec="milliseconds").replace("+00:00", "Z"),
            completed.isoformat(timespec="milliseconds").replace("+00:00", "Z"))


def atomic_success_outputs(worlds: list[dict[str, Any]], summaries: list[dict[str, Any]], nulls: list[dict[str, Any]],
                           checkpoints: list[Path], metadata: dict[str, dict[str, Any]], contract: dict[str, Any]) -> None:
    created = now_utc()
    # Stage next to OUT so the final os.replace stays on one filesystem
    # (Windows refuses cross-drive renames, e.g. C: temp -> E: tables).
    with tempfile.TemporaryDirectory(prefix="p0_2_outputs_", dir=OUT.parent) as td:
        temp = Path(td)
        write_tsv(temp / "p0_2_powercurve_worlds.tsv", worlds)
        write_tsv(temp / "p0_2_powercurve_summary.tsv", summaries)
        write_tsv(temp / "p0_2_null_reference.tsv", nulls)
        (temp / "p0_2_powercurve_summary.md").write_text(markdown(summaries, nulls, sha256(PAPER)), encoding="utf-8")
        input_started, input_completed = checkpoint_time_range(metadata)
        manifest = {
            "status": "complete", "created_at_utc": created, "script_sha256": sha256(Path(__file__)),
            "paper_contract_sha256": sha256(PAPER), "parent_historical_contract_sha256": sha256(PARENT),
            "input_started_at_utc": input_started,
            "input_finished_at_utc": input_completed,
            "python_version": sys.version, "h5py_version": h5py.__version__, "scipy_version": scipy.__version__,
            "completeness": {"expected_worlds": sum(contract["additive"]["grid_worlds"]), "validated_worlds": len(checkpoints), "validated_world_run_rows": len(worlds), "errors": 0},
            "checkpoints": [{"filename": p.name, "sha256": sha256(p)} for p in checkpoints],
        }
        (temp / "p0_2_aggregation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        OUT.mkdir(parents=True, exist_ok=True)
        for p in temp.iterdir():
            os.replace(p, OUT / p.name)


def reader_smoke(contract: dict[str, Any]) -> int:
    files = [V2 / f"world_g{gi:02d}_n001.mat" for gi in (1, 6, 7)]
    missing = [str(p) for p in files if not p.is_file()]
    if missing:
        raise AggregationError("reader smoke checkpoint(s) missing: " + ", ".join(missing))
    metadata = matlab_metadata(files)
    parent_sha = sha256(PARENT)
    result = []
    for gi, path in zip((1, 6, 7), files, strict=True):
        # read_world also tests HDF5 strings/numerics, all-plus direction, exact p, and inversion.
        rows = read_world(gi, 1, path, metadata[str(path.resolve())], contract, parent_sha,
                          HISTORICAL_V2_PROFILING_CONTRACT_SHA)
        result.append({"file": path.name, "grid_id": f"G{gi:02d}", "runs_verified": [r["run"] for r in rows], "result": "profiling reader smoke, not a paper result"})
    print(json.dumps({"status": "PASS", "result": "profiling reader smoke, not a paper result", "checkpoints": result}, indent=2))
    return 0


def self_test(contract: dict[str, Any]) -> int:
    # Binary-exact values avoid artificial round-off ties in the independently
    # evaluated direct dot product while still exercising nonzero translations.
    signs = np.array([[1, 1, 1, 1, -1, -1, -1, -1],
                      [1, 1, -1, -1, 1, 1, -1, -1],
                      [1, -1, 1, -1, 1, -1, 1, -1]], dtype=float)
    x = np.array([0.125, -0.0625, 0.03125])
    stats = signs.T @ x / x.size
    observed = float(stats[0])
    grid = np.array([-0.0625, -0.03125, 0.0, 0.03125, 0.0625])
    fast = p_curve_from_arrays(signs, stats, observed, grid)
    direct = np.array([np.mean(np.abs(signs.T @ (x - theta) / x.size) >= abs(observed - theta)) for theta in grid])
    if not np.array_equal(fast, direct):
        raise AssertionError("translation formula disagrees with direct sign enumeration")
    if cp_interval(0, 8)[0] != 0.0 or cp_interval(8, 8)[1] != 1.0:
        raise AssertionError("Clopper-Pearson endpoint rule failed")
    if not math.isfinite(cp_interval(1, 8)[0]) or not math.isfinite(cp_interval(7, 8)[1]):
        raise AssertionError("Clopper-Pearson finite endpoint rule failed")
    print("PASS self-test: translation=direct enumeration; CP endpoints correct; incomplete-check logic is exercised by default mode.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reader-smoke", action="store_true", help="validate only the three v2 profiling checkpoints")
    parser.add_argument("--self-test", action="store_true", help="run synthetic translation and CP tests only")
    args = parser.parse_args()
    if args.reader_smoke and args.self_test:
        parser.error("--reader-smoke and --self-test are mutually exclusive")
    contract = json.loads(PAPER.read_text(encoding="utf-8"))
    if args.self_test:
        return self_test(contract)
    if args.reader_smoke:
        return reader_smoke(contract)
    files, errors = list_main_files(contract)
    errors.extend(main_manifest_errors(contract))
    if errors:
        print("P0-2 aggregation refused: main run is incomplete or has invalid file cardinality.", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 2
    if len(files) != sum(contract["additive"]["grid_worlds"]):
        print("P0-2 aggregation refused: total checkpoint count is not 986.", file=sys.stderr)
        return 2
    try:
        paths = [p for _, _, p in files]
        metadata = matlab_metadata(paths)
        parent_sha = sha256(PARENT)
        worlds: list[dict[str, Any]] = []
        for gi, si, path in files:
            worlds.extend(read_world(gi, si, path, metadata[str(path.resolve())], contract, parent_sha, sha256(PAPER)))
        summaries = summarize(worlds, contract)
        if len(summaries) != len(contract["additive"]["grid_id"]) * len(contract["sign_flip"]["runs"]):
            raise AggregationError("incomplete grid/run summary after validation")
        nulls = null_reference(contract)
        atomic_success_outputs(worlds, summaries, nulls, paths, metadata, contract)
        print(f"P0-2 aggregation complete: {len(paths)} worlds, {len(worlds)} world/run rows; output={OUT}")
        return 0
    # h5py failures surface as OSError subclasses; there is no h5py.H5Error.
    except (AggregationError, OSError, ValueError, KeyError) as exc:
        print(f"P0-2 aggregation refused: {exc}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    sys.exit(main())
