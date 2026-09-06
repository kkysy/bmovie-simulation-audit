#!/usr/bin/env python3
"""Read-only P0-1 calibration summary from completed ACC00-SIM pools."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import h5py
from scipy.special import betaincinv


ROOT = Path(__file__).resolve().parents[3]
NULL = ROOT / "processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/nullpool_edgeguard_20260829"
REC = ROOT / "processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/recovery_mpower_only_20260829"
OUT = ROOT / "processed/subject/group/ieeg_methods_paper/p0_1_calibration_summary"
ALLOWED_POOL_SHAS = {
    "E0930B71CE91CD1641EDD9B1F7A1012C7DBB7769458ACB96386CA5892EAF8CEF",
    "CF9C9B989644C346544AAE5F982F3644C2718A66E645CB9C196F48B8D89690BE",
}
ALPHA = 0.05


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def cp_interval(k: int, n: int, alpha: float = ALPHA) -> tuple[float, float]:
    """Clopper-Pearson: [I^-1(alpha/2;k,n-k+1), I^-1(1-alpha/2;k+1,n-k)]."""
    lo = 0.0 if k == 0 else float(betaincinv(k, n - k + 1, alpha / 2.0))
    hi = 1.0 if k == n else float(betaincinv(k + 1, n - k, 1.0 - alpha / 2.0))
    return lo, hi


def mean_sd(values: list[float]) -> tuple[float, float]:
    m = sum(values) / len(values)
    s = math.sqrt(sum((x - m) ** 2 for x in values) / (len(values) - 1))
    return m, s


def read_scalar(h5: h5py.File, key: str) -> float:
    return float(h5[key][()][0][0] if h5[key].shape == (1, 1) else h5[key][()])


def main() -> None:
    if "pre_mphase_repair" in str(NULL).lower():
        raise RuntimeError("refusing historical pre-Mphase-repair pool")
    complete = json.loads((NULL / "nullpool_complete.json").read_text(encoding="utf-8"))
    null_sha = str(complete.get("contract_sha256", ""))
    if complete.get("n_worlds") != 400 or null_sha not in ALLOWED_POOL_SHAS:
        raise RuntimeError("null pool count or contract SHA mismatch")
    rec_complete = json.loads((REC / "recovery_complete.json").read_text(encoding="utf-8"))
    rec_sha = str(rec_complete.get("contract_sha256", ""))
    if rec_complete.get("n_worlds") != 300 or rec_sha not in ALLOWED_POOL_SHAS:
        raise RuntimeError("additive pool is not 300 worlds")

    fpr = read_tsv(NULL / "nullpool_fpr.tsv")
    dist = {(r["run"], r["metric"]): r for r in read_tsv(NULL / "nullpool_distribution_summary.tsv")}
    rows: list[dict[str, str | int | float]] = []
    for r in fpr:
        run, family = r["run"], r["test"]
        n, k = int(r["n_worlds"]), int(r["trigger_count_alpha_0_05"])
        lo, hi = cp_interval(k, n)
        if family == "maxstat_family":
            vals: list[float] = []
            for i in range(1, n + 1):
                with h5py.File(NULL / f"world_{i:04d}.mat", "r") as h:
                    vals.append(read_scalar(h, f"checkpoint/inference/{run}/max_observed"))
            mean, sd = mean_sd(vals)
            metric = "family_max_observed"
        else:
            d = dist[(run, family)]
            mean, sd, metric = float(d["mean"]), float(d["std_sigma_null"]), family
        rows.append({
            "run": run, "estimator_family": metric, "n_worlds": n,
            "rejections_alpha_0_05": k, "fpr": k / n, "ci95_low": lo, "ci95_high": hi,
            "mcse": math.sqrt((k / n) * (1 - k / n) / n), "null_mean": mean,
            "null_sd": sd, "null_mean_over_sd": mean / sd if sd else math.nan,
            "seed_count": n, "contract_sha256": null_sha,
        })

    additive: list[dict[str, str | int | float]] = []
    rec_rows = read_tsv(REC / "recovery_worlds.tsv")
    for effect in sorted({int(r["effect_index"]) for r in rec_rows}):
        subset = [r for r in rec_rows if int(r["effect_index"]) == effect]
        if len(subset) != 100:
            raise RuntimeError(f"effect {effect} has {len(subset)} worlds, expected 100")
        amp = float(subset[0]["a_erp_uV"])
        for run in ("R1", "R2"):
            pcol = f"{run}_Mpower_p_exact"
            scol = f"{run}_Mpower_slope_log_residual"
            k = sum(float(r[pcol]) <= ALPHA for r in subset)
            lo, hi = cp_interval(k, len(subset))
            slopes = [float(r[scol]) for r in subset]
            sm, ss = mean_sd(slopes)
            additive.append({
                "effect_index": effect, "a_erp_uV": amp, "run": run, "n_worlds": 100,
                "rejections_alpha_0_05": k, "rejection_rate": k / 100,
                "ci95_low": lo, "ci95_high": hi,
                "mcse": math.sqrt((k / 100) * (1 - k / 100) / 100),
                "slope_mean_log10_per_z": sm, "slope_sd_log10_per_z": ss,
                "contract_sha256": rec_sha,
            })

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "p0_1_null_calibration.tsv").open("w", newline="", encoding="utf-8") as f:
        fields = list(rows[0].keys()); w = csv.DictWriter(f, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(rows)
    with (OUT / "p0_1_additive_recovery.tsv").open("w", newline="", encoding="utf-8") as f:
        fields = list(additive[0].keys()); w = csv.DictWriter(f, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(additive)

    md = ["# P0-1 calibration summary", "", "Read-only summary of completed pools; this is a paper result table, not a gate.", "",
          "## Null calibration", "", "|run|family|n|rejections|FPR|CP95% CI|MCSE|null mean|null SD|mean/SD|contract SHA|", "|---|---|---:|---:|---:|---|---:|---:|---:|---:|---|"]
    for r in rows:
        md.append(f"|{r['run']}|{r['estimator_family']}|{r['n_worlds']}|{r['rejections_alpha_0_05']}|{r['fpr']:.8f}|[{r['ci95_low']:.8f}, {r['ci95_high']:.8f}]|{r['mcse']:.8f}|{r['null_mean']:.10g}|{r['null_sd']:.10g}|{r['null_mean_over_sd']:.6f}|{r['contract_sha256']}|")
    md += ["", "CP interval uses the exact beta inversion `betaincinv(k, n-k+1, alpha/2)` and `betaincinv(k+1, n-k, 1-alpha/2)`; MCSE is `sqrt(p_hat*(1-p_hat)/n)`. Family max-stat null mean/SD are computed from `checkpoint/inference/{R1,R2}/max_observed` in all 400 MAT checkpoints.", "", "## Additive recovery (Mpower only)", "", "|effect|A_ERP (uV)|run|n|rejections|rate|CP95% CI|MCSE|slope mean|slope SD|contract SHA|", "|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---|"]
    for r in additive:
        md.append(f"|{r['effect_index']}|{r['a_erp_uV']:.15g}|{r['run']}|{r['n_worlds']}|{r['rejections_alpha_0_05']}|{r['rejection_rate']:.8f}|[{r['ci95_low']:.8f}, {r['ci95_high']:.8f}]|{r['mcse']:.8f}|{r['slope_mean_log10_per_z']:.10g}|{r['slope_sd_log10_per_z']:.10g}|{r['contract_sha256']}|")
    md += ["", "Mphase fields are descriptive only and are intentionally excluded because Mphase calibration failed.", "", "## Lineage", "", f"- Null pool: `{NULL}`; `nullpool_fpr.tsv`, `nullpool_distribution_summary.tsv`, 400 `world_*.mat`; contract SHA `{null_sha}`.", f"- Additive pool: `{REC}`; `recovery_worlds.tsv`, `recovery_complete.json`; contract SHA `{rec_sha}`.", "- Excluded by design: `nullpool_stage2_pre_mphase_repair_20260828` and all pre-repair historical pools.", "- Source implementation evidence: `scripts/ieeg/run_acc00_nullpool.m:84-116` (aggregation fields and family p-values); `scripts/ieeg/acc00_sim_contract.json:additive`, `:world_plan`, and `:implementation_boundary` (frozen parameters and worker rationale)."]
    (OUT / "p0_1_calibration_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
