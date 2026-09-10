# Release port of the analysis-workspace tool of the same name; see README in this folder.
"""Independent Stage-C validator; reads bridge TSVs plus the read-only Stage B reference."""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import beta, skew

def _repo_root():
    for _q in Path(__file__).resolve().parents:
        if (_q / "MANIFEST.md").is_file():
            return _q
    raise RuntimeError("release repository root not found")
ROOT = _repo_root()
sys.path.insert(0, str(ROOT / "scripts" / "ieeg" / "paper"))
from verify_p0_6_paired_contrast import newcombe_method10  # noqa: E402


def cp(k, n):
    return (0.0 if k == 0 else float(beta.ppf(.025, k, n - k + 1)),
            1.0 if k == n else float(beta.ppf(.975, k + 1, n - k)))


def signflip(x):
    n = len(x)
    if not n:
        return math.nan, math.nan, math.nan
    rows = np.arange(2 ** n, dtype=np.uint32)[:, None]
    bits = (rows >> np.arange(n - 1, -1, -1, dtype=np.uint32)) & 1
    signs = 1.0 - 2.0 * bits
    signs[0, :] = 1.0
    stats = signs @ x / n
    return float(stats[0]), float(np.mean(np.abs(stats) >= abs(stats[0]))), float(np.std(stats, ddof=1))


def close(a, b, tol=1e-12):
    return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= tol


def validate(output: Path, mode: str):
    pairs = pd.read_csv(output / "validator_pair_rows.tsv", sep="\t")
    subs = pd.read_csv(output / "validator_subject_rows.tsv", sep="\t")
    worlds = pd.read_csv(output / "validator_world_rows.tsv", sep="\t")
    ref_worlds = pd.read_csv(ROOT / "tools" / "stageb_mphase_power" / "output" / "validator_world_rows.tsv", sep="\t")
    ref_pairs = pd.read_csv(ROOT / "tools" / "stageb_mphase_power" / "output" / "validator_pair_rows.tsv", sep="\t")
    errors = []
    flip_stats = {}  # per (reducer, run): one-flip budget and straddle counts
    contract = json.loads((Path(__file__).resolve().parent / "stagec_mphase_gainflat_contract.json").read_text(encoding="utf-8"))
    eps_cfg = contract["injection"]["epsilon_resultant_check"]
    expected = 2 if mode == "smoke" else 60
    seeds = sorted(pairs.seed_index.unique().tolist())
    if seeds != list(range(1, expected + 1)):
        errors.append(f"seed set is not 1:{expected}")
    if pairs.world_seed.min() < 22360828 or pairs.world_seed.max() > 22360887:
        errors.append("world seed outside shared Stage-B namespace")
    if set(pairs.arm) != {"P0F"} or not (pairs.mphase_surrogate_reducer == "mean").all() or not (pairs.refit_surrogate_count == 32).all():
        errors.append("arm/reducer/refit identity failure")
    if not (pairs.gain_slope_b == 0).all() or not (worlds.gain_slope_b == 0).all():
        errors.append("gain slope is not flattened to 0 everywhere")
    # CRN pairing proof: base fingerprints must equal Stage B's per (seed, subject,
    # session, pair). Smoke uses synthetic pairs with no Stage B counterpart; skip.
    if mode == "formal":
        rb = ref_pairs[ref_pairs.arm == "P0"][["seed_index", "subject", "session_id", "pair_id", "base_fingerprint"]]
        j = pairs.merge(rb, on=["seed_index", "subject", "session_id", "pair_id"], how="left",
                        suffixes=("", "_stageb"), validate="one_to_one")
        if j.base_fingerprint_stageb.isna().any():
            errors.append("pair rows without a Stage B fingerprint counterpart")
        elif (j.base_fingerprint != j.base_fingerprint_stageb).any():
            errors.append("base fingerprint differs from Stage B (CRN pairing broken)")
    # Injection / stream provenance checks.
    if pairs.epsilon_resultant.isna().any():
        errors.append("epsilon_resultant has NaN rows")
    big = pairs.epsilon_count >= 1000
    if mode == "formal" and not big.all():
        errors.append("formal run has pairs with <1000 injected events")
    pb = pairs[big]
    if (pb.epsilon_resultant > float(eps_cfg["P0_pair_max_abs"])).any():
        errors.append("P0F per-pair epsilon resultant above corruption cap")
    wm = pb.groupby("seed_index").epsilon_resultant.mean()
    if (wm > float(eps_cfg["P0_world_mean_max"])).any():
        errors.append("P0F world-mean epsilon resultant above uniform tolerance")
    for seed in seeds:
        pw = pairs[pairs.seed_index == seed]
        dup = pw.groupby(["subject", "session_id", "pair_id"]).injection_stream_seed.nunique()
        if (dup != 1).any() or pw.injection_stream_seed.duplicated().any():
            errors.append(f"seed {seed}: injection stream seed not unique per pair")
    if (worlds.truth_a_g06_uV != 4.61799281479342).any() or (worlds.truth_grid_id != "G06").any():
        errors.append("injection truth amplitude/grid_id mismatch")
    if (pairs.n_injected_phase_intersection != pairs.n_phase_events).any() or (pairs.n_injected_events < pairs.n_phase_events).any():
        errors.append("injected/phase event count identity violated")
    # Independent reducer + sign-flip recomputation.
    rows = []
    for seed in seeds:
        for run in ("R1", "R2"):
            pr = pairs[(pairs.seed_index == seed) & (pairs.run == run)]
            for reducer in ("median", "mean"):
                vals = []
                for subject, g in pr.groupby("subject", sort=True):
                    v = g.corrected_mphase.to_numpy(float)
                    v = v[np.isfinite(v)]
                    if not len(v):
                        continue
                    stat = float(np.median(v) if reducer == "median" else np.mean(v))
                    ps = float(skew(v, bias=False)) if len(v) >= 3 and not np.all(v == v[0]) else math.nan
                    sr = subs[(subs.seed_index == seed) & (subs.run == run) & (subs.reducer == reducer) & (subs.subject == subject)]
                    if len(sr) != 1:
                        errors.append(f"{seed}/{run}/{reducer}/{subject}: subject row")
                    else:
                        r = sr.iloc[0]
                        if not close(stat, float(r.Mphase)) or int(r.n_pairs_finite) != len(v) or not close(ps, float(r.pair_corrected_skewness)):
                            errors.append(f"{seed}/{run}/{reducer}/{subject}: reducer mismatch")
                    vals.append(stat)
                obs, p, sd = signflip(np.asarray(vals, float))
                wr = worlds[(worlds.seed_index == seed) & (worlds.run == run) & (worlds.reducer == reducer)]
                if len(wr) != 1:
                    errors.append(f"{seed}/{run}/{reducer}: world row")
                else:
                    r = wr.iloc[0]
                    fs = flip_stats.setdefault((reducer, run), {"budget": 0.0, "s005": 0, "s050": 0, "s075": 0})
                    pcv = float(r.p_two_sided)
                    fs["budget"] += abs(p - pcv)
                    fs["s005"] += int((p <= .05) != (pcv <= .05))
                    fs["s050"] += int((p < .5) != (pcv < .5))
                    fs["s075"] += int((p > .75) != (pcv > .75))
                    if not close(obs, float(r.observed)) or abs(p - pcv) > 1 / (2 ** max(len(vals), 1)) + 1e-12 or not close(sd, float(r.studentization_sd)):
                        errors.append(f"{seed}/{run}/{reducer}: signflip mismatch")
                rows.append(dict(seed_index=seed, arm="P0F", reducer=reducer, run=run,
                                 observed=obs, p_two_sided=p, n_subjects=len(vals)))
    rw = pd.DataFrame(rows)
    fpr, dist = [], []
    for reducer in ("median", "mean"):
        for run in ("R1", "R2"):
            g = rw[(rw.reducer == reducer) & (rw.run == run)]
            n = len(g); k = int((g.p_two_sided <= .05).sum()); lo, hi = cp(k, n) if n else (math.nan, math.nan)
            fpr.append(dict(arm="P0F", reducer=reducer, run=run, n_worlds=n, rejections=k,
                            fpr=k / n if n else math.nan, cp_low=lo, cp_high=hi))
            below = int((g.p_two_sided < .5).sum()); above = int((g.p_two_sided > .75).sum())
            bl, bh = cp(below, n) if n else (math.nan, math.nan)
            al, ah = cp(above, n) if n else (math.nan, math.nan)
            dist.append(dict(arm="P0F", reducer=reducer, run=run, n_worlds=n,
                             mean_p=float(g.p_two_sided.mean()) if n else math.nan,
                             p_lt_0p5_count=below, p_lt_0p5_fraction=below / n if n else math.nan,
                             p_lt_0p5_cp_low=bl, p_lt_0p5_cp_high=bh, p_gt_0p75_count=above,
                             p_gt_0p75_fraction=above / n if n else math.nan,
                             p_gt_0p75_cp_low=al, p_gt_0p75_cp_high=ah))
    paired = []
    for reducer in ("mean", "median"):
        for run in ("R1", "R2"):
            x = rw[(rw.reducer == reducer) & (rw.run == run)][["seed_index", "p_two_sided"]]
            y = ref_worlds[(ref_worlds.arm == "P0") & (ref_worlds.reducer == reducer) & (ref_worlds.run == run)][["seed_index", "p_two_sided"]].rename(columns={"p_two_sided": "p_ref"})
            z = x.merge(y, on="seed_index", validate="one_to_one")
            n = len(z); pf = z.p_two_sided <= .05; pbv = z.p_ref <= .05
            both = int((pf & pbv).sum()); flat_only = int((pf & ~pbv).sum())
            ref_only = int((~pf & pbv).sum()); neither = int((~pf & ~pbv).sum())
            d, lo, hi = newcombe_method10(flat_only, ref_only, both, n) if n else (math.nan,) * 3
            paired.append(dict(reducer=reducer, run=run, n=n, both_reject=both,
                               p0f_only_reject=flat_only, p0stageb_only_reject=ref_only, neither_reject=neither,
                               fpr_p0f=(both + flat_only) / n if n else math.nan,
                               fpr_p0_stageb=(both + ref_only) / n if n else math.nan,
                               delta_p0f_minus_stageb=d, ci_low=lo, ci_high=hi))
    # Summary agreement with per-group flip-budget tolerances.
    for name, frame in (("main_fpr.tsv", pd.DataFrame(fpr)), ("main_p_distribution.tsv", pd.DataFrame(dist)), ("main_paired_contrast.tsv", pd.DataFrame(paired))):
        got = pd.read_csv(output / name, sep="\t")
        if list(got.columns) != list(frame.columns) or got.shape != frame.shape:
            errors.append(f"{name}: schema mismatch")
            continue
        keycols = ["arm", "reducer", "run"] if "arm" in frame.columns else ["reducer", "run"]
        for _, frow in frame.iterrows():
            sel = got
            for kc in keycols:
                sel = sel[sel[kc] == frow[kc]]
            if len(sel) != 1:
                errors.append(f"{name}: row key mismatch for {tuple(frow[kc] for kc in keycols)}")
                continue
            mrow = sel.iloc[0]
            if "arm" in keycols:
                fs = flip_stats[(frow["reducer"], frow["run"])]
                n_g = max(int(frow["n_worlds"]), 1)
                atol_map = {"mean_p": fs["budget"] / n_g + 1e-12, "rejections": float(fs["s005"]),
                            "p_lt_0p5_count": float(fs["s050"]), "p_gt_0p75_count": float(fs["s075"]),
                            "fpr": fs["s005"] / n_g + 1e-12,
                            "p_lt_0p5_fraction": fs["s050"] / n_g + 1e-12,
                            "p_gt_0p75_fraction": fs["s075"] / n_g + 1e-12}
                cp_src = {"cp_low": ("rejections", 0), "cp_high": ("rejections", 1),
                          "p_lt_0p5_cp_low": ("p_lt_0p5_count", 0), "p_lt_0p5_cp_high": ("p_lt_0p5_count", 1),
                          "p_gt_0p75_cp_low": ("p_gt_0p75_count", 0), "p_gt_0p75_cp_high": ("p_gt_0p75_count", 1)}
            else:
                fs = flip_stats[(frow["reducer"], frow["run"])]
                n_g = max(int(frow["n"]), 1)
                # the reference arm's own one-flip straddles are not tracked here; allow one
                # straddle per side as the documented cross-stage tolerance
                atol_map = {"both_reject": float(fs["s005"] + 1), "p0f_only_reject": float(fs["s005"] + 1),
                            "p0stageb_only_reject": float(fs["s005"] + 1), "neither_reject": float(fs["s005"] + 1),
                            "fpr_p0f": fs["s005"] / n_g + 1e-12, "fpr_p0_stageb": (fs["s005"] + 1) / n_g + 1e-12}
                cp_src = {}
            for col in frame.columns:
                if col in keycols:
                    continue
                if col in cp_src:
                    ccol, cidx = cp_src[col]
                    okc = close(float(mrow[col]), float(cp(int(mrow[ccol]), int(mrow["n_worlds"]))[cidx]))
                elif col in ("delta_p0f_minus_stageb", "ci_low", "ci_high"):
                    d2, l2, h2 = newcombe_method10(int(mrow["p0f_only_reject"]), int(mrow["p0stageb_only_reject"]),
                                                  int(mrow["both_reject"]), int(mrow["n"]))
                    okc = close(float(mrow[col]), {"delta_p0f_minus_stageb": d2, "ci_low": l2, "ci_high": h2}[col])
                elif pd.api.types.is_numeric_dtype(frame[col]):
                    okc = close(float(mrow[col]), float(frow[col]), atol_map.get(col, 1e-12))
                else:
                    okc = str(mrow[col]) == str(frow[col])
                if not okc:
                    errors.append(f"{name}: mismatch in {col} for {tuple(frow[kc] for kc in keycols)}")
    d = json.loads((output / "main_decision.json").read_text(encoding="utf-8"))
    if mode == "formal":
        complete = seeds == list(range(1, 61)) and all(
            int(rw[rw.run == run].n_subjects.min()) == (16 if run == "R1" else 13) for run in ("R1", "R2"))
        mf = pd.DataFrame(fpr); mp = pd.DataFrame(paired)
        mf = mf[mf.reducer == "mean"]; mp = mp[mp.reducer == "mean"]
        c1 = bool(all(r.ci_high < 0 for r in mp.itertuples()))
        c2 = bool(all(r.cp_low <= .05 <= r.cp_high for r in mf.itertuples()))
        if c1 and c2:
            mech = "gain_coupling_dominant"
        elif c1 and not c2:
            mech = "both_contribute"
        elif (not c1) and (not c2):
            mech = "envelope_or_other_dominant"
        else:
            mech = "indeterminate"
        expected_status = "COMPLETE" if complete else "INCOMPLETE"
        if d.get("status") != expected_status or d.get("rules", {}).get("C1_gain_coupling_contribution") != c1 \
                or d.get("rules", {}).get("C2_residual_elevation_absent") != c2 \
                or (complete and d.get("mechanism_class") != mech):
            errors.append("decision mismatch")
    elif d.get("status") != "INCOMPLETE_SMOKE" or d.get("scientific_rules_evaluated") is not False:
        errors.append("smoke decision improperly evaluates scientific rules")
    report = {"status": "FAIL" if errors else "PASS", "mode": mode, "checkpoints": len(seeds),
              "errors": errors,
              "checks": ["CRN base-fingerprint pairing vs Stage B", "epsilon uniformity (world-mean)",
                         "injection/stream provenance", "injection truth amplitude+grid_id+gain",
                         "event count identities", "pair reducers", "exact sign-flip column 2",
                         "FPR/CP", "paired Newcombe vs Stage B", "summary agreement", "decision"]}
    (output / "validator_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("output", type=Path)
    p.add_argument("--mode", choices=("smoke", "formal"), required=True)
    a = p.parse_args()
    r = validate(a.output, a.mode)
    print(json.dumps(r, indent=2))
    raise SystemExit(0 if r["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
