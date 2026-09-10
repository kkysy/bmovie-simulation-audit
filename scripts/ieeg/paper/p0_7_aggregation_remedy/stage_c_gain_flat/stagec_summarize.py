# Release port of the analysis-workspace tool of the same name; see README in this folder.
"""Stage C summary: P0F FPR/dist + paired Newcombe vs the Stage B P0 reference arm."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import pandas as pd
from scipy.stats import beta

def _repo_root():
    for _q in Path(__file__).resolve().parents:
        if (_q / "MANIFEST.md").is_file():
            return _q
    raise RuntimeError("release repository root not found")
ROOT = _repo_root()
sys.path.insert(0, str(ROOT / "scripts" / "ieeg" / "paper"))
from verify_p0_6_paired_contrast import newcombe_method10  # noqa: E402


def cp(k: int, n: int):
    return (0.0 if k == 0 else float(beta.ppf(.025, k, n - k + 1)),
            1.0 if k == n else float(beta.ppf(.975, k + 1, n - k)))


def summarize(output: Path, mode: str) -> dict:
    w = pd.read_csv(output / "validator_world_rows.tsv", sep="\t")
    ref = pd.read_csv(ROOT / "tools" / "stageb_mphase_power" / "output" / "validator_world_rows.tsv", sep="\t")
    ref = ref[ref.arm == "P0"][["seed_index", "reducer", "run", "p_two_sided"]].rename(columns={"p_two_sided": "p_ref"})
    fpr, dist = [], []
    for reducer in ("median", "mean"):
        for run in ("R1", "R2"):
            g = w[(w.reducer == reducer) & (w.run == run)].sort_values("seed_index")
            n = len(g); k = int((g.p_two_sided <= .05).sum()); lo, hi = cp(k, n) if n else (float("nan"),) * 2
            fpr.append(dict(arm="P0F", reducer=reducer, run=run, n_worlds=n, rejections=k,
                            fpr=k / n if n else float("nan"), cp_low=lo, cp_high=hi))
            below = int((g.p_two_sided < .5).sum()); above = int((g.p_two_sided > .75).sum())
            bl, bh = cp(below, n) if n else (float("nan"),) * 2
            al, ah = cp(above, n) if n else (float("nan"),) * 2
            dist.append(dict(arm="P0F", reducer=reducer, run=run, n_worlds=n,
                             mean_p=float(g.p_two_sided.mean()) if n else float("nan"),
                             p_lt_0p5_count=below, p_lt_0p5_fraction=below / n if n else float("nan"),
                             p_lt_0p5_cp_low=bl, p_lt_0p5_cp_high=bh, p_gt_0p75_count=above,
                             p_gt_0p75_fraction=above / n if n else float("nan"),
                             p_gt_0p75_cp_low=al, p_gt_0p75_cp_high=ah))
    paired = []
    for reducer in ("mean", "median"):
        for run in ("R1", "R2"):
            x = w[(w.reducer == reducer) & (w.run == run)][["seed_index", "p_two_sided"]]
            y = ref[(ref.reducer == reducer) & (ref.run == run)][["seed_index", "p_ref"]]
            z = x.merge(y, on="seed_index", validate="one_to_one")
            n = len(z); pf = z.p_two_sided <= .05; pb = z.p_ref <= .05
            both = int((pf & pb).sum()); flat_only = int((pf & ~pb).sum())
            ref_only = int((~pf & pb).sum()); neither = int((~pf & ~pb).sum())
            d, lo, hi = newcombe_method10(flat_only, ref_only, both, n) if n else (float("nan"),) * 3
            paired.append(dict(reducer=reducer, run=run, n=n, both_reject=both,
                               p0f_only_reject=flat_only, p0stageb_only_reject=ref_only, neither_reject=neither,
                               fpr_p0f=(both + flat_only) / n if n else float("nan"),
                               fpr_p0_stageb=(both + ref_only) / n if n else float("nan"),
                               delta_p0f_minus_stageb=d, ci_low=lo, ci_high=hi))
    fprdf = pd.DataFrame(fpr); distdf = pd.DataFrame(dist); pairdf = pd.DataFrame(paired)
    fprdf.to_csv(output / "main_fpr.tsv", sep="\t", index=False)
    distdf.to_csv(output / "main_p_distribution.tsv", sep="\t", index=False)
    pairdf.to_csv(output / "main_paired_contrast.tsv", sep="\t", index=False)
    complete = sorted(w.seed_index.unique().tolist()) == list(range(1, 61)) and \
        all(len(w[(w.reducer == r) & (w.run == u)]) == 60 for r in ("mean", "median") for u in ("R1", "R2"))
    if mode != "formal":
        decision = {"analysis_id": "Bmovie-MPHASE-GAINFLAT-DIAGNOSTIC-STAGE-C-v1.0.0",
                    "status": "INCOMPLETE_SMOKE", "scientific_rules_evaluated": False,
                    "smoke_observed_difference_role": "implementation evidence only"}
    else:
        mf = fprdf[fprdf.reducer == "mean"]; mp = pairdf[pairdf.reducer == "mean"]
        c1 = bool(all(row.ci_high < 0 for row in mp.itertuples()))
        c2 = bool(all(row.cp_low <= .05 <= row.cp_high for row in mf.itertuples()))
        if c1 and c2:
            mech = "gain_coupling_dominant"
        elif c1 and not c2:
            mech = "both_contribute"
        elif (not c1) and (not c2):
            mech = "envelope_or_other_dominant"
        else:
            mech = "indeterminate"
        decision = {"analysis_id": "Bmovie-MPHASE-GAINFLAT-DIAGNOSTIC-STAGE-C-v1.0.0",
                    "status": "COMPLETE" if complete else "INCOMPLETE", "complete": complete,
                    "rules": {"C1_gain_coupling_contribution": c1, "C2_residual_elevation_absent": c2},
                    "mechanism_class": mech if complete else None,
                    "scope": "mechanism diagnostic only; enters no Stage A/B verdict"}
    (output / "main_decision.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")
    return decision


def main():
    p = argparse.ArgumentParser()
    p.add_argument("output", type=Path)
    p.add_argument("--mode", choices=("smoke", "formal"), required=True)
    a = p.parse_args()
    print(json.dumps(summarize(a.output, a.mode), indent=2))


if __name__ == "__main__":
    main()
