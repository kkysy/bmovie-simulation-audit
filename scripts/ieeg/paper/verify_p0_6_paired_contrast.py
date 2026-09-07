"""Independent Python verification of P0-6 paired FPR contrast and verdict.

Implements Newcombe (1998b, Stat Med 17:2635-2650) method 10 for the paired
difference of binomial proportions from the paper's closed form (phi-correlation
correction with the method-10 continuity correction), recomputes every row of
paired_fpr_contrast.tsv from per_world_metrics.tsv, and re-derives the verdict
labels from the frozen thresholds.

Usage: python verify_p0_6_paired_contrast.py <root> [partial]
Exit 0 = all checks pass.
"""
import json
import math
import sys
from pathlib import Path

Z = 1.959963984540054
ALPHA = 0.025


def wilson(k, n):
    p = k / n
    den = 1 + Z * Z / n
    cen = p + Z * Z / (2 * n)
    rad = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n))
    return (cen - rad) / den, (cen + rad) / den


def newcombe_method10(b, c, e, n):
    """Newcombe 1998b method 10 paired-difference CI.

    Cells of the paired 2x2 table: e=both reject, b=new-only reject,
    c=old-only reject, h=neither; n=e+b+c+h.
    """
    h = n - b - c - e
    p1 = (e + b) / n
    p2 = (e + c) / n
    lb, ub = wilson(e + b, n)
    lc, uc = wilson(e + c, n)
    dl2, du2 = p1 - lb, ub - p1
    dl3, du3 = p2 - lc, uc - p2
    den = math.sqrt((e + b) * (c + h) * (e + c) * (b + h))
    num = e * h - b * c
    if den > 0:
        if num > 0:
            num = max(num - n / 2, 0)  # method-10 continuity correction
        phi = num / den
    else:
        phi = 0.0
    d = (b - c) / n
    lo = d - math.sqrt(dl2 ** 2 - 2 * phi * dl2 * du3 + du3 ** 2)
    hi = d + math.sqrt(du2 ** 2 - 2 * phi * du2 * dl3 + dl3 ** 2)
    return d, lo, hi


def clopper_pearson(n, k):
    from scipy.stats import beta
    lo = 0.0 if k == 0 else beta.ppf(ALPHA, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - ALPHA, k + 1, n - k)
    return lo, hi


def load_tsv(path):
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            z = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, z)))
    return rows


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    od = root / "processed" / "subject" / "group" / "ieeg_p0_6_perturbation"

    # per-world metrics: key = (arm, run, seed_index)
    pw = {}
    for r in load_tsv(od / "per_world_metrics.tsv"):
        pw[(r["arm"], r["run"], int(r["seed_index"]))] = {
            "Mpower_reject": r["Mpower_reject"] == "1",
            "Mphase_reject": r["Mphase_reject"] == "1",
        }

    # arm candidates per (family, density)
    cand = {
        ("Mpower", "1p5x"): "Mpower_1p5",
        ("Mphase", "1p5x"): "Mphase_1p5",
        ("Mphase", "0p5x"): "Mphase_0p5",
    }
    errors = []

    for row in load_tsv(od / "paired_fpr_contrast.tsv"):
        fam, run, density = row["family"], row["run"], row["density"]
        col = "Mpower_reject" if fam == "Mpower" else "Mphase_reject"
        arm = cand.get((fam, density))
        if arm is None:
            errors.append(f"{fam}/{density}/{run}: unknown contrast row")
            continue
        seeds = sorted(k[2] for k in pw if k[0] == arm and k[1] == run)
        xs = [pw[(arm, run, s)][col] for s in seeds]
        ys = [pw[("historical_1p0", run, s)][col] for s in seeds]
        b = sum(x and not y for x, y in zip(xs, ys))
        c = sum(not x and y for x, y in zip(xs, ys))
        e = sum(x and y for x, y in zip(xs, ys))
        n_ = len(xs)
        if n_ != int(row["n"]):
            errors.append(f"{fam}/{run} ({arm}): n row={row['n']} recomputed={n_}")
        if (b, c) != (int(row["b"]), int(row["c"])):
            errors.append(f"{fam}/{run} ({arm}): discordant counts row=({row['b']},{row['c']}) recomputed=({b},{c})")
        d, lo, hi = newcombe_method10(b, c, e, n_)
        for name, got, want in (("Delta", float(row["Delta"]), d),
                                ("CI_lower", float(row["CI_lower"]), lo),
                                ("CI_upper", float(row["CI_upper"]), hi)):
            if abs(got - want) > 1e-12:
                errors.append(f"{fam}/{run} ({arm}): {name} row={got:.15f} recomputed={want:.15f}")
        if abs(sum(xs) / n_ - float(row["FPR_1p5"])) > 1e-12 or abs(sum(ys) / n_ - float(row["FPR_1p0"])) > 1e-12:
            errors.append(f"{fam}/{density}/{run} ({arm}): FPR point estimates mismatch")
        cp_lo, cp_hi = clopper_pearson(n_, sum(xs))
        cp0_lo, cp0_hi = clopper_pearson(n_, sum(ys))
        if abs(cp_lo - float(row["FPR_1p5_CP_low"])) > 1e-9 or abs(cp_hi - float(row["FPR_1p5_CP_high"])) > 1e-9 \
                or abs(cp0_lo - float(row["FPR_1p0_CP_low"])) > 1e-9 or abs(cp0_hi - float(row["FPR_1p0_CP_high"])) > 1e-9:
            errors.append(f"{fam}/{run} ({arm}): CP intervals mismatch")

    verdict = json.loads((od / "primary_verdict.json").read_text(encoding="utf-8"))
    for item in verdict["items"]:
        lo, hi = item["CI_lower"], item["CI_upper"]
        want = "attenuation" if hi <= -0.05 else ("persistence" if lo >= -0.02 else "indeterminate")
        if item["label"] != want:
            errors.append(f"verdict {item['run']}: label {item['label']} != {want}")
        if item["text"] not in (
                "Attenuation: the paired 95% Newcombe confidence interval for FPR at the high-density schedule minus FPR at the original schedule has an upper bound no greater than -0.05.",
                "Persistence: the paired 95% Newcombe confidence interval for FPR at the high-density schedule minus FPR at the original schedule has a lower bound no less than -0.02.",
                "The density perturbation did not yield a decisive attenuation or persistence verdict under the preregistered precision rule."):
            errors.append(f"verdict {item['run']}: non-preregistered text")

    if errors:
        print("FAIL")
        for e in errors:
            print(" -", e)
        sys.exit(1)
    print("PASS: paired contrast rows, Newcombe method 10, CP intervals, and verdict all verified")


if __name__ == "__main__":
    main()
