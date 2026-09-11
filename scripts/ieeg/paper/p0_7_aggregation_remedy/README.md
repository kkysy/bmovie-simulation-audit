# P0-7 — Aggregation remedy and injection-matched validation (Mphase route)

Three sub-stages, each in its own folder with its frozen contract, MATLAB
runners, Python summarizer, and independent validator. These are release ports
of the analysis-workspace tools of the same names; the only changes are the
path bootstrap lines (the workspace originals resolved the workspace root from
the `tools/<stage>/` layout; the ports expect to be run from the repository
root, and the Python summarizers/validators locate `scripts/ieeg/paper/` by
walking up to `MANIFEST.md`). Contracts and world generators ship
byte-identical to the frozen workspace versions.

| sub-stage | question | verdict |
|---|---|---|
| `stage_a_confirmation/` | Does replacing the participant-median reducer with the participant mean restore pure-null Mphase calibration? (120 fresh pure-null worlds, both reducers on identical pair-level values) | `CONFIRMED` |
| `stage_b_injection_matched/` | Does the repaired route stay calibrated under an energy-matched random-phase injection (P0), and does it detect phase concentration (P4, κ = 4)? (60 matched worlds × 2 arms) | `NOT_CONFIRMED` — P0 over-rejected (0.200/0.117); all rejections in both arms were negative-direction, so the P4−P0 contrast does not index phase-detection power |
| `stage_c_gain_flat/` | Is the residual inflation driven by the coupling between injection gain and the `z_e`-maxima event selection? (flat gain `b = 0`, same 60 backgrounds) | `indeterminate` — FPR returned to 0.067/0.067 (CI contains 0.05); paired reduction vs P0 significant in R1 (−0.133 [−0.253, −0.017]), inconclusive in R2 |
| `stage_a_jrobust/` | Is the Stage A reducer conclusion driven by Monte Carlo variability from the production 32 refitted surrogates? (the same 120 null worlds rescored at J = 256 with pair-level common-random-number states) | `ROBUST_TO_SURROGATE_COUNT` — all three preregistered rule families passed in both runs; the J = 32 arm reproduced Stage A exactly |

Result tables live under
`processed/subject/group/ieeg_p0_7_aggregation_remedy/<sub-stage>/`, including
`validator_pair_rows.tsv.gz` (the atomic pair-level records underlying the
direction audit) and the world/subject-level validator bridges.
