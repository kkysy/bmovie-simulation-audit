# Stage A surrogate-count sensitivity (J = 32 vs J = 256)

Matched sensitivity check on the Stage A confirmation worlds: the same 120
pure-null worlds (identical world seeds, pair-level common-random-number
states) were rescored with the refitted-surrogate target raised from J = 32
to J = 256, with exhaustive enumeration for pairs whose feasible candidate
count did not exceed 256.

Verdict: `ROBUST_TO_SURROGATE_COUNT` — all three preregistered rule families
(mean-statistic stability within the ±0.10 Stage-A SD equivalence band with
RMSE ratio ≤ 0.25; J = 256 calibration with the mean interval containing
0.05 and the median lower bound above 0.05; paired mean-minus-median
Newcombe upper bound below zero) passed in both runs, and the J = 32 arm
reproduced the Stage A reference tables exactly.

These are release ports of the analysis-workspace tools in
`tools/confirm_mphase_jrobust/`; the only changes are the path bootstrap
lines (the workspace originals resolved the workspace root from the
`tools/<stage>/` layout; the ports expect to be run from the repository
root, and the Python summarizer/validator locate the repository root by
walking up to `MANIFEST.md`). The contract ships byte-identical to the
frozen workspace version.

Two workspace-layout references inside the frozen contract do not resolve
literally in this repository: the Stage A reference contract and runner are
pinned by their workspace hashes, while this repository ships the Stage A
port under `stage_a_confirmation/` (byte-patched bootstrap, hence a
different contract hash), and the Stage A reference tables pinned in the
contract as `tools/confirm_mphase_reducer/output/validator_*_rows.tsv` are
the same bytes shipped under `stage_a_confirmation/` (the gzipped pair table
decompresses to the pinned SHA-256 `DC69015D…`). World-level MATLAB
computation additionally requires the analysis workspace's full input tree;
the shipped artifacts support table- and hash-level auditing, which is what
the manuscript's numerical claims rely on.

Result tables live under
`processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/`,
including `validator_pair_rows.tsv.gz` (the atomic pair-level records for
both arms) and the independent validator report (12 checks, PASS).
