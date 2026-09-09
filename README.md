# Audit-first simulation evaluation of event-related iEEG estimators (Bmovie schedule)

Code, frozen simulation contracts, and summary tables for the methods paper on
null calibration, empirical power, mechanism ablation, a blinded diagnostic
panel, and surrogate feasible-domain auditing for event-related intracranial
EEG analysis under a dense, irregular naturalistic event schedule.

The audit is fully synthetic: **no real neural signals are used or analyzed**.
Worlds preserve the Bmovie event schedule, clean support, and hierarchical
nesting (participant → session → pair) but contain only simulated 1/f
background plus injected transients.

## Programme

| phase | component | entry points |
|---|---|---|
| P0-1 | null calibration (400 worlds) + Mpower recovery pool (300 worlds) | `run_acc00_nullpool.m`, `run_acc00_recovery_safe.m`, `summarize_p0_1_calibration.py` |
| P0-2 | Mpower power curve, 986 worlds over seven additive cells | `run_acc00_powercurve_main.m`, `summarize_p0_2_powercurve.py`, `validate_acc00_powercurve.py` |
| P0-3 | mechanism ablations, 768 worlds (Mpower arm + Mphase guard replays) | `run_acc00_p0_3_ablation_main.m`, `summarize_p0_3_ablation.py`, `validate_acc00_p0_3_ablation.m` |
| P0-4 | preregistered blinded four-member diagnostic panel, 480 worlds | `run_acc00_p0_4_panel_main.m`, `summarize_p0_4_panel.py`, `validate_acc00_p0_4_panel.m` |
| P0-5 | surrogate-selection record and feasible-domain census | `acc01_build_empirical_surrogates.m`, census table under `processed/.../edgeguard_fix_20260829/validator/` |
| P0-6 | schedule-density perturbation with paired FPR contrasts and census sensitivity checks | `run_acc00_p0_6_main.m`, `summarize_acc00_p0_6.m`, `validate_acc00_p0_6_perturbation.m`, `verify_p0_6_paired_contrast.py` |

## Repository layout

Paths mirror the production analysis workspace so that the contracts' internal
path references and the paper's source map stay valid.

```
scripts/ieeg/                      frozen contracts + runners + validators + summarizers
  acc00_sim_contract.json          parent contract (schedule, guards, surrogate domain, seeds)
  acc01/                           synthetic-world generator and shared estimator kernels
  paper/                           P0-2..P0-6 contracts, runners, validators, summarizers
processed/subject/group/           summary tables underlying all figures and numerical claims
manuscript/figure/                 figure generation script; rendered manuscript figures are not shipped
MANIFEST.md                        SHA-256 of every contract and shipped table
```

## Environment

- MATLAB with the Statistics and Machine Learning Toolbox (`robustfit`) and the
  Parallel Computing Toolbox (`parfor`). Runners are launched from the
  repository root; each entry point resolves the root itself.
- Python ≥ 3.10: `pip install -r requirements.txt` (numpy, scipy, pandas,
  h5py, matplotlib).
- The P0-6 contract writer accepts `BMOVIE_SOURCE_ROOT` for local external
  inputs and writes only repository-relative paths into the public contract.

## Reproduction chain

1. Contracts are frozen first and hashed (`*.sha256`, `MANIFEST.md`); every
   runner asserts contract and input hashes before computing.
2. World generation and scoring are deterministic (Threefry streams; world
   seeds recorded in every checkpoint). P0-4 pins `maxNumCompThreads(1)` in
   client and workers so recomputed arithmetic reproduces production bitwise.
3. Summarizers are read-only over completed checkpoints and write the tables
   under `processed/subject/group/`.
4. Independent validators re-derive FPRs, Clopper–Pearson intervals, and
   pairing identities from the world-level tables.
5. Figures: `python manuscript/figure/generate_figures.py` regenerates all six
   figures from the shipped tables (runs its own layout checks). The rendered
   manuscript figures are intentionally kept outside this public release.

## Data source and citation

The event schedule, clean support, and recording hierarchy analyzed here come
from the Bmovie naturalistic movie-viewing iEEG dataset:

> Keles, U., Dubois, J., Le, K. J. M., Tyszka, J. M., Kahn, D. A., Reed, C. M.,
> Chung, J. M., Mamelak, A. N., Adolphs, R., & Rutishauser, U. (2024).
> *Multimodal single-neuron, intracranial EEG, and fMRI brain responses during
> movie watching in human patients.* Scientific Data, 11.
> https://doi.org/10.1038/s41597-024-03029-1

Dataset access: DANDI 000623 (https://dandiarchive.org/dandiset/000623);
NWB BIDS release: https://github.com/rutishauserlab/bmovie-release-NWB-BIDS.

Contract paths are repository-relative. The raw NWB files and other
external-input documents referenced by the contracts are not redistributed;
their expected relative locations are retained only for provenance and
independent re-use when those inputs are available locally.

## Data boundary

This repository contains no video, no physiological recordings, and no
per-event schedule: the movie stimulus is not redistributed (it was obtained
from the original authors by request and is not part of this audit), and the
frozen contracts pin their real-data-derived inputs (the causal
prediction-difficulty event label, event covariates, and bad-segment support)
by path and SHA-256 only. Those inputs are not redistributed here; summary
tables and contracts are sufficient to reproduce every figure and numerical
claim in the paper.

## Status

The blind seal for P0-4 (`paper_p0_4_blind_seal.json`) is published as
unblinded; scoring used opaque member IDs only, per the preregistered protocol.
License: to be added on release.
