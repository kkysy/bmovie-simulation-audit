#!/usr/bin/env python3
"""Mechanically derive the frozen P0-4 null-panel contract from P0-3."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PARENT = HERE / "acc00_sim_p0_3_ablation_contract.json"
OUT = HERE / "acc00_sim_p0_4_panel_contract.json"
PARENT_SHA256 = "C27DD782384AEAED32995E9C2C508A38D57349BFEA3D3EC7BC2E4D08D983F949"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def build(parent: dict) -> dict:
    """Apply only the v3 §6 parent-merge overrides/removals."""
    contract = copy.deepcopy(parent)
    for key in ("cells", "event_mask", "spectral_shape", "seed"):
        contract.pop(key, None)
    contract.update({
        "analysis_id": "ACC00-SIM-P0-4-DIAGNOSTIC-PANEL-v1.0.0",
        "status": "paper_only_frozen_implementation_not_authorized",
        "purpose": "pre-registered blinded null-FPR diagnostic panel; implementation and smoke only, no main run",
    })
    contract["additive"]["A_uV"] = 0.0
    contract["additive"]["A_source"] = "P0-4 pure-null panel"
    contract["surrogate"].update({
        "domain_guard_edge_lower_s": 1.5,
        "domain_guard_edge_upper_s": 1.5,
        "guard_is_cell_specific_only_for": ["M-D2"],
    })
    contract["panel"] = {
        "null_worlds": 120,
        "base_seed": 91000000,
        "scenario_index": 5,
        "effect_index": 0,
        "world_seed_range": [91500001, 91500120],
        "crn_paired_across_members": True,
        "member_id_in_seed": False,
        "holdout_relative_to_historical_pools": True,
        "opaque_members": ["PanelMember-01", "PanelMember-02", "PanelMember-03", "PanelMember-04"],
        "pre_registered_directional_contrast_opaque_ids": ["PanelMember-03", "PanelMember-02"],
        "pre_registered_exploratory_opaque_id": "PanelMember-04",
        "member_semantics_seal": "paper_p0_4_blind_seal.json",
    }
    contract["outputs"] = {
        "root": "processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel",
        "checkpoints": "processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/checkpoints/<opaque_id>/world_<seed_index>.mat",
        "tables": "processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables",
        "checkpoint_format": "MATLAB -v7.3; write .tmp.mat then atomic move",
        "overwrite_policy": "identity-verified resume only; mismatch stops",
    }
    contract["execution"] = {
        "implementation_authorized": False,
        "smoke_authorized": False,
        "main_run_authorized": True,
        "parallel_workers": 12,
        "worker_rationale": "single-world peak about 3.1GB; run members sequentially; do not co-schedule member classes",
        "total_new_world_evaluations": 480,
        "no_new_gate": True,
    }
    contract["parent_p0_3"] = {"path": "scripts/ieeg/paper/acc00_sim_p0_3_ablation_contract.json", "sha256": PARENT_SHA256, "immutable": True}
    return contract


def main() -> int:
    actual = sha256(PARENT)
    if actual != PARENT_SHA256:
        raise RuntimeError(f"immutable P0-3 parent SHA mismatch: {actual}")
    OUT.write_text(json.dumps(build(json.loads(PARENT.read_text(encoding="utf-8"))), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {OUT.name} SHA256={sha256(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
