# Manifest

SHA-256 of every frozen contract and shipped summary table, so the frozen-input chain stays checkable after publication. Table paths are relative to the repository root and mirror the production tree of the analysis workspace.

## Contracts

| file | sha256 |
|---|---|
| `scripts/ieeg/acc00_sim_contract.json` | `51F06B9A733D6FF88534FAECE148E5AB89EDB17B6372657313BF3D4D3BF53260` |
| `scripts/ieeg/acc01_contract.json` | `F2C874719CA9F83010666538017570D4980E0E7D622A0BC0A1E5F4501368D2C9` |
| `scripts/ieeg/paper/acc00_sim_powercurve_contract.json` | `E3C35D065F6819AA611B7A60573689BFDCED06E399F80B6DAAE40F46BE3E1863` |
| `scripts/ieeg/paper/acc00_sim_p0_3_ablation_contract.json` | `C27DD782384AEAED32995E9C2C508A38D57349BFEA3D3EC7BC2E4D08D983F949` |
| `scripts/ieeg/paper/acc00_sim_p0_4_panel_contract.json` | `F6D6A81C7347D7261899F66FC7147E61FFA9E36D24B0EAC9C3B3FFD578F9D661` |
| `scripts/ieeg/paper/paper_p0_4_blind_seal.json` | `9F96BDED8E492BD61A15AD9A2B84F9B4A9F4595B740143C49FAA8276656235F0` |

Pinned parent-contract digests recorded inside the paper contracts (`parent_*` fields) refer to the `acc00_sim_contract.json` / `acc00_sim_powercurve_contract.json` hashes above.

## Summary tables

| table | sha256 |
|---|---|
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/edgeguard_fix_20260829/validator/task-bangyouredead_desc-acc00sim-feasible-domain-census.tsv` | `FBB8C5976B39C5AACD8B574D1C9A9419A722C9F62A8424F7FBF1C5C536C0E6EB` |
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/edgeguard_fix_20260829/validator/task-bangyouredead_desc-acc00sim-feasible-domain-sessions.tsv` | `DD1C724D52AACD71F4CA9EA9A349A0D0613CA61F5475805A818BF2E901B4BEE1` |
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/edgeguard_fix_20260829/validator/task-bangyouredead_desc-acc00sim-validation.tsv` | `F1F5B8E5CBB8A191BB8BFDC3790AED7E66102ED66BB5897B31F5C67E1A94D1F8` |
| `processed/subject/group/ieeg_methods_paper/p0_1_calibration_summary/p0_1_additive_recovery.tsv` | `36C4DDF957A3BC95E57AFBF9576F7AA2A9BCD800E1A794746AD8C61784CB3342` |
| `processed/subject/group/ieeg_methods_paper/p0_1_calibration_summary/p0_1_null_calibration.tsv` | `884982E3A58B8B73CD75140020F3CD309287720B8E2856DAF48FBF4F34D2D7FF` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/independent_coverage_check.tsv` | `03A5F3832F73E851AD3EA8037A2073319F18A951D38BBC83D9C2CF4FFEC834F4` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/independent_power_summary.tsv` | `EB681A3CBBF2A9F7A1053A9A6A692E1E1B0DE8A73D836BC89EB10C9563CA1D29` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/p0_2_null_reference.tsv` | `A693E6E48941A58CB1DF9AD427C62068F15F70406C3FE130F96D312046ECF5BF` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/p0_2_powercurve_summary.tsv` | `8A357D58915CDD7433530B9EB9AEDBC648B000DCFB8EA2DB586A2808043D2A5D` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/p0_2_powercurve_worlds.tsv` | `4C89EA30E4C45A1E86C04F5A9B2A74BB16F8D9970F349C81EADA1F82D69206EA` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_independent_verification.tsv` | `09656B30E9C30289D3E5569D30B69E0A8C45150DA3708B5FEBB31A56C9E261E3` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_manifest.tsv` | `88E7B374F4A566CBD58C2AD20816F6F3FB3B8CE9790011047C9BF200E9F50EB6` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_mphase_historical_guard_summary.tsv` | `35B7BE41653831909A0344E3E825C4BEF70503E579CBD51EA1274CA20513A4F1` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_mphase_replay_shape_summary.tsv` | `BE7D8800A0A676AEAA032D21805067B29F171F51F6E7094A7742F25EDA5975B3` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_mphase_worlds.tsv` | `728C1EE4FF2EA13EDCF7DE52A2DF8FDB2897B7712813BB1D93C6D863D4498933` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_mpower_cell_summary.tsv` | `45D828112DE2F01844C33BAEADCBF41505BCEB8E5A83A6EA2456A6A456CADFD8` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_mpower_factor_contrasts.tsv` | `9294B6D1FED283D243F8AF34D2EB9D5A1939F4EE8CD314A0D977C777ED7C1E70` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_mpower_pair_diagnostics.tsv.gz` | `9562086EAA1E768E8BD89605184D21A5B1133DCC94775F910DBE6D177298AD62` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_mpower_worlds.tsv` | `5FF4C5EBAA26A4806F585E985165F056ECF1093F0F2CB4D8C6F40D07C0DC3515` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_subject_vectors_bridge.tsv` | `0CA72162E1FF33D239038551B09CD293EF19C1AFBDFD06E2B8385AC979E217CD` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_validator.tsv` | `C40E10880974CCEA2F2A4F7FEB05A8AFCF046329CD10962970B46EF44F2A00E2` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p10_fft_leakage_description.tsv` | `28674F2D93B9D682180FD9802A2BFE999794082181C2CE3BA5B7C2B27CDB0E84` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_blind_decode.tsv` | `91F2C6A27863FA4A299B57DF8C417F81D2B2F8490B4FDC1F906FCDEA5A7DE150` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_exploratory_MD1.tsv` | `F809E19A9184934CCF44E16E130913C309788E1B4089E4294305F40E2E6DF9F7` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_manifest.tsv` | `8DA8D20573E3184767422F994C5FF6DA9C16529D211E4734937D4E97BD57AD8D` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_member_fpr_summary.tsv` | `CC90E7C2D0FDC28F05AC86F3DD67D4CC70CFCF91812B1D568CDA1478138BB4EE` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_member_worlds.tsv` | `F73EFFF2E51D797C20F561FB7AC13EED58EC378F19B93BA82E42221D4782CF09` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_paired_contrast_MD2.tsv` | `5A800C98EC6173DD2252D83B315E5C1B5A3321665ACF082B6717FD0CF5B4AA8F` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_validator.tsv` | `4B426C744315335EDD8747D4B01B3375D601905B2272E6D21C34CE6CAFC60F99` |
