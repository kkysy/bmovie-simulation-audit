# Manifest

SHA-256 of every frozen contract and shipped summary table, so the frozen-input chain stays checkable after publication. Table paths are relative to the repository root and mirror the production tree of the analysis workspace.

## Contracts

| file | sha256 |
|---|---|
| `scripts/ieeg/acc00_sim_contract.json` | `AE9B0EEE540CC07AA3D3FE45D445B0C71815E95D7D461C6AB615903532705D09` |
| `scripts/ieeg/acc01_contract.json` | `D15F89687CCB1FEB4D0DFD9C902B12255876CA54551F8014229FD65E9FE32ABD` |
| `scripts/ieeg/paper/acc00_sim_powercurve_contract.json` | `59B22764181C475D824B3847F1FC4F2E3EA726C62618AB02D8691B028031CC55` |
| `scripts/ieeg/paper/acc00_sim_p0_3_ablation_contract.json` | `C7BB5C0024D60D31C398FCC5513E203341B4E5A4B318E1656A91AB33C8DC4354` |
| `scripts/ieeg/paper/acc00_sim_p0_4_panel_contract.json` | `6F8A9518D469322019D73B9FD43165AF60A57B4F0B3DAE9FA928867B7EE0C567` |
| `scripts/ieeg/paper/acc00_sim_p0_6_perturbation_contract.json` | `BB7CB3D53CE724C4349A9ADF874B7AA43BA41D5DC8B1E0B95E4205619AAA916D` |
| `scripts/ieeg/paper/paper_p0_4_blind_seal.json` | `9F96BDED8E492BD61A15AD9A2B84F9B4A9F4595B740143C49FAA8276656235F0` |

Pinned parent-contract digests recorded inside the paper contracts (`parent_*` fields) refer to the `acc00_sim_contract.json` / `acc00_sim_powercurve_contract.json` hashes above.

## Summary tables

| table | sha256 |
|---|---|
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/edgeguard_fix_20260829/validator/task-bangyouredead_desc-acc00sim-feasible-domain-census.tsv` | `FBB8C5976B39C5AACD8B574D1C9A9419A722C9F62A8424F7FBF1C5C536C0E6EB` |
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/edgeguard_fix_20260829/validator/task-bangyouredead_desc-acc00sim-feasible-domain-sessions.tsv` | `DD1C724D52AACD71F4CA9EA9A349A0D0613CA61F5475805A818BF2E901B4BEE1` |
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/edgeguard_fix_20260829/validator/task-bangyouredead_desc-acc00sim-validation.tsv` | `8604732601EB1ABD62A3BEDD27AF5986EB827DBF69A6BB0A6856D284FD4A575B` |
| `processed/subject/group/ieeg_methods_paper/p0_1_calibration_summary/p0_1_additive_recovery.tsv` | `36C4DDF957A3BC95E57AFBF9576F7AA2A9BCD800E1A794746AD8C61784CB3342` |
| `processed/subject/group/ieeg_methods_paper/p0_1_calibration_summary/p0_1_null_calibration.tsv` | `884982E3A58B8B73CD75140020F3CD309287720B8E2856DAF48FBF4F34D2D7FF` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/independent_coverage_check.tsv` | `03A5F3832F73E851AD3EA8037A2073319F18A951D38BBC83D9C2CF4FFEC834F4` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/independent_power_summary.tsv` | `EB681A3CBBF2A9F7A1053A9A6A692E1E1B0DE8A73D836BC89EB10C9563CA1D29` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/p0_2_null_reference.tsv` | `A693E6E48941A58CB1DF9AD427C62068F15F70406C3FE130F96D312046ECF5BF` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/p0_2_powercurve_summary.tsv` | `D4ABDEACA7F6A0FA47C21BA4FFCDB49E034E6D8AEA779A376E0460E45E403263` |
| `processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/tables/p0_2_powercurve_worlds.tsv` | `DF958FECB53B09A272EEC4B1BA3E7C395646C0CDB4A9A83B40513DF88F776E6F` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_independent_verification.tsv` | `09656B30E9C30289D3E5569D30B69E0A8C45150DA3708B5FEBB31A56C9E261E3` |
| `processed/subject/group/ieeg_methods_paper/p0_3_ablation/tables/p0_3_manifest.tsv` | `38850FCD59DD71C29B894C150FF1D57ECDCB98F83B20FC4F184DEDABA348F97A` |
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
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_manifest.tsv` | `C4E089D92A99F42F425A52DA06F4A684EFF8DFC076E3E24B4EC42D72A26B7063` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_member_fpr_summary.tsv` | `CC90E7C2D0FDC28F05AC86F3DD67D4CC70CFCF91812B1D568CDA1478138BB4EE` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_member_worlds.tsv` | `F73EFFF2E51D797C20F561FB7AC13EED58EC378F19B93BA82E42221D4782CF09` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_paired_contrast_MD2.tsv` | `5A800C98EC6173DD2252D83B315E5C1B5A3321665ACF082B6717FD0CF5B4AA8F` |
| `processed/subject/group/ieeg_methods_paper/p0_4_diagnostic_panel/tables/p0_4_validator.tsv` | `4B426C744315335EDD8747D4B01B3375D601905B2272E6D21C34CE6CAFC60F99` |
| `processed/subject/group/ieeg_p0_6_perturbation/census_sensitivity.tsv` | `CBDE49155D7F47ADFEA2010866D62BE260D8647BBE6501AA42D20E9D3C17C800` |
| `processed/subject/group/ieeg_p0_6_perturbation/corner_probe_G06.tsv` | `C8E7A6047890F67B3E9B6F420E4668F11CABCAA38E25D06E27431F7C1B7C5969` |
| `processed/subject/group/ieeg_p0_6_perturbation/paired_fpr_contrast.tsv` | `B5ADA9414D6498384F649A7F978C7CB5087DF09C3218DB6EDC48EB767E4BE457` |
| `processed/subject/group/ieeg_p0_6_perturbation/per_world_metrics.tsv` | `ED564D779A2141980738C0100B846431261F26D28201CA1DFA7CBE3FAFEE6088` |
| `processed/subject/group/ieeg_p0_6_perturbation/primary_verdict.json` | `62AFC3FC9389D5B9EEF51C7728C24A89F52414AFD5FD6192D28622EEF33823D9` |
| `processed/subject/group/ieeg_p0_6_perturbation/run_manifest.json` | `77F04E547BD2ED6451CCEBB63B4EE05B4D8AB2DDDF78EB29B283121BD0A699C1` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_0p5x.tsv` | `88DF4AF4D52C807689A381EAC52C83A703F97F74FAC33F5698334A06A58253B5` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_1p5x.tsv` | `0187ABC3D09148F3C2AD46531ACDD9B2183A654BBAEFD34D2A1974BFC12598A8` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_census_eps_0.025xSD.tsv` | `BF1D84462EA32BDB83B8C2C6007F711B4C1DF023C105A076D759DCA48004B27D` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_census_eps_0.10xSD.tsv` | `35B8B473B341C831676B3EDFA30DEA5C70F1F4086B9F59040E6E68EDAE3DAB69` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_census_intervalperm_1.tsv` | `02F8FF52ADE6B4881735871208FBD6F4C1AC7F977D466F5FC2B279968C1E0130` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_census_intervalperm_2.tsv` | `E9559F3716FB9386138B57959C4581810BABB9770B23B74DCA31151826D2EE9C` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_manifest.tsv` | `7958FF82333A232C0315702B885CA98EE0468DF978AC903A5E0C0FE47F362794` |
| `processed/subject/group/ieeg_p0_6_perturbation/validator_report.json` | `4032623136B64AA40B08F7DA8D565010633FF1E727B99CBEF09704253C866B60` |
