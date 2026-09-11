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
| `scripts/ieeg/paper/p0_7_aggregation_remedy/stage_a_confirmation/confirm_mphase_reducer_contract.json` | `2D1181CD1AFBE43BE01435598F22F793760613BA881C55DB1B1A725B8253401E` |
| `scripts/ieeg/paper/p0_7_aggregation_remedy/stage_b_injection_matched/stageb_mphase_power_contract.json` | `0B10CF92B727D4EE7196AB3D071B2ACD94A032BE63BC3A7F1D543E91FC19A577` |
| `scripts/ieeg/paper/p0_7_aggregation_remedy/stage_c_gain_flat/stagec_mphase_gainflat_contract.json` | `1E66EB5835C263260559A72391F60D6A7CF69CEA911447B015BAE6C49487752D` |
| `scripts/ieeg/paper/p0_7_aggregation_remedy/stage_a_jrobust/confirm_mphase_jrobust_contract.json` | `B0421FF08D5B8807E66B8E73E146B6DCE1BE85D1D041A6D1257BA573B8AC9E1F` |

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
| `processed/subject/group/ieeg_p0_6_perturbation/corner_probe_G06.tsv` | `A3AE21AF026C2A22ACA4AE7E9C4F868DF2F291C0C573995EFAC83ADD16CACB59` |
| `processed/subject/group/ieeg_p0_6_perturbation/paired_fpr_contrast.tsv` | `B5ADA9414D6498384F649A7F978C7CB5087DF09C3218DB6EDC48EB767E4BE457` |
| `processed/subject/group/ieeg_p0_6_perturbation/per_world_metrics.tsv` | `190673115E1C167A8A9988F2B68D1FDFCAA21D0D6688FF7144E989567E2E2271` |
| `processed/subject/group/ieeg_p0_6_perturbation/primary_verdict.json` | `62AFC3FC9389D5B9EEF51C7728C24A89F52414AFD5FD6192D28622EEF33823D9` |
| `processed/subject/group/ieeg_p0_6_perturbation/run_manifest.json` | `77F04E547BD2ED6451CCEBB63B4EE05B4D8AB2DDDF78EB29B283121BD0A699C1` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_0p5x.tsv` | `88DF4AF4D52C807689A381EAC52C83A703F97F74FAC33F5698334A06A58253B5` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_1p5x.tsv` | `0187ABC3D09148F3C2AD46531ACDD9B2183A654BBAEFD34D2A1974BFC12598A8` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_census_eps_0.025xSD.tsv` | `BF1D84462EA32BDB83B8C2C6007F711B4C1DF023C105A076D759DCA48004B27D` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_census_eps_0.10xSD.tsv` | `35B8B473B341C831676B3EDFA30DEA5C70F1F4086B9F59040E6E68EDAE3DAB69` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_census_intervalperm_1.tsv` | `02F8FF52ADE6B4881735871208FBD6F4C1AC7F977D466F5FC2B279968C1E0130` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_census_intervalperm_2.tsv` | `E9559F3716FB9386138B57959C4581810BABB9770B23B74DCA31151826D2EE9C` |
| `processed/subject/group/ieeg_p0_6_perturbation/schedules/schedule_manifest.tsv` | `7958FF82333A232C0315702B885CA98EE0468DF978AC903A5E0C0FE47F362794` |
| `processed/subject/group/ieeg_p0_6_perturbation/validator_report.json` | `5B0F211099530C91B223161CBBBF1E5D41853D826DB5EAFE126723B579B79F68` |
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/task-bangyouredead_desc-acc00sim-covariate-manifest.tsv` | `85AAB391C57DE4FD2B144720ABC80E3C69E662ABE4DAF8B1FBF94BF29082AB3A` |
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/task-bangyouredead_desc-acc00sim-feasible-domain-census.tsv` | `FBB8C5976B39C5AACD8B574D1C9A9419A722C9F62A8424F7FBF1C5C536C0E6EB` |
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/task-bangyouredead_desc-acc00sim-feasible-domain-sessions.tsv` | `DD1C724D52AACD71F4CA9EA9A349A0D0613CA61F5475805A818BF2E901B4BEE1` |
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/task-bangyouredead_desc-acc00sim-label-covariates.tsv` | `FB371D5DB3F23439BB2569DD8CB2E5159475B5D8D3CDE6B4806FF74D8BC44CA2` |
| `processed/subject/group/ieeg_acc00_sim/BangYoureDead/task-bangyouredead_desc-acc00sim-validation.tsv` | `F8B96FAAC0D5F5295677F70D999D80406F5DB1A2D0ABAC4037782A297B270504` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/bridge_report.json` | `15EA1B94B4F939570010140E666B53B2453505193D711566978E1C89B132904F` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/main_decision.json` | `8C1658DB53FF154C7794B637B2A4E56EBDABD63581C73FCE9D6463E7AF2346F0` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/main_fpr.tsv` | `98F56FBF3B956873FD06D4EEF672F4E26FF41A740B02141A601422A19E448D91` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/main_p_distribution.tsv` | `FC2F13426CB6FA96DB24A730259E23A60B1CCB4AD7B6694E1E7AF4D9605DBC6B` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/main_paired_risk_difference.tsv` | `55E44F37872EF9F3828F5910E15E0309021093A28AC17EA11A036D856F4B7A4A` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/run_manifest.json` | `4C5CCC4107DDEF5C7144D95C56F3705AD200DC82B3490F345465C025FAF33611` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/validator_pair_rows.tsv.gz` | `FBD365B474517D76DB77D5E698B075EF5D9715ACD769BFCBEE8496ED4B2C709F` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/validator_report.json` | `4CC9AF5089B61E64328603CF7F5AD4ACBCE4513237E395B25CB0B4A5CECB869E` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/validator_subject_rows.tsv` | `F9F89E1D9D63B794CFFF18FAD8785B6E3A39899180088F7223802AA669A0AAA4` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_confirmation/validator_world_rows.tsv` | `233C73B1F806E77C9A50F4E92C2C021FA5A5B6A8FC0B5D4DBBA134FA8E46B03B` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/bridge_report.json` | `0641F3ECDBC2A4E44F93A322D28007DA9D0417BD9B80A4EEB611943A789C61D3` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/main_decision.json` | `C1255375AFB99FFF45F7B29876203C7073AE752F3129385387B0161EAD40E055` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/main_fpr.tsv` | `9B4FAFAE58A8B502F95F091AEFD301734A12E39E97CC1048BE8578C363B7FDA3` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/main_p_distribution.tsv` | `587DF87F64454710F87C045E035F19D3C58DD6DE192BFC3ED30F94FDC71E7438` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/main_paired_contrast.tsv` | `E922D1967F9C44EB4B7F453F1498E148E6A79666A68C5291F9ED5394FB191F48` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/run_manifest.json` | `3D4347C0BBF0914F452DF6BBDEAD131A56D5FC4F36A233333652B0B8537BFCF7` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/validator_pair_rows.tsv.gz` | `9089F39DFE81018D825D6A34BCDE54FB6FD71A89AEC1FF0F04B7F15F8B930E3F` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/validator_report.json` | `60FFE98222888DD459793403E279FE1CE271738933EAC6120666A9C735438E70` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/validator_subject_rows.tsv` | `193C39D29C46BAD0D612E57E9AA40528E1149F82694D241D772059E007D286A4` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/validator_world_rows.tsv` | `99D58E3428D2A41FCA5F9FF2DFBE32A21EF8301CD30E4E6B8633CE0B8EFA2DD1` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_b_injection_matched/world_ppc_means.tsv` | `D30AC550DEE7C08155D2973896AC2C6D79AEA6414E6702F1E0104CF4614782B5` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/bridge_report.json` | `901E87C8BF5226F7664D8D069D2EE63643315C9AF9E7EBA9AF9DBEC529345336` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/main_decision.json` | `94EF3EA5B9062E30B542CE48389AEBBBA44CFFD2B6E9508665DC55A282204121` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/main_fpr.tsv` | `8EDF62E36513E9C487C023A78F7A8D4A39F9816AFC28FD53187033ABD2DD1B49` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/main_p_distribution.tsv` | `EA9D8876A6B93085F7D5D02033061B668EA97CE4A65C374E0DDBE3942322CD81` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/main_paired_contrast.tsv` | `7F8687F3A7CBB0A8CF8B2046B47D06BAA59ECE3F31DA657C2589080E9A44557A` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/run_manifest.json` | `235C94A8EDF36E6BB58B092B5AB4CFF4DEECED2BE0B5DAF6BC0069BDEC76290F` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/validator_pair_rows.tsv.gz` | `842E093C8CA385802D3F6AC0C5B0EA78B56953129291FF2F846A9B45E05C5A6B` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/validator_report.json` | `BE30A78CC267C4F0FB87841ED6DB75A17CA8631064994160F0B4BEB03EB0850F` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/validator_subject_rows.tsv` | `D7376EF0C647FB8B158046D73432F15D8A8D090C221759DBBB8DB8BFFA150B55` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_c_gain_flat/validator_world_rows.tsv` | `9CEFB5EBA7C8A7FDE6F0FE64E49B6AFEFA551E5767CD478664D3CF2FF8092338` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/bridge_report.json` | `0333F82F8CF00E7CD5249725E95C343208AC154B70BC02B38B217FAA4618EABA` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/main_decision.json` | `7149B451052CF6E559A45B89DDCCF249AFB0E62FEA8057FF22FE66172853FCC9` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/main_fpr.tsv` | `A8D6630499305759DA0F23F1221660C4E6A3A4EE4B9CBCFE16AC8090E76A8678` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/main_j_paired_difference.tsv` | `07CCEA6C14287D2D3AB835486B4C118BCB63D81E22512C6292585A2FAF4FFE4E` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/main_mean_statistic_stability.tsv` | `EC8BC7A6280345DD650C058CDDBDA5FE4C24AC20094590CA3A69ABEFDBD13A81` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/main_p_distribution.tsv` | `B3D11A7141A66452B2557A3C6F24F6E635477FE78895A08820FDA98DB29BE3D4` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/main_reducer_paired_difference.tsv` | `B4F93F45537EB1351A119BC6987910DF8491069DFBE5633D3B1712708BB39879` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/run_manifest.json` | `8F527DD4A2565FDB213648466D5331810F1A178F675E0B236EE63349351B6D2A` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/validated_decision.json` | `5CCAD73301D23E6DEE316441E1C88E77F4630134A2991E359A100F0BCCAA294D` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/validator_pair_rows.tsv.gz` | `E8D4943C8002739FE23FD7771082DC3B6DD7B04C51ACADA1CF9764A6236AC036` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/validator_report.json` | `BAD5CCC6702B4A8E068FCA1A0B998E2EC0BBDCAAD367BF38A9A69A2F577405B4` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/validator_subject_rows.tsv` | `3AC839792465FA28231BA835E1CD8ED271A76758B4E46354221D4BACDEDF15B6` |
| `processed/subject/group/ieeg_p0_7_aggregation_remedy/stage_a_jrobust/validator_world_rows.tsv` | `7AF042B5D0D4EC483918DC1278308601E79675A62FC38752E5129AD54A1F9533` |
