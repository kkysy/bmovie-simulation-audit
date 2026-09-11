function report = jrobust_export_bridge(root, outputRoot, runKind)
%JROBUST_EXPORT_BRIDGE Export v7.3 table/string fields to flat TSV files.
arguments
    root (1,1) string
    outputRoot (1,1) string
    runKind (1,1) string {mustBeMember(runKind,["smoke" "formal"])}
end
root = string(char(java.io.File(char(root)).getCanonicalPath()));
toolDir = fileparts(mfilename("fullpath"));
c = jsondecode(fileread(fullfile(toolDir,"confirm_mphase_jrobust_contract.json")));
assert(isfile(fullfile(root,string(c.parent_contract.path))),"Parent contract is missing.");
files = dir(fullfile(outputRoot,"checkpoints",string(c.world.arm_id),"world_*.mat"));
[~,order] = sort({files.name});
files = files(order);
assert(~isempty(files),"No J-robust checkpoints found.");
pairParts = cell(numel(files)*2,1);
subjectParts = cell(numel(files)*2,1);
worldParts = cell(numel(files)*8,1);
np = 0;
ns = 0;
nw = 0;
for i = 1:numel(files)
    path = fullfile(files(i).folder,files(i).name);
    z = load(path,"checkpoint");
    q = z.checkpoint;
    assert(string(q.status)=="complete" && string(q.analysis_id)==string(c.analysis_id) ...
        && string(q.run_kind)==runKind && string(q.generator_mode)==string(c.world.generator_mode), ...
        "Bridge checkpoint identity mismatch: %s",path);
    assert(isequal(string(q.arm_order),string(c.arms.order)) && numel(q.arm_results)==2, ...
        "Bridge matched-arm identity mismatch: %s",path);
    for a = 1:numel(q.arm_results)
        arm = q.arm_results(a);
        np = np+1;
        p = arm.pair_rows;
        p.seed_index(:) = q.seed_index;
        p.world_seed(:) = q.world_seed;
        p.run_kind(:) = string(q.run_kind);
        p.generator_mode(:) = string(q.generator_mode);
        p.background_world_sha256(:) = string(q.background_world_sha256);
        pairParts{np} = p(:,["seed_index" "world_seed" "run_kind" "generator_mode" ...
            "background_world_sha256" "arm" "J_contract" "subject" "run" "session_id" "pair_id" ...
            "corrected_mphase" "observed_ppc" "mean_refit_surrogate_ppc" "n_phase_events" ...
            "quantized_feasible_count" "domain_fraction" "surrogate_count" "surrogate_family_mode" ...
            "refit_selection_mode" "refit_surrogate_count" "mphase_surrogate_reducer" ...
            "candidate_delta_sha256" "refit_delta_sha256" "pair_stream_initial_state_sha256" ...
            "pair_stream_final_state_sha256" "pair_signal_sha256" "pair_elapsed_s" "surrogate_ppc_s"]);

        ns = ns+1;
        s = arm.subject_rows;
        s.seed_index(:) = q.seed_index;
        s.world_seed(:) = q.world_seed;
        s.run_kind(:) = string(q.run_kind);
        s.arm(:) = string(arm.arm);
        s.J_contract(:) = double(arm.J_contract);
        subjectParts{ns} = s(:,["seed_index" "world_seed" "run_kind" "arm" "J_contract" ...
            "subject" "run" "reducer" "Mphase" "n_pairs_total" "n_pairs_finite" ...
            "pair_corrected_skewness" "pair_positive_fraction" "pair_zero_count"]);

        for reducer = ["median" "mean"]
            for run = ["R1" "R2"]
                sf = arm.inference.(char(reducer)).(char(run));
                d = arm.run_diagnostics(arm.run_diagnostics.reducer==reducer & arm.run_diagnostics.run==run,:);
                assert(height(d)==1,"Missing run diagnostic in %s.",path);
                nw = nw+1;
                worldParts{nw} = table(q.seed_index,q.world_seed,string(q.run_kind),string(arm.arm), ...
                    double(arm.J_contract),reducer,run,double(sf.n_subjects),double(sf.observed(2)), ...
                    double(sf.p_two_sided(2)),double(sf.studentization_sd(2)), ...
                    double(d.subject_statistic_skewness),double(d.n_positive),double(d.n_negative), ...
                    double(d.n_zero),double(d.positive_fraction),string(arm.effective_contract_sha256), ...
                    string(q.background_world_sha256), ...
                    'VariableNames',["seed_index" "world_seed" "run_kind" "arm" "J_contract" ...
                    "reducer" "run" "n_subjects" "observed" "p_two_sided" "studentization_sd" ...
                    "subject_statistic_skewness" "n_positive" "n_negative" "n_zero" ...
                    "positive_fraction" "effective_contract_sha256" "background_world_sha256"]);
            end
        end
    end
end
pairRows = vertcat(pairParts{1:np});
subjectRows = vertcat(subjectParts{1:ns});
worldRows = vertcat(worldParts{1:nw});
writetable(pairRows,fullfile(outputRoot,"validator_pair_rows.tsv"),"FileType","text","Delimiter",char(9));
writetable(subjectRows,fullfile(outputRoot,"validator_subject_rows.tsv"),"FileType","text","Delimiter",char(9));
writetable(worldRows,fullfile(outputRoot,"validator_world_rows.tsv"),"FileType","text","Delimiter",char(9));
report = struct("status","PASS","run_kind",runKind,"checkpoints",numel(files), ...
    "pair_rows",height(pairRows),"subject_rows",height(subjectRows),"world_rows",height(worldRows));
writeJsonAtomic(fullfile(outputRoot,"bridge_report.json"),report);
disp(jsonencode(report,"PrettyPrint",true));
end

function writeJsonAtomic(path,value)
tmp=path+".tmp";
fid=fopen(tmp,"w");assert(fid>=0,"Cannot write %s",tmp);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,"%s\n",jsonencode(value,"PrettyPrint",true));
clear cleanup
movefile(tmp,path,"f");
end
