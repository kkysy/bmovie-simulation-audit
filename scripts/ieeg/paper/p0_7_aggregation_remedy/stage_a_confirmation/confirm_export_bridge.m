% Release port of the analysis-workspace tool of the same name; see README in this folder.
function report = confirm_export_bridge(root, outputRoot, mode)
%CONFIRM_EXPORT_BRIDGE Zero-statistics bridge for MATLAB v7.3 tables/strings.
arguments
    root (1,1) string = string(pwd) % release port: run from the repository root
        % (the workspace original resolved the root from the tools/ layout) %#ok<INUSA>
    outputRoot (1,1) string = fullfile(fileparts(mfilename("fullpath")),"output")
    mode (1,1) string {mustBeMember(mode,["smoke" "benchmark"])} = "benchmark"
end
toolDir=fileparts(mfilename("fullpath"));
c=jsondecode(fileread(fullfile(toolDir,"confirm_mphase_reducer_contract.json")));
files=dir(fullfile(outputRoot,"checkpoints",string(c.arm.id),"world_*.mat"));
[~,ix]=sort({files.name});files=files(ix);
assert(~isempty(files),"No checkpoints found for bridge export.");
pairParts=cell(numel(files),1);subjectParts=cell(numel(files),1);worldParts=cell(numel(files)*4,1);w=0;
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);z=load(path,"checkpoint");q=z.checkpoint;
    assert(string(q.analysis_id)==string(c.analysis_id) && string(q.arm)==string(c.arm.id) ...
        && string(q.mode)==mode,"Bridge checkpoint identity mismatch: %s",path);
    p=q.pair_rows;p.seed_index(:)=q.seed_index;p.world_seed(:)=q.world_seed;p.mode(:)=string(q.mode);
    p=p(:,["seed_index" "world_seed" "mode" "subject" "run" "session_id" "pair_id" ...
        "corrected_mphase" "observed_ppc" "mean_refit_surrogate_ppc" "n_phase_events" ...
        "domain_fraction" "surrogate_count" "refit_surrogate_count" "mphase_surrogate_reducer"]);
    pairParts{i}=p;
    s=q.subject_rows;s.seed_index(:)=q.seed_index;s.world_seed(:)=q.world_seed;s.mode(:)=string(q.mode);
    s=s(:,["seed_index" "world_seed" "mode" "subject" "run" "reducer" "Mphase" ...
        "n_pairs_total" "n_pairs_finite" "pair_corrected_skewness" "pair_positive_fraction" "pair_zero_count"]);
    subjectParts{i}=s;
    for reducer=["median" "mean"]
        for run=["R1" "R2"]
            sf=q.inference.(char(reducer)).(char(run));
            d=q.run_diagnostics(q.run_diagnostics.reducer==reducer & q.run_diagnostics.run==run,:);
            assert(height(d)==1,"Missing run diagnostic.");w=w+1;
            worldParts{w}=table(q.seed_index,q.world_seed,string(q.mode),reducer,run, ...
                double(sf.n_subjects),double(sf.observed(2)),double(sf.p_two_sided(2)), ...
                double(sf.studentization_sd(2)),double(d.subject_statistic_skewness), ...
                double(d.n_positive),double(d.n_negative),double(d.n_zero),double(d.positive_fraction), ...
                'VariableNames',["seed_index" "world_seed" "mode" "reducer" "run" "n_subjects" ...
                "observed" "p_two_sided" "studentization_sd" "subject_statistic_skewness" ...
                "n_positive" "n_negative" "n_zero" "positive_fraction"]);
        end
    end
end
pairRows=vertcat(pairParts{:});subjectRows=vertcat(subjectParts{:});worldRows=vertcat(worldParts{1:w});
writetable(pairRows,fullfile(outputRoot,"validator_pair_rows.tsv"),"FileType","text","Delimiter",char(9));
writetable(subjectRows,fullfile(outputRoot,"validator_subject_rows.tsv"),"FileType","text","Delimiter",char(9));
writetable(worldRows,fullfile(outputRoot,"validator_world_rows.tsv"),"FileType","text","Delimiter",char(9));
report=struct("status","PASS","mode",mode,"checkpoints",numel(files), ...
    "pair_rows",height(pairRows),"subject_rows",height(subjectRows),"world_rows",height(worldRows));
fid=fopen(fullfile(outputRoot,"bridge_report.json"),"w");assert(fid>=0);cl=onCleanup(@()fclose(fid));
fprintf(fid,"%s\n",jsonencode(report,"PrettyPrint",true));disp(jsonencode(report,"PrettyPrint",true));
end
