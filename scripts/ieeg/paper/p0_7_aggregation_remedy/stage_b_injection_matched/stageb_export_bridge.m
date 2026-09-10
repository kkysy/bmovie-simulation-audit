function report=stageb_export_bridge(root,outputRoot,mode)
%STAGEB_EXPORT_BRIDGE MATLAB-only reader for v7.3 tables; no statistics.
arguments
    root (1,1) string
    outputRoot (1,1) string = fullfile(fileparts(mfilename("fullpath")),"output")
    mode (1,1) string {mustBeMember(mode,["smoke" "benchmark"])} = "benchmark"
end
toolDir=fileparts(mfilename("fullpath"));c=jsondecode(fileread(fullfile(toolDir,"stageb_mphase_power_contract.json")));
files=dir(fullfile(outputRoot,"checkpoints","world_*.mat"));[~,ix]=sort({files.name});files=files(ix);assert(~isempty(files),"No checkpoints found.");
pairParts=cell(numel(files)*2,1);subjectParts=pairParts;worldParts=cell(numel(files)*4*2,1);np=0;ns=0;nw=0;
for i=1:numel(files)
    z=load(fullfile(files(i).folder,files(i).name),"checkpoint");q=z.checkpoint;
    assertIdentity(q,c,root,toolDir,mode);
    for a=1:2
        ar=q.arm_results(a);np=np+1;pr=ar.pair_rows;pr.seed_index(:)=q.seed_index;pr.world_seed(:)=q.world_seed;pr.mode(:)=string(q.mode);pr.arm(:)=string(ar.arm);pr.kappa(:)=ar.kappa;
        pairParts{np}=pr(:,["seed_index" "world_seed" "mode" "arm" "kappa" "subject" "run" "session_id" "pair_id" ...
            "corrected_mphase" "observed_ppc" "mean_refit_surrogate_ppc" "n_phase_events" "n_injected_events" ...
            "n_injected_phase_intersection" "epsilon_count" "epsilon_resultant" "injection_draw_count" "injection_stream_seed" "base_fingerprint" "domain_fraction" ...
            "surrogate_count" "refit_surrogate_count" "mphase_surrogate_reducer"]);
        sr=ar.subject_rows;sr.seed_index(:)=q.seed_index;sr.world_seed(:)=q.world_seed;sr.mode(:)=string(q.mode);sr.arm(:)=string(ar.arm);sr.kappa(:)=ar.kappa;ns=ns+1;
        subjectParts{ns}=sr(:,["seed_index" "world_seed" "mode" "arm" "kappa" "subject" "run" "reducer" "Mphase" ...
            "n_pairs_total" "n_pairs_finite" "pair_corrected_skewness" "pair_positive_fraction" "pair_zero_count"]);
        for reducer=["median" "mean"]
            for run=["R1" "R2"]
                sf=ar.inference.(char(reducer)).(char(run));d=ar.run_diagnostics(ar.run_diagnostics.reducer==reducer&ar.run_diagnostics.run==run,:);nw=nw+1;
                worldParts{nw}=table(q.seed_index,q.world_seed,string(q.mode),string(ar.arm),ar.kappa,reducer,run,double(sf.n_subjects),double(sf.observed(2)),double(sf.p_two_sided(2)),double(sf.studentization_sd(2)),double(d.subject_statistic_skewness),double(d.n_positive),double(d.n_negative),double(d.n_zero),double(d.positive_fraction),string(ar.metric_stream_initial_state_hash),string(q.base_stream_state_hash),double(q.truth.A_G06_uV),string(q.truth.grid_id),'VariableNames',["seed_index" "world_seed" "mode" "arm" "kappa" "reducer" "run" "n_subjects" "observed" "p_two_sided" "studentization_sd" "subject_statistic_skewness" "n_positive" "n_negative" "n_zero" "positive_fraction" "metric_stream_initial_state_hash" "base_stream_state_hash" "truth_a_g06_uV" "truth_grid_id"]);
            end
        end
    end
end
pairs=vertcat(pairParts{1:np});subjects=vertcat(subjectParts{1:ns});worlds=vertcat(worldParts{1:nw});
writetable(pairs,fullfile(outputRoot,"validator_pair_rows.tsv"),"FileType","text","Delimiter",char(9));
writetable(subjects,fullfile(outputRoot,"validator_subject_rows.tsv"),"FileType","text","Delimiter",char(9));
writetable(worlds,fullfile(outputRoot,"validator_world_rows.tsv"),"FileType","text","Delimiter",char(9));
report=struct("status","PASS","mode",mode,"matched_checkpoints",numel(files),"pair_rows",height(pairs),"subject_rows",height(subjects),"world_rows",height(worlds));
fid=fopen(fullfile(outputRoot,"bridge_report.json"),"w");assert(fid>=0);cl=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",jsonencode(report,"PrettyPrint",true));disp(jsonencode(report,"PrettyPrint",true));
end
function assertIdentity(q,c,root,toolDir,mode)
assert(string(q.analysis_id)==string(c.analysis_id)&&string(q.mode)==mode&&q.status=="complete"&&numel(q.arm_results)==2,"Bridge identity mismatch.");
assert(string(q.schedule_sha256)==string(c.schedule.sha256),"Bridge schedule mismatch.");
assert(string(q.contract_sha256)==sha256File(fullfile(toolDir,"stageb_mphase_power_contract.json")),"Bridge contract hash mismatch.");
assert(string(q.parent_contract_sha256)==string(c.parent_contract.sha256),"Bridge parent hash mismatch.");
assert(string(q.generator_sha256)==string(c.local_generator_sha256),"Bridge generator hash mismatch.");
assert(all(string({q.arm_results.arm})==["P0" "P4"]),"Bridge arm mismatch.");
assert(isfile(fullfile(root,string(c.schedule.path))),"Schedule missing.");
end
function h=sha256File(path)
fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance("SHA-256");
while true,b=fread(fid,1048576,"*uint8");if isempty(b),break,end;md.update(typecast(b,"int8"));end
h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end
