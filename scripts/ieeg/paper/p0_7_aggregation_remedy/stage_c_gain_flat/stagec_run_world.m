function checkpoint = stagec_run_world(root, seedIndex, mode, outputRoot)
%STAGEC_RUN_WORLD Run one Stage C gain-flattened P0F world.
arguments
    root (1,1) string
    seedIndex (1,1) double {mustBeInteger,mustBePositive}
    mode (1,1) string {mustBeMember(mode,["smoke" "benchmark"])} = "benchmark"
    outputRoot (1,1) string = fullfile(fileparts(mfilename("fullpath")),"output")
end
maxNumCompThreads(1);root=canonicalPath(root);toolDir=fileparts(mfilename("fullpath"));
addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(toolDir);
contractPath=fullfile(toolDir,"stagec_mphase_gainflat_contract.json");
parentPath=fullfile(root,"scripts","ieeg","acc00_sim_contract.json");
c=jsondecode(fileread(contractPath));base=jsondecode(fileread(parentPath));
validateContract(c,base,root,parentPath,seedIndex);
worldSeed=double(c.seed_namespace.base)+seedIndex;
outDir=fullfile(outputRoot,"checkpoints");if ~isfolder(outDir),mkdir(outDir);end
final=fullfile(outDir,sprintf("world_%03d.mat",seedIndex));tmp=final+".tmp.mat";
if isfile(final)
    assert(~isfile(tmp),"Complete checkpoint has stale tmp residue.");z=load(final,"checkpoint");
    assertIdentity(z.checkpoint,c,contractPath,parentPath,root,seedIndex,worldSeed,mode);checkpoint=z.checkpoint;return
end
if isfile(tmp)
    z=load(tmp,"checkpoint");assertIdentity(z.checkpoint,c,contractPath,parentPath,root,seedIndex,worldSeed,mode);
    movefile(tmp,final,"f");checkpoint=z.checkpoint;return
end
started=utcNow();clock=tic;world=stagec_generate_gainflat_world(root,c,mode,worldSeed);
arm="P0F";pairs=world.pairs.P0F;metricStream=RandStream("Threefry","Seed",worldSeed);
metricStream.State=world.base_stream_state;pairRows=repmat(pairTemplate(),numel(pairs),1);
for i=1:numel(pairs)
    pair=pairs(i);[~,events]=acc01_build_support_and_events(pair,world.covariates,base);
    % Contract intersection_rule (same as Stage B): phaseGood & mx is a subset of
    % the source_complete injection set and its wider filter support implies
    % that every selected event has a valid [0.10,0.50] injection window.
    assert(all(events.source_complete==1),"Phase event is not source_complete.");
    assert(all(double(events.label_time_s)-0.10>=-1e-12 & ...
        double(events.label_time_s)+0.50<=double(pair.T)+1e-12), ...
        "Phase event has no valid injection window.");
    ph=acc01_compute_phase_metric_optimized(pair,events,base,metricStream);
    pairRows(i)=struct("subject",string(pair.subject),"run",string(pair.run), ...
        "session_id",string(pair.session_id),"pair_id",string(pair.pair_id), ...
        "corrected_mphase",double(ph.difference),"observed_ppc",double(ph.observed_ppc), ...
        "mean_refit_surrogate_ppc",double(ph.surrogate_ppc),"n_phase_events",height(events), ...
        "n_injected_events",double(pair.n_injected_events), ...
        "n_injected_phase_intersection",height(events),"epsilon_count",numel(pair.epsilon_event_rad), ...
        "epsilon_resultant",epsilonResultant(pair.epsilon_event_rad),"injection_draw_count",double(pair.injection_draw_count), ...
        "injection_stream_seed",double(pair.injection_stream_seed), ...
        "base_fingerprint",string(pair.base_fingerprint),"domain_fraction",double(ph.surrogate.domain_fraction), ...
        "surrogate_count",double(ph.surrogate.sample_count),"refit_surrogate_count",double(ph.refit_surrogate_count), ...
        "mphase_surrogate_reducer",string(ph.diagnostics.reducer));
end
[subjectRows,inference,runDiagnostics]=reduceBoth(pairRows);
armResult=struct("status","complete","arm",arm,"kappa",0.0,"gain_slope_b",0.0, ...
    "pair_rows",struct2table(pairRows),"subject_rows",subjectRows,"inference",inference, ...
    "run_diagnostics",runDiagnostics,"metric_stream_initial_state_hash",stateHash(world.base_stream_state), ...
    "metric_stream_final_state_hash",stateHash(metricStream.State));
checkpoint=struct("status","complete","analysis_id",string(c.analysis_id),"mode",mode, ...
    "seed_index",seedIndex,"world_seed",worldSeed,"arms",[arm], ...
    "contract_sha256",sha256File(contractPath),"parent_contract_sha256",sha256File(parentPath), ...
    "schedule_sha256",sha256File(fullfile(root,string(c.schedule.path))), ...
    "generator_sha256",sha256File(fullfile(toolDir,"stagec_generate_gainflat_world.m")), ...
    "input_hashes",inputHashes(base.inputs),"truth",world.truth,"arm_results",armResult, ...
    "base_stream_state_hash",stateHash(world.base_stream_state),"started_at_utc",started, ...
    "completed_at_utc",utcNow(),"worker_threads",maxNumCompThreads,"matlab_version",string(version), ...
    "timing",struct("elapsed_s",toc(clock),"synthesis_s",world.timing.synthesis_s, ...
    "peak_memory_bytes",world.timing.peak_memory_bytes));
save(tmp,"checkpoint","-v7.3");movefile(tmp,final,"f");
end

function [t,inference,diagTable]=reduceBoth(rows)
subject=string({rows.subject})';run=string({rows.run})';value=double([rows.corrected_mphase])';
keys=unique(table(subject,run),'rows','sorted');reducers=["median" "mean"];out=cell(height(keys)*2,9);k=0;
for i=1:height(keys)
    keep=subject==keys.subject(i)&run==keys.run(i);v=value(keep);vf=v(isfinite(v));
    assert(~isempty(vf),"Subject has no finite Mphase pairs.");pairSkew=finiteSkewness(vf);
    for reducer=reducers,k=k+1;if reducer=="median",statistic=median(vf);else,statistic=mean(vf);end
        out(k,:)={keys.subject(i),keys.run(i),reducer,statistic,numel(v),numel(vf),pairSkew,mean(vf>0),nnz(vf==0)};
    end
end
t=cell2table(out,'VariableNames',["subject" "run" "reducer" "Mphase" "n_pairs_total" ...
    "n_pairs_finite" "pair_corrected_skewness" "pair_positive_fraction" "pair_zero_count"]);
inference=struct();diagRows=cell(4,10);diagRow=0;
for reducer=reducers
    for r=["R1" "R2"]
        z=sortrows(t(t.reducer==reducer&t.run==r,:),"subject");x=double(z.Mphase);
        if isempty(x),sf=emptySignflip();else,q=acc01_exact_signflip_family([zeros(numel(x),1),x]);sf=reduceSignflip(q,string(z.subject));end
        inference.(char(reducer)).(char(r))=sf;
        diagRow=diagRow+1;diagRows(diagRow,:)={reducer,r,numel(x),finiteSkewness(x),nnz(x>0),nnz(x<0),nnz(x==0), ...
            safeMean(x>0),safeMean(x),safeStd(x)};
    end
end
diagRows=diagRows(1:diagRow,:);diagTable=cell2table(diagRows,'VariableNames',["reducer" "run" "n_subjects" "subject_statistic_skewness" ...
    "n_positive" "n_negative" "n_zero" "positive_fraction" "subject_mean" "subject_sample_sd"]);
end
function s=reduceSignflip(q,subjects),s=struct("subjects",subjects,"n_subjects",numel(subjects),"observed",q.observed,"p_two_sided",q.p_two_sided,"p_maxstat",q.p_maxstat,"studentization_sd",q.studentization_sd,"studentized_observed",q.studentized_observed,"max_observed",q.max_observed);end
function s=emptySignflip(),s=struct("subjects",strings(0,1),"n_subjects",0,"observed",[NaN NaN],"p_two_sided",[NaN NaN],"p_maxstat",[NaN NaN],"studentization_sd",[NaN NaN],"studentized_observed",[NaN NaN],"max_observed",NaN);end
function q=pairTemplate(),q=struct("subject","","run","","session_id","","pair_id","","corrected_mphase",NaN,"observed_ppc",NaN,"mean_refit_surrogate_ppc",NaN,"n_phase_events",0,"n_injected_events",0,"n_injected_phase_intersection",0,"epsilon_count",0,"epsilon_resultant",NaN,"injection_draw_count",0,"injection_stream_seed",0,"base_fingerprint","","domain_fraction",NaN,"surrogate_count",0,"refit_surrogate_count",0,"mphase_surrogate_reducer","");end
function validateContract(c,b,root,parentPath,idx)
assert(string(c.parent_contract.sha256)==sha256File(parentPath)&&string(b.analysis_id)==string(c.parent_contract.required_analysis_id),"Parent identity mismatch.");
assert(idx>=double(c.seed_namespace.seed_index(1))&&idx<=double(c.seed_namespace.seed_index(2)),"seed_index outside 1:60.");
assert(c.injection.A_G06_uV==4.61799281479342&&string(c.injection.grid_id)=="G06"&&c.injection.frequency_hz==6,"Injection constants changed.");
assert(c.injection.gain_intercept_a==0&&c.injection.gain_slope_b==0,"Gain must be flattened (b=0).");
assert(b.phase.refit_surrogates.J==32&&string(b.phase.mphase_surrogate_reducer)=="mean","Production Mphase changed.");
assert(sha256File(fullfile(root,string(c.schedule.path)))==string(c.schedule.sha256),"Schedule identity mismatch.");
end
function assertIdentity(q,c,contractPath,parentPath,root,idx,seed,mode)
assert(string(q.status)=="complete"&&string(q.analysis_id)==string(c.analysis_id)&&q.seed_index==idx&&q.world_seed==seed&&string(q.mode)==mode,"Checkpoint identity mismatch.");
assert(isequal(string(q.arms),["P0F"])&&numel(q.arm_results)==1&&string(q.arm_results.status)=="complete","Incomplete arm.");
assert(string(q.contract_sha256)==sha256File(contractPath)&&string(q.parent_contract_sha256)==sha256File(parentPath),"Contract hash mismatch.");
assert(string(q.schedule_sha256)==sha256File(fullfile(root,string(c.schedule.path))),"Schedule hash mismatch.");
assert(string(q.generator_sha256)==string(c.local_generator_sha256),"Stage-C generator hash mismatch.");
end
function x=epsilonResultant(v),if isempty(v),x=NaN;else,x=abs(mean(exp(1i*double(v(:)))));end,end
function x=finiteSkewness(v),if numel(v)<3||all(v==v(1)),x=NaN;else,x=skewness(v,0);end,end
function x=safeMean(v),if isempty(v),x=NaN;else,x=mean(v);end,end
function x=safeStd(v),if numel(v)<2,x=NaN;else,x=std(v,0);end,end
function h=stateHash(s),h=sha256Text(sprintf('%u,',uint32(s)));end
function h=inputHashes(inputs),h=struct();for n=string(fieldnames(inputs))',h.(n)=upper(string(inputs.(n).sha256));end,end
function s=utcNow(),s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));end
function p=canonicalPath(p),p=string(char(java.io.File(char(p)).getCanonicalPath()));end
function h=sha256Text(t),md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(uint8(char(t)),"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance("SHA-256");while true,b=fread(fid,1048576,"*uint8");if isempty(b),break,end;md.update(typecast(b,"int8"));end;h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
