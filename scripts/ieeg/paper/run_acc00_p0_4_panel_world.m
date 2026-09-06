function checkpoint = run_acc00_p0_4_panel_world(root, opaqueId, seedIndex, mode, outputRoot)
%RUN_ACC00_P0_4_PANEL_WORLD One atomic/resumable pure-null panel checkpoint.
arguments
    root (1,1) string
    opaqueId (1,1) string
    seedIndex (1,1) double {mustBeInteger,mustBePositive}
    mode (1,1) string {mustBeMember(mode,["smoke" "synthetic"])} = "synthetic"
    outputRoot (1,1) string = ""
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));paperPath=fullfile(root,"scripts","ieeg","paper");addpath(fullfile(root,"scripts","ieeg"));addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(paperPath);
contract=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_4_panel_contract.json")));seal=jsondecode(fileread(fullfile(paperPath,string(contract.panel.member_semantics_seal))));assert(any(string(contract.panel.opaque_members)==opaqueId),"Unknown opaque P0-4 member.");assert(seedIndex<=double(contract.panel.null_worlds),"seed_index exceeds frozen panel size.");
member=seal.opaque_to_member.(matlab.lang.makeValidName(char(opaqueId)));seed=worldSeed(contract,seedIndex);if strlength(outputRoot)==0,outputRoot=fullfile(root,string(contract.outputs.root));end
outDir=fullfile(outputRoot,"checkpoints",opaqueId);if ~isfolder(outDir),mkdir(outDir),end;final=fullfile(outDir,sprintf("world_%03d.mat",seedIndex));tmp=final+".tmp.mat";
if isfile(final)
    assert(~isfile(tmp),"P0-4 final checkpoint has stale tmp residue; stop rather than overwrite.");s=load(final,"checkpoint");assertIdentity(s.checkpoint,contract,root,opaqueId,seedIndex,seed,mode);checkpoint=s.checkpoint;return
end
if isfile(tmp)
    s=load(tmp,"checkpoint");assertIdentity(s.checkpoint,contract,root,opaqueId,seedIndex,seed,mode);movefile(tmp,final,"f");checkpoint=s.checkpoint;return
end
started=utcNow();clock=tic;base=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));world=acc00_sim_generate_world(root,base,legacyMode(mode),seed,"null",struct());
if string(member.estimator)=="Mpower_aperiodic_adjusted"
    [subjectMetrics,signflip,pairDiagnostics,domain,diagnostics]=runPower(world,base);
    checkpoint=baseCheckpoint("Mpower",subjectMetrics,signflip,pairDiagnostics,domain,diagnostics);
else
    sim=base;sim.surrogate.domain_guard_edge_lower_s=double(member.guard.edge_lower_s);sim.surrogate.domain_guard_edge_upper_s=double(member.guard.edge_upper_s);sim.surrogate.domain_guard_bad_dilation_s=double(member.guard.bad_dilation_s);
    [subjectMetrics,mphaseShape,pairDiagnostics,domain,diagnostics]=runPhase(world,sim,string(member.crossfit_coefficient_mode));
    checkpoint=baseCheckpoint("Mphase",subjectMetrics,struct(),pairDiagnostics,domain,diagnostics);checkpoint.mphase_shape=mphaseShape;
end
checkpoint.status="complete";checkpoint.analysis_id=string(contract.analysis_id);checkpoint.arm=checkpoint.arm;checkpoint.member_id=string(member.member_id);checkpoint.opaque_id=opaqueId;checkpoint.world_seed=seed;checkpoint.seed_index=seedIndex;checkpoint.scenario_index=double(contract.panel.scenario_index);checkpoint.effect_index=double(contract.panel.effect_index);checkpoint.contract_sha256=sha256File(fullfile(paperPath,"acc00_sim_p0_4_panel_contract.json"));checkpoint.parent_p0_3_contract_sha256=string(contract.parent_p0_3.sha256);checkpoint.input_hashes=inputHashes(contract.inputs);checkpoint.started_at_utc=started;checkpoint.completed_at_utc=utcNow();checkpoint.timing=struct("elapsed_s",toc(clock),"synthesis_s",world.timing.synthesis_s,"peak_memory_bytes",world.timing.peak_memory_bytes);checkpoint.mutation=string(member.mutation);checkpoint.guard=guardOf(member);checkpoint.semantic_identity=semanticIdentity(contract,member);checkpoint.mode=mode;
save(tmp,"checkpoint","-v7.3");movefile(tmp,final,"f");
    function c=baseCheckpoint(arm,metrics,sf,pairs,domainSummary,diag)
        c=struct("arm",arm,"subject_metrics",metrics,"signflip",sf,"pair_diagnostics",pairs,"domain",domainSummary,"diagnostics",diag,"world_summary",struct("n_pairs",numel(world.pairs),"null_only",true));
    end
end
function [agg,sf,pairs,domain,diag]=runPower(world,sim)
rows=repmat(struct("subject","","run","","power",NaN,"phase",0),numel(world.pairs),1);pairs=repmat(struct("subject","","run","","slope",NaN,"n_events",NaN),numel(world.pairs),1);fit=0;lim=0;
for i=1:numel(world.pairs)
    [events,~,~]=acc01_build_support_and_events(world.pairs(i),world.covariates,sim);pw=acc01_compute_power_metric(world.pairs(i),events,sim);rows(i)=struct("subject",world.pairs(i).subject,"run",world.pairs(i).run,"power",pw.slope,"phase",0);pairs(i)=struct("subject",world.pairs(i).subject,"run",world.pairs(i).run,"slope",pw.slope,"n_events",height(events));fit=fit+pw.robustfit_fit_count;lim=lim+pw.robustfit_iteration_limit_count;
end
x=acc01_aggregate_subject_metrics(rows);agg=x.table;sf=struct();for r=["R1" "R2"],q=acc01_exact_signflip_family(x.matrix(x.runs==r,:));sf.(char(r))=reduceSignflip(q);end;domain=struct("mode","not_applicable_Mpower","quantized_points",NaN,"interval_count",NaN,"domain_fraction",NaN);diag=struct("robustfit_fit_count",fit,"robustfit_iteration_limit_count",lim);
end
function [agg,shape,pairs,domain,diag]=runPhase(world,sim,coefficientMode)
rows=repmat(struct("subject","","run","","power",0,"phase",NaN),numel(world.pairs),1);pairs=repmat(struct("subject","","run","","domain_fraction",NaN,"surrogate_count",NaN,"refit_surrogate_count",NaN,"coefficient_mode","","observed_coefficient_solves",NaN,"shift_epoch_coefficient_solves",NaN),numel(world.pairs),1);
for i=1:numel(world.pairs)
    [~,events]=acc01_build_support_and_events(world.pairs(i),world.covariates,sim);ph=acc01_compute_phase_metric_optimized(world.pairs(i),events,sim,world.stream,"CrossfitCoefficientMode",coefficientMode);rows(i)=struct("subject",world.pairs(i).subject,"run",world.pairs(i).run,"power",0,"phase",ph.difference);d=ph.diagnostics.observed_path.crossfit;pairs(i)=struct("subject",world.pairs(i).subject,"run",world.pairs(i).run,"domain_fraction",ph.surrogate.domain_fraction,"surrogate_count",ph.surrogate.sample_count,"refit_surrogate_count",ph.refit_surrogate_count,"coefficient_mode",string(d.coefficient_mode),"observed_coefficient_solves",d.observed_coefficient_solves,"shift_epoch_coefficient_solves",d.shift_epoch_coefficient_solves);
end
x=acc01_aggregate_subject_metrics(rows);agg=x.table;shape=struct();for r=["R1" "R2"],keep=x.runs==r&isfinite(x.matrix(:,2));shape.(char(r))=acc00_p0_4_shape_for_run(x.subjects(keep),x.matrix(keep,2));end
domain=struct("mode","empirical_surrogate","quantized_points",sum([pairs.surrogate_count]),"interval_count",NaN,"domain_fraction",mean([pairs.domain_fraction],"omitnan"));diag=struct("coefficient_mode",coefficientMode,"observed_coefficient_solves",sum([pairs.observed_coefficient_solves]),"shift_epoch_coefficient_solves",sum([pairs.shift_epoch_coefficient_solves]));
end
function s=reduceSignflip(q),s=struct("observed",q.observed,"p_two_sided",q.p_two_sided,"p_maxstat",q.p_maxstat,"studentization_sd",q.studentization_sd,"studentized_observed",q.studentized_observed,"max_observed",q.max_observed);end
function g=guardOf(member),if isfield(member,"guard")&&isstruct(member.guard),g=member.guard;else,g=struct();end,end
function s=semanticIdentity(c,m),s=struct("parent_p0_3_contract_sha256",string(c.parent_p0_3.sha256),"input_hashes",inputHashes(c.inputs),"sampling",c.sampling,"background",c.background,"windows_s",c.windows_s,"power",c.power,"phase",c.phase,"sign_flip",c.sign_flip,"estimator",string(m.estimator));end
function s=worldSeed(c,i),s=double(c.panel.base_seed)+100000*double(c.panel.scenario_index)+1000*double(c.panel.effect_index)+i;end
function m=legacyMode(m),if m=="smoke",m="smoke";else,m="benchmark";end,end
function assertIdentity(q,c,root,id,si,seed,mode),assert(string(q.status)=="complete"&&string(q.analysis_id)==string(c.analysis_id)&&string(q.opaque_id)==id&&q.seed_index==si&&q.world_seed==seed&&string(q.mode)==mode&&(mode=="smoke"||double(q.world_summary.n_pairs)==351)&&string(q.contract_sha256)==sha256File(fullfile(root,"scripts","ieeg","paper","acc00_sim_p0_4_panel_contract.json")),"P0-4 checkpoint identity mismatch (mode/scale mismatch; stop rather than resume or overwrite).");end
function h=inputHashes(x),h=struct();for k=string(fieldnames(x))',h.(k)=string(x.(k).sha256);end,end
function s=utcNow(),s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));b=fread(fid,Inf,"*uint8");md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(b,"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
