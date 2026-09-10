function checkpoint=run_acc00_p0_6_world(root,arm,seedIndex,outputRoot)
%RUN_ACC00_P0_6_WORLD One atomic/resumable P0-6 checkpoint.
% M-REF route replicates the M02 replay path (base-contract surrogate guards
% 1.5/1.5/0.75 + acc01_compute_phase_metric_optimized); the schedule swap is
% the only mutation versus the historical pool, so a matched seed shares the
% synthesized background bitwise (eta/polarity/phi/background draw outside
% the event loop; additive injection consumes no stream draws).
arguments
    root (1,1) string
    arm (1,1) string
    seedIndex (1,1) double {mustBeInteger,mustBePositive}
    outputRoot (1,1) string = ""
end
maxNumCompThreads(1);
root=string(char(java.io.File(char(root)).getCanonicalPath()));paperPath=fullfile(root,"scripts","ieeg","paper");
addpath(fullfile(root,"scripts","ieeg"));addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(paperPath);
c=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_6_perturbation_contract.json")));
if strlength(outputRoot)==0,outputRoot=fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation");end
outDir=fullfile(outputRoot,"checkpoints",arm);if ~isfolder(outDir),mkdir(outDir),end
final=fullfile(outDir,sprintf("world_%03d.mat",seedIndex));tmp=final+".tmp.mat";
if isfile(final),s=load(final,"checkpoint");assertIdentity(s.checkpoint,c,root,arm,seedIndex);checkpoint=s.checkpoint;return,end
if isfile(tmp),s=load(tmp,"checkpoint");assertIdentity(s.checkpoint,c,root,arm,seedIndex);movefile(tmp,final,"f");checkpoint=s.checkpoint;return,end
switch arm
    case "Mphase_0p5",den="x0p5x";estimator="Mphase";
    case "Mphase_1p5",den="x1p5x";estimator="Mphase";
    case "Mpower_1p5",den="x1p5x";estimator="Mpower";
    case "corner_probe_G06",den="x1p5x";estimator="corner";
    otherwise,error("Unknown arm '%s'.",arm);
end
sch=fullfile(root,c.density_levels.(den).schedule_tsv);assert(isfile(sch),"Missing frozen schedule %s.",sch);
base=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
sim=base;sim.inputs.label_covariates.path=sch;% only the schedule changes; guards stay at the base-contract 1.5/1.5/0.75 reference
started=utcNow();clock=tic;
if estimator=="corner"
    pc=jsondecode(fileread(fullfile(paperPath,"acc00_sim_powercurve_contract.json")));
    gi=6;assert(string(pc.additive.grid_id(gi))=="G06","Powercurve grid_id(6) is not G06.");
    A=double(pc.additive.a_erp_grid_uV(gi+1));% grid carries a leading A=0 reference: G06 = a_erp_grid_uV(7) = 4.61799281479342
    seed=91000000+100000*2+1000*gi+seedIndex;scenario="additive";
    effect=struct("a_erp_uV",A,"lambda_a_uV",0,"kappa",0,"provisional",false);
    world=acc00_sim_generate_world(root,sim,"benchmark",seed,scenario,effect);
    metrics=runPowerRoute(world,sim);
    corner=struct("grid_id",string(pc.additive.grid_id(gi)),"a_erp_uV",A,"gain_intercept_a",double(pc.additive.gain_intercept_a),"gain_slope_b",double(pc.additive.gain_slope_b),"kernel",string(pc.additive.kernel),"scenario_index",2);
else
    seed=20260827+100000+seedIndex;scenario="null";effect=struct();
    world=acc00_sim_generate_world(root,sim,"benchmark",seed,scenario,effect);
    if estimator=="Mpower",metrics=runPowerRoute(world,sim);corner=struct();else,metrics=runPhaseRoute(world,sim);corner=struct();end
end
elapsed=toc(clock);
checkpoint=struct("status","complete","analysis_id",string(c.analysis_id),"arm",arm,"density",den(2:end),"seed_index",seedIndex,"world_seed",seed, ...
    "schedule_path",sch,"schedule_sha256",sha256File(sch),"contract_sha256",sha256File(fullfile(paperPath,"acc00_sim_p0_6_perturbation_contract.json")), ...
    "parent_acc00_contract_sha256",string(c.parent.acc00_sim_contract.sha256),"subject_metrics",metrics.agg,"signflip",metrics.signflip, ...
    "pair_diagnostics",metrics.pairs,"corner",corner,"worker_threads",maxNumCompThreads,"started_at_utc",started,"completed_at_utc",utcNow(), ...
    "timing",struct("elapsed_s",elapsed,"synthesis_s",world.timing.synthesis_s,"peak_memory_bytes",world.timing.peak_memory_bytes));
save(tmp,"checkpoint","-v7.3");movefile(tmp,final,"f");
end
function m=runPowerRoute(world,sim)
rows=repmat(struct("subject","","run","","power",NaN,"phase",0,"n_events",0),numel(world.pairs),1);fit=0;lim=0;
for i=1:numel(world.pairs)
    [events,~,~]=acc01_build_support_and_events(world.pairs(i),world.covariates,sim);pw=acc01_compute_power_metric(world.pairs(i),events,sim);
    rows(i)=struct("subject",world.pairs(i).subject,"run",world.pairs(i).run,"power",pw.slope,"phase",0,"n_events",height(events));
    fit=fit+pw.robustfit_fit_count;lim=lim+pw.robustfit_iteration_limit_count;
end
x=acc01_aggregate_subject_metrics(rows);sf=struct();for r=["R1" "R2"],keep=x.runs==r&all(isfinite(x.matrix(:,[1 2])),2);sf.(char(r))=acc01_exact_signflip_family(x.matrix(keep,:));end
m=struct("agg",x.table,"signflip",reduceSignflip(sf),"pairs",rows,"diag",struct("robustfit_fit_count",fit,"robustfit_iteration_limit_count",lim));
end
function m=runPhaseRoute(world,sim)
% M-REF route: both metrics are computed for every pair so the subject
% inclusion rule (both columns finite, as in acc01_run_metrics_optimized)
% is bit-identical to the historical reference pool; the schedule swap is
% the only mutation versus the historical pool.
rows=repmat(struct("subject","","run","","power",NaN,"phase",NaN,"n_events",0,"domain_fraction",NaN,"surrogate_count",NaN),numel(world.pairs),1);
fit=0;lim=0;
for i=1:numel(world.pairs)
    [pe,he]=acc01_build_support_and_events(world.pairs(i),world.covariates,sim);pw=acc01_compute_power_metric(world.pairs(i),pe,sim);ph=acc01_compute_phase_metric_optimized(world.pairs(i),he,sim,world.stream);
    rows(i)=struct("subject",world.pairs(i).subject,"run",world.pairs(i).run,"power",pw.slope,"phase",ph.difference,"n_events",height(he), ...
        "domain_fraction",ph.surrogate.domain_fraction,"surrogate_count",ph.surrogate.sample_count);
    fit=fit+pw.robustfit_fit_count;lim=lim+pw.robustfit_iteration_limit_count;
end
x=acc01_aggregate_subject_metrics(rows);sf=struct();for r=["R1" "R2"],keep=x.runs==r&all(isfinite(x.matrix(:,[1 2])),2);sf.(char(r))=acc01_exact_signflip_family(x.matrix(keep,:));end
m=struct("agg",x.table,"signflip",reduceSignflip(sf),"pairs",rows,"diag",struct("robustfit_fit_count",fit,"robustfit_iteration_limit_count",lim));
end
function s=reduceSignflip(q),s=struct();for r=["R1" "R2"],s.(char(r))=struct("observed",q.(r).observed,"p_two_sided",q.(r).p_two_sided,"p_maxstat",q.(r).p_maxstat,"studentization_sd",q.(r).studentization_sd,"studentized_observed",q.(r).studentized_observed,"max_observed",q.(r).max_observed);end,end
function assertIdentity(q,c,root,arm,seedIndex) %#ok<INUSD>
assert(string(q.status)=="complete"&&string(q.arm)==arm&&q.seed_index==seedIndex,"P0-6 checkpoint identity mismatch.");
assert(string(q.contract_sha256)==sha256File(fullfile(root,"scripts","ieeg","paper","acc00_sim_p0_6_perturbation_contract.json")),"P0-6 contract hash changed.");
switch arm
    case "Mphase_0p5",den="x0p5x";
    case {"Mphase_1p5","Mpower_1p5","corner_probe_G06"},den="x1p5x";
    otherwise,error("Unknown arm '%s'.",arm);
end
assert(string(q.schedule_sha256)==sha256File(fullfile(root,c.density_levels.(den).schedule_tsv)),"P0-6 schedule hash changed for arm %s.",arm);
end
function s=utcNow(),s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));b=fread(fid,Inf,"*uint8");md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(b,"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
