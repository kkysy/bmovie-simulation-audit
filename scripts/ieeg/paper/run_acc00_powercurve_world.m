function checkpoint = run_acc00_powercurve_world(root, gridIndex, seedIndex, version)
% Execute one paper-only additive Mpower world with resumable identity checks.
arguments
    root (1,1) string
    gridIndex (1,1) double {mustBeInteger,mustBePositive}
    seedIndex (1,1) double {mustBeInteger,mustBePositive}
    version (1,1) string = "v2"
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
paperPath=fullfile(root,"scripts","ieeg","paper"); addpath(fullfile(root,"scripts","ieeg","acc01"));
paper=jsondecode(fileread(fullfile(paperPath,"acc00_sim_powercurve_contract.json")));
base=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
assert(gridIndex<=numel(paper.additive.grid_worlds),"gridIndex outside frozen grid");
assert(seedIndex<=double(paper.additive.grid_worlds(gridIndex)),"seedIndex exceeds frozen cell size");
seed=double(paper.seed.paper_base_seed)+100000*double(paper.seed.scenario_index)+1000*gridIndex+seedIndex;
outDir=fullfile(root,"processed","subject","group","ieeg_methods_paper","p0_2_sim_powercurve","profiling",char(version));
if ~isfolder(outDir),mkdir(outDir);end
final=fullfile(outDir,sprintf("world_g%02d_n%03d.mat",gridIndex,seedIndex));
% Resume-skip: a completed world with fully verified identity is returned
% as-is; identity mismatch stops the run rather than silently overwriting.
if isfile(final)
    s=load(final,"checkpoint"); assertIdentity(s.checkpoint,paper,root,gridIndex,seedIndex,seed);
    checkpoint=s.checkpoint; return;
end
started=utcNow(); clock=tic; pid=feature("getpid"); memMat=0; memProc=0;
effect=struct("a_erp_uV",double(paper.additive.a_erp_grid_uV(gridIndex+1)),"provisional",false);
w=acc00_sim_generate_world(root,base,"benchmark",seed,"additive",effect); [memMat,memProc]=sampleMemory(memMat,memProc,pid);
rows=repmat(struct("subject","","run","","power",NaN,"phase",NaN),numel(w.pairs),1);
diagCount=0; fitCount=0; nPower=0; nPhase=0;
for p=1:numel(w.pairs)
    [ev,ph,evDiag]=acc01_build_support_and_events(w.pairs(p),w.covariates,base);
    q=acc01_compute_power_metric(w.pairs(p),ev,base);
    rows(p)=struct("subject",w.pairs(p).subject,"run",w.pairs(p).run,"power",q.slope,"phase",0);
    diagCount=diagCount+q.robustfit_iteration_limit_count; fitCount=fitCount+q.robustfit_fit_count;
    nPower=nPower+evDiag.n_power; nPhase=nPhase+evDiag.n_phase;
    if mod(p,10)==0 || p==numel(w.pairs), [memMat,memProc]=sampleMemory(memMat,memProc,pid); end
end
agg=acc01_aggregate_subject_metrics(rows);
sfR1=acc01_exact_signflip_family(agg.matrix(agg.runs=="R1",:)); sfR2=acc01_exact_signflip_family(agg.matrix(agg.runs=="R2",:));
A=double(paper.additive.a_erp_grid_uV(gridIndex+1)); beta=double(paper.additive.beta_gen_uV2_per_z(gridIndex)); expected=A^2*double(paper.bridge.expected_slope_gain_G_log10_per_uV2_per_z);
timing=struct("elapsed_s",toc(clock),"synthesis_s",w.timing.synthesis_s,"peak_memory_bytes",max([memMat,memProc]),"peak_mem_used_matlab_bytes",memMat,"peak_process_working_set_bytes",memProc,"memory_sampling","segment-end memory() MemUsedMATLAB plus PowerShell Get-Process PeakWorkingSet64; maxima retained");
checkpoint=struct("status","complete","analysis_id",paper.analysis_id,"scenario","additive","scenario_index",2,"grid_index",gridIndex,"grid_id",paper.additive.grid_id(gridIndex),"seed_index",seedIndex,"seed",seed,"A_uV",A,"beta_gen",beta,"expected_slope",expected,"truth_bridge_slope",expected,"parent_historical_contract_sha256",upper(string(sha256File(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")))),"input_hashes",inputHashManifest(paper.inputs),"subject_matrix",agg.table,"inference",struct("R1",sfR1,"R2",sfR2),"counts",struct("pairs",numel(w.pairs),"subjects",numel(unique(string({w.pairs.subject}))),"subject_rows",height(agg.table),"power_events",nPower,"phase_events",nPhase),"truth",w.truth,"started_at_utc",started,"completed_at_utc",utcNow(),"diagnostics",struct("robustfit_iteration_limit_count",diagCount,"robustfit_fit_count",fitCount),"timing",timing,"contract_sha256",upper(string(sha256File(fullfile(paperPath,"acc00_sim_powercurve_contract.json")))));
tmp=final+".tmp"; save(tmp,"checkpoint","-v7.3"); movefile(tmp,final,"f");
end

function assertIdentity(q,paper,root,gi,si,seed)
assert(string(q.status)=="complete"&&q.grid_index==gi&&q.seed_index==si&&q.seed==seed,"checkpoint identity mismatch");
assert(string(q.grid_id)==string(paper.additive.grid_id(gi))&&string(q.contract_sha256)==sha256File(fullfile(root,"scripts","ieeg","paper","acc00_sim_powercurve_contract.json")),"checkpoint grid/contract mismatch");
assert(string(q.parent_historical_contract_sha256)==sha256File(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")),"checkpoint parent contract mismatch");
end
function [mm,mp]=sampleMemory(mm,mp,pid)
try, u=memory; mm=max(mm,double(u.MemUsedMATLAB)); catch, end
try
    [status,out]=system(sprintf('powershell -NoProfile -Command "(Get-Process -Id %d).PeakWorkingSet64"',pid));
    if status==0, v=str2double(strtrim(out)); if isfinite(v), mp=max(mp,v); end, end
catch, end
end
function s=utcNow(),s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));end
function h=inputHashManifest(inputs),h=struct();for k=string(fieldnames(inputs))',h.(k)=upper(string(inputs.(k).sha256));end,end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance("SHA-256");while true,q=fread(fid,1048576,"*uint8");if isempty(q),break,end;md.update(typecast(q,"int8"));end;h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
