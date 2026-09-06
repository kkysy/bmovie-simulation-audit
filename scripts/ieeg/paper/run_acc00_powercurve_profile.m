function out = run_acc00_powercurve_profile(root, gridIndex, seedIndex)
% Run exactly one authorized paper-only Mpower profiling world.
arguments
    root (1,1) string
    gridIndex (1,1) double {mustBeMember(gridIndex,[1 6 7])}
    seedIndex (1,1) double {mustBeInteger,mustBePositive} = 1
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
addpath(fullfile(root,"scripts","ieeg","acc01"));
contract=jsondecode(fileread(fullfile(root,"scripts","ieeg","paper","acc00_sim_powercurve_contract.json")));
assert(seedIndex<=double(contract.additive.grid_worlds(gridIndex)),"seedIndex exceeds frozen cell size");
seed=double(contract.seed.paper_base_seed)+100000*double(contract.seed.scenario_index)+1000*gridIndex+seedIndex;
outDir=fullfile(root,"processed","subject","group","ieeg_methods_paper","p0_2_design","profiling"); if ~isfolder(outDir),mkdir(outDir);end
final=fullfile(outDir,sprintf("world_g%02d_n%03d.mat",gridIndex,seedIndex)); tmp=final+".tmp";
if isfile(final), s=load(final,"checkpoint"); out=s.checkpoint; return; end
clock=tic; w=acc00_sim_generate_world(root,jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json"))),"benchmark",seed,"additive",struct("a_erp_uV",contract.additive.a_erp_grid_uV(gridIndex+1),"provisional",false));
rows=repmat(struct("subject","","run","","power",NaN,"phase",NaN),numel(w.pairs),1); diagCount=0; fitCount=0;
for p=1:numel(w.pairs), [ev,~,~]=acc01_build_support_and_events(w.pairs(p),w.covariates,jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")))); q=acc01_compute_power_metric(w.pairs(p),ev,jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")))); rows(p)=struct("subject",w.pairs(p).subject,"run",w.pairs(p).run,"power",q.slope,"phase",0); diagCount=diagCount+q.robustfit_iteration_limit_count; fitCount=fitCount+q.robustfit_fit_count; end
agg=acc01_aggregate_subject_metrics(rows); sfR1=acc01_exact_signflip_family(agg.matrix(agg.runs=="R1",:)); sfR2=acc01_exact_signflip_family(agg.matrix(agg.runs=="R2",:));
truthSlope=double(contract.additive.a_erp_grid_uV(gridIndex+1))^2*double(contract.bridge.expected_slope_gain_G_log10_per_uV2_per_z);
checkpoint=struct("status","complete","analysis_id",contract.analysis_id,"grid_index",gridIndex,"grid_id",contract.additive.grid_id(gridIndex),"seed_index",seedIndex,"seed",seed,"A_uV",contract.additive.a_erp_grid_uV(gridIndex+1),"truth_bridge_slope",truthSlope,"subject_matrix",agg.table,"inference",struct("R1",sfR1,"R2",sfR2),"diagnostics",struct("robustfit_iteration_limit_count",diagCount,"robustfit_fit_count",fitCount),"timing",struct("elapsed_s",toc(clock),"peak_memory_bytes",NaN),"contract_sha256",sha256File(fullfile(root,"scripts","ieeg","paper","acc00_sim_powercurve_contract.json")));
save(tmp,"checkpoint","-v7.3"); movefile(tmp,final,"f"); out=checkpoint;
end
function h=sha256File(path),fid=fopen(path,"rb");cl=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance("SHA-256");while true,q=fread(fid,1048576,"*uint8");if isempty(q),break,end;md.update(typecast(q,"int8"));end;h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
