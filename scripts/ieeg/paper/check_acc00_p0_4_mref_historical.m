function report = check_acc00_p0_4_mref_historical(root, seedIndex)
%CHECK_ACC00_P0_4_MREF_HISTORICAL Direct P0-4 M-REF route vs M02 history.
arguments
    root (1,1) string
    seedIndex (1,1) double {mustBeInteger,mustBePositive} = 1
end
maxNumCompThreads(1);% M02 historical route is bitwise only in the frozen single-threaded context.
root=string(char(java.io.File(char(root)).getCanonicalPath()));addpath(fullfile(root,"scripts","ieeg"));addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(fullfile(root,"scripts","ieeg","paper"));p=jsondecode(fileread(fullfile(root,"scripts","ieeg","paper","acc00_sim_p0_4_panel_contract.json")));base=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));ref=p.historical_references.Mphase_guard_1p5;seed=20260827+100000+seedIndex;world=acc00_sim_generate_world(root,base,"benchmark",seed,"null",struct());sim=base;sim.surrogate.domain_guard_edge_lower_s=1.5;sim.surrogate.domain_guard_edge_upper_s=1.5;sim.surrogate.domain_guard_bad_dilation_s=.75;rows=repmat(struct("subject","","run","","power",0,"phase",NaN),numel(world.pairs),1);
for i=1:numel(world.pairs),[~,events]=acc01_build_support_and_events(world.pairs(i),world.covariates,sim);q=acc01_compute_phase_metric_optimized(world.pairs(i),events,sim,world.stream);rows(i)=struct("subject",world.pairs(i).subject,"run",world.pairs(i).run,"power",0,"phase",q.difference);end
a=acc01_aggregate_subject_metrics(rows);h=load(fullfile(root,string(ref.path),sprintf("world_%04d.mat",seedIndex)),"checkpoint");ok=true;maxDiff=0;
for r=["R1" "R2"],t=h.checkpoint.subject_matrix;z=t(t.run==r,:);[subjects,ix]=sort(string(z.subject));got=a.table(a.table.run==r,:);[gs,gix]=sort(string(got.subject));d=max(abs(double(got.Mphase(gix))-double(z.Mphase(ix))),[],"omitnan");if isempty(d),d=0;end;maxDiff=max(maxDiff,d);ok=ok&&isequal(subjects,gs)&&isequaln(double(got.Mphase(gix)),double(z.Mphase(ix)));end
report=struct("status",ternary(ok,"PASS","FAIL"),"seed_index",seedIndex,"world_seed",seed,"route","M02_bitwise","max_abs_difference",maxDiff,"historical_contract_sha256",string(ref.contract_sha256));assert(ok,"P0-4 M-REF historical M02 route did not reproduce bitwise.");
end
function x=ternary(c,a,b),if c,x=a;else,x=b;end,end
