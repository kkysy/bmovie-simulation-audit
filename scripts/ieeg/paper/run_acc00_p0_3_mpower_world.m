function checkpoint = run_acc00_p0_3_mpower_world(root, cellId, seedIndex, mode, outputRoot)
%RUN_ACC00_P0_3_MPOWER_WORLD One atomic/resumable synthetic Mpower checkpoint.
arguments
    root (1,1) string
    cellId (1,1) string
    seedIndex (1,1) double {mustBeInteger,mustBePositive}
    mode (1,1) string {mustBeMember(mode,["smoke" "synthetic"])} = "synthetic"
    outputRoot (1,1) string = ""
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));paperPath=fullfile(root,"scripts","ieeg","paper");addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(paperPath);
paper=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_3_ablation_contract.json")));base=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));cell=powerCell(paper,cellId);assert(seedIndex<=double(cell.worlds),"seed_index exceeds frozen cell size.");
if strlength(outputRoot)==0,outputRoot=fullfile(root,string(paper.outputs.root));end;outDir=fullfile(outputRoot,"mpower","checkpoints",cellId);if ~isfolder(outDir),mkdir(outDir),end
final=fullfile(outDir,sprintf("world_%03d.mat",seedIndex));seed=worldSeed(paper,seedIndex);
if isfile(final),s=load(final,"checkpoint");assertIdentity(s.checkpoint,paper,root,cellId,seedIndex,seed);checkpoint=s.checkpoint;return,end
started=utcNow();clock=tic;[world,inj]=acc00_sim_p0_3_generate_world(root,paper,base,cellId,seedIndex,mode);rows=repmat(struct("subject","","run","","power",NaN,"phase",0),numel(world.pairs),1);pairDiag=repmat(pairTemplate(),numel(world.pairs),1);fitCount=0;limitCount=0;nPower=0;overlapN=zeros(numel(world.pairs),1);overlapD=zeros(numel(world.pairs),1);overlapF=nan(numel(world.pairs),1);
for i=1:numel(world.pairs)
    [events,~,eventDiag]=acc01_build_support_and_events(world.pairs(i),world.covariates,base);[overlapN(i),overlapD(i),overlapF(i)]=predecessorOverlap(events,world.covariates,inj,cell);pw=acc01_compute_power_metric(world.pairs(i),events,base,"ReturnAperiodicDiagnostics",true);rows(i)=struct("subject",world.pairs(i).subject,"run",world.pairs(i).run,"power",pw.slope,"phase",0);pairDiag(i)=makePairDiag(i,events,pw,world,inj,overlapF(i));fitCount=fitCount+pw.robustfit_fit_count;limitCount=limitCount+pw.robustfit_iteration_limit_count;nPower=nPower+eventDiag.n_power;
end
if string(cell.kernel_position)=="pre",inj.predecessor_in_pre_overlap_n=NaN;inj.predecessor_in_pre_overlap_denominator=NaN;inj.predecessor_in_pre_overlap_fraction=NaN;inj.predecessor_in_pre_overlap_reason="no_post_component";else,inj.predecessor_in_pre_overlap_n=sum(overlapN);inj.predecessor_in_pre_overlap_denominator=sum(overlapD);inj.predecessor_in_pre_overlap_fraction=sum(overlapN)/max(sum(overlapD),1);inj.predecessor_in_pre_overlap_reason="pooled_analysis_power_events";end
agg=acc01_aggregate_subject_metrics(rows);sfR1=acc01_exact_signflip_family(agg.matrix(agg.runs=="R1",:));sfR2=acc01_exact_signflip_family(agg.matrix(agg.runs=="R2",:));
checkpoint=struct("status","complete","analysis_id",paper.analysis_id,"arm","Mpower","cell_id",cellId,"world_seed",seed,"seed_index",seedIndex,"scenario_index",double(paper.seed.P01_to_P10.scenario_index),"effect_index",double(paper.seed.P01_to_P10.effect_index), ...
    "contract_sha256",sha256File(fullfile(paperPath,"acc00_sim_p0_3_ablation_contract.json")),"parent_p0_2_contract_sha256",string(paper.parent_p0_2.sha256),"input_hashes",inputHashes(paper.inputs),"started_at_utc",started,"completed_at_utc",utcNow(), ...
    "timing",struct("elapsed_s",toc(clock),"synthesis_s",world.timing.synthesis_s,"peak_memory_bytes",world.timing.peak_memory_bytes),"diagnostics",struct("robustfit_iteration_limit_count",limitCount,"robustfit_fit_count",fitCount),"injection_provenance",inj,"pair_diagnostics",pairDiag,"subject_metrics",agg.table,"signflip",reduceInference(sfR1,sfR2), ...
    "world_summary",summarizeWorld(pairDiag,agg,nPower),"mode",mode);
tmp=final+".tmp.mat";save(tmp,"checkpoint","-v7.3");movefile(tmp,final,"f");
end
function q=pairTemplate(),q=struct("pair_row_index",NaN,"n_power_events",NaN,"offset_pre_mean",NaN,"offset_post_mean",NaN,"offset_delta_mean",NaN,"exponent_pre_mean",NaN,"exponent_post_mean",NaN,"exponent_delta_mean",NaN,"offset_pre_zbeta",NaN,"offset_post_zbeta",NaN,"offset_delta_zbeta",NaN,"exponent_pre_zbeta",NaN,"exponent_post_zbeta",NaN,"exponent_delta_zbeta",NaN,"raw_theta_pre_mean",NaN,"raw_theta_post_mean",NaN,"raw_theta_delta_mean",NaN,"residual_theta_delta_mean",NaN,"residual_theta_delta_zbeta",NaN,"injection_kernel_attempted_event_count",NaN,"injection_kernel_skipped_event_count",NaN,"predecessor_in_pre_overlap_fraction_post_component",NaN);end
function q=makePairDiag(i,events,pw,world,inj,overlapFraction)
a=pw.aperiodic;X=design(events);q=pairTemplate();q.pair_row_index=i;q.n_power_events=height(events);q.offset_pre_mean=mean(a.offset_pre_log10,"omitnan");q.offset_post_mean=mean(a.offset_post_log10,"omitnan");q.offset_delta_mean=mean(a.delta_offset_log10,"omitnan");q.exponent_pre_mean=mean(a.exponent_pre,"omitnan");q.exponent_post_mean=mean(a.exponent_post,"omitnan");q.exponent_delta_mean=mean(a.delta_exponent,"omitnan");q.offset_pre_zbeta=olsZ(X,a.offset_pre_log10);q.offset_post_zbeta=olsZ(X,a.offset_post_log10);q.offset_delta_zbeta=olsZ(X,a.delta_offset_log10);q.exponent_pre_zbeta=olsZ(X,a.exponent_pre);q.exponent_post_zbeta=olsZ(X,a.exponent_post);q.exponent_delta_zbeta=olsZ(X,a.delta_exponent);q.raw_theta_pre_mean=mean(a.raw_theta_pre_log10,"omitnan");q.raw_theta_post_mean=mean(a.raw_theta_post_log10,"omitnan");q.raw_theta_delta_mean=mean(a.raw_theta_delta_log10,"omitnan");q.residual_theta_delta_mean=mean(a.residual_theta_delta_log10,"omitnan");q.residual_theta_delta_zbeta=olsZ(X,a.residual_theta_delta_log10);q.injection_kernel_attempted_event_count=inj.per_pair_attempted_event_count(i);q.injection_kernel_skipped_event_count=inj.per_pair_skipped_out_of_bounds_count(i);
q.predecessor_in_pre_overlap_fraction_post_component=overlapFraction;
end
function X=design(events),X=[ones(height(events),1),double(events.z_e),double(events{:,endsWith(events.Properties.VariableNames,"_delta_z")|endsWith(events.Properties.VariableNames,"movie_time_linear_z")|endsWith(events.Properties.VariableNames,"movie_time_quadratic_z")})];end
function [n,d,f]=predecessorOverlap(events,covariates,inj,pCell)
if string(pCell.kernel_position)=="pre",n=NaN;d=NaN;f=NaN;return,end
d=height(events);n=0;injectionTimes=double(covariates.label_time_s(double(inj.sorted_selected_source_row_indices)));
for e=1:d,t=double(events.label_time_s(e));n=n+any(injectionTimes+.10<t-.10 & injectionTimes+.50>t-.50);end
f=n/max(d,1);
end
function out=reduceInference(r1,r2),out=struct("R1",reduceOne(r1),"R2",reduceOne(r2));end
function out=reduceOne(s),out=struct("observed",s.observed,"p_two_sided",s.p_two_sided,"p_maxstat",s.p_maxstat,"studentization_sd",s.studentization_sd,"studentized_observed",s.studentized_observed,"max_observed",s.max_observed);end
function b=olsZ(X,y),keep=all(isfinite([X y]),2);if nnz(keep)<=size(X,2),b=NaN;else,q=X(keep,:)\y(keep);b=q(2);end,end
function s=summarizeWorld(q,agg,nPower),v=struct2table(q);s=struct("n_pairs",numel(q),"n_power_events",nPower,"pair_equal_means",struct("offset_delta",mean(v.offset_delta_mean,"omitnan"),"exponent_delta",mean(v.exponent_delta_mean,"omitnan"),"raw_theta_delta",mean(v.raw_theta_delta_mean,"omitnan"),"residual_theta_delta",mean(v.residual_theta_delta_mean,"omitnan")),"subject_median_by_run",agg.table);end
function cell=powerCell(paper,id),z=paper.cells.Mpower;ix=find(string({z.id})==id,1);assert(~isempty(ix),"Unknown cell.");cell=z(ix);end
function s=worldSeed(paper,si),s=double(paper.seed.P01_to_P10.paper_base_seed)+100000*double(paper.seed.P01_to_P10.scenario_index)+1000*double(paper.seed.P01_to_P10.effect_index)+si;end
function assertIdentity(q,paper,root,id,si,seed),assert(string(q.status)=="complete"&&string(q.analysis_id)==string(paper.analysis_id)&&string(q.arm)=="Mpower"&&string(q.cell_id)==id&&q.seed_index==si&&q.world_seed==seed&&string(q.contract_sha256)==sha256File(fullfile(root,"scripts","ieeg","paper","acc00_sim_p0_3_ablation_contract.json")),"P0-3 Mpower checkpoint identity mismatch.");end
function h=inputHashes(x),h=struct();for k=string(fieldnames(x))',h.(k)=string(x.(k).sha256);end,end
function s=utcNow(),s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));b=fread(fid,Inf,"*uint8");md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(b,"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
