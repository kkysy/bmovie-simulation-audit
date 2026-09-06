function report = validate_acc00_p0_3_ablation(root, outputRoot, requireFull, poolSize)
%VALIDATE_ACC00_P0_3_ABLATION Independent contract/checkpoint assertions.
% poolSize=0 runs the serial client path. poolSize>0 runs per-checkpoint
% validations in a parfor pool of that size: each worker collects its own
% assertion rows, rows are concatenated in sorted file order, so serial and
% parallel runs emit byte-identical tables (validation plan B, frozen
% 2026-09-02). Same revision closes frozen-spec gaps 5.2/5.3/5.7/5.9
% (strengthening only; no production artifact is touched).
arguments
    root (1,1) string
    outputRoot (1,1) string = ""
    requireFull (1,1) logical = false
    poolSize (1,1) double = 0
end
maxNumCompThreads(1);% Recompute must reproduce production parfor-worker arithmetic (workers run single-threaded BLAS): multithreaded client recompute diverges up to 1.4e-5 in the ~3% of robustfit fits hitting the iteration limit; single-threaded recompute verified bitwise on the P01 profiling world (2026-09-02).
root=string(char(java.io.File(char(root)).getCanonicalPath()));paperPath=fullfile(root,"scripts","ieeg","paper");addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(paperPath);paper=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_3_ablation_contract.json")));base=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));if strlength(outputRoot)==0,outputRoot=fullfile(root,string(paper.outputs.root));end;outputRoot=string(char(java.io.File(char(outputRoot)).getCanonicalPath()));names=strings(0,1);status=names;detail=names;% canonicalized: manifest stores canonical-root paths; cosmetic outputRoot variants must not cause false FAIL (round-5 review)
contractSha=sha256File(fullfile(paperPath,"acc00_sim_p0_3_ablation_contract.json"));
add("contract schema, execution flags, 768 plan",string(paper.analysis_id)=="ACC00-SIM-P0-3-ABLATION-v1.0.0"&&~logical(paper.execution.implementation_authorized)&&~logical(paper.execution.smoke_authorized)&&logical(paper.execution.main_run_authorized)&&double(paper.execution.total_new_world_evaluations)==768&&numel(paper.cells.Mpower)==10&&numel(paper.cells.Mphase)==2,"frozen contract; main_run_authorized flipped true 2026-09-02 after user-adjudicated route b");
add("parent and input SHA256",sha256File(fullfile(root,string(paper.parent_p0_2.path)))==string(paper.parent_p0_2.sha256)&&string(paper.parent_acc00_null.edgeguard_contract_sha256)==string(paper.historical_references.Mphase_guard_1p5.contract_sha256)&&isfolder(fullfile(root,string(paper.historical_references.Mphase_guard_1p5.path)))&&allInputHashes(root,paper.inputs),"P0-2 parent, historical edgeguard reference, input hash manifest");
add("P10 deterministic FFT specification",runP10(root,outputRoot),"P10 formula/support/RMS/peak-bin");
mpowerFiles=sortByPath(dir(fullfile(outputRoot,"mpower","checkpoints","P*","world_*.mat")));mphaseFiles=sortByPath(dir(fullfile(outputRoot,"mphase","checkpoints","M*","world_*.mat")));
if poolSize>0
    pool=gcp('nocreate');if ~isempty(pool)&&pool.NumWorkers~=poolSize,delete(pool);pool=[];end
    if isempty(pool),pool=parpool('local',poolSize);end
    pctRunOnAll(sprintf('addpath(''%s'');addpath(''%s'');maxNumCompThreads(1);',char(fullfile(root,"scripts","ieeg","acc01")),char(paperPath)));% workers: pinned path + single-threaded BLAS so recompute reproduces production arithmetic
end
mpRes=runOverFiles(mpowerFiles,"Mpower",paper,base,root,contractSha,poolSize);
phRes=runOverFiles(mphaseFiles,"Mphase",paper,base,root,contractSha,poolSize);
threads=[cellfun(@(r)double(r.threads),mpRes);cellfun(@(r)double(r.threads),phRes)];
for k=1:numel(mpRes),names=[names;mpRes{k}.names];status=[status;mpRes{k}.status];detail=[detail;mpRes{k}.detail];end %#ok<AGROW>
for k=1:numel(phRes),names=[names;phRes{k}.names];status=[status;phRes{k}.status];detail=[detail;phRes{k}.detail];end %#ok<AGROW>
if requireFull,add("all 768 expected checkpoint identities",numel(mpowerFiles)==640&&numel(mphaseFiles)==128&&perCellCountsOK(outputRoot,paper)&&isempty(dir(fullfile(outputRoot,"**","*.tmp.mat"))),sprintf("mpower=%d mphase=%d; per-cell counts pinned to contract worlds",numel(mpowerFiles),numel(mphaseFiles)));else,add("no temporary checkpoint residue",isempty(dir(fullfile(outputRoot,"**","*.tmp.mat"))),"smoke subset permits incomplete cells");end
[okSeeds,seeds]=m01m02SemanticIdentity(outputRoot,paper);for j=1:numel(seeds),add("M01/M02 semantic field identity seed "+string(seeds(j)),okSeeds(j),"addendum 5.7: all semantic fields/inputs/seeds equal across M01/M02; only guard edge lower/upper differs; downstream numerics excluded by design");end
[okM,dM]=runManifestCrossCheck(outputRoot,paper,contractSha,mpowerFiles,mphaseFiles,requireFull);add("run_manifest checkpoint SHA256 cross-check",okM,dM);
if poolSize>0,threadOK=all(threads==1)&&pool.NumWorkers==poolSize;else,nT=maxNumCompThreads;threadOK=(nT==1);end
add("validator compute context single-threaded",threadOK,"maxNumCompThreads==1 in every compute context; pool size cross-checked in parallel mode and logged by driver");
checks=table(names,status,detail,'VariableNames',["check" "status" "detail"]);tableDir=fullfile(outputRoot,"tables");if ~isfolder(tableDir),mkdir(tableDir),end;writetable(checks,fullfile(tableDir,"p0_3_validator.tsv"),"FileType","text","Delimiter","\t");report=struct("checks",checks,"pass",nnz(status=="PASS"),"fail",nnz(status=="FAIL"));assert(report.fail==0,"P0-3 validator failed %d assertions.",report.fail);
    function add(n,ok,d),names(end+1,1)=string(n);status(end+1,1)=ternary(ok,"PASS","FAIL");detail(end+1,1)=string(d);end
end
function results=runOverFiles(files,arm,paper,base,root,contractSha,poolSize)
results=cell(numel(files),1);
if poolSize>0
    parfor k=1:numel(files)
        results{k}=validateOneCheckpoint(fullfile(files(k).folder,files(k).name),arm,paper,base,root,contractSha);
    end
else
    for k=1:numel(files)
        results{k}=validateOneCheckpoint(fullfile(files(k).folder,files(k).name),arm,paper,base,root,contractSha);
    end
end
end
function r=validateOneCheckpoint(path,arm,paper,base,root,contractSha)
if arm=="Mpower",r=validateOneMpower(path,paper,base,root,contractSha);else,r=validateOneMphase(path,paper,root,contractSha);end
end
function r=validateOneMpower(path,paper,base,root,contractSha)
names=strings(0,1);status=strings(0,1);detail=strings(0,1);
q=load(path,"checkpoint");c=q.checkpoint;idseed=string(c.cell_id)+"/"+c.seed_index;
ok=string(c.status)=="complete"&&string(c.arm)=="Mpower"&&string(c.contract_sha256)==contractSha&&filenameIdentity(path,c);
add("Mpower identity "+idseed,ok,"atomic checkpoint identity; filename cell/seed cross-checked against stored fields");
if ok
    expected=double(paper.seed.P01_to_P10.paper_base_seed)+100000*double(paper.seed.P01_to_P10.scenario_index)+1000*double(paper.seed.P01_to_P10.effect_index)+c.seed_index;add("Mpower seed/CRN "+idseed,c.world_seed==expected&&c.scenario_index==4&&c.effect_index==1,"cell id intentionally excluded from seed");
    inj=c.injection_provenance;source=readtable(fullfile(root,string(paper.inputs.label_covariates.path)),"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");if modeOf(c)=="smoke",source=source(1:12,:);source.label_time_s=linspace(2,13,12)';source.z_e=sin(2*pi*(0:11)'/3);source.source_complete(:)=true;end;cell=powerCell(paper,string(c.cell_id));[mask,sourceHash]=independentMask(source,cell,paper,c.world_seed);add("mask exact count/hash "+idseed,inj.N_source_complete==mask.N&&inj.K_retained==mask.K&&string(inj.mask_sha256)==sourceHash&&isequal(double(inj.sorted_selected_source_row_indices(:)),double(mask.rows(:))),"Threefry substream 2 mask");
    [w,winj]=acc00_sim_p0_3_generate_world(root,paper,base,string(c.cell_id),c.seed_index,modeOf(c));n=acc00_sim_generate_world(root,base,legacyMode(modeOf(c)),c.world_seed,"null",struct());
    add("analysis-event hash equals canonical "+idseed,analysisEventHash(w,base)==analysisEventHash(n,base),"addendum 5.3: analysis event selection invariant to injection mask/kernel");
    if string(c.cell_id)=="P01"
        w1=acc00_sim_generate_world(root,base,legacyMode(modeOf(c)),c.world_seed,"additive",struct("a_erp_uV",double(paper.additive.A_uV),"provisional",false));
        add("P01 canonical bitwise signal "+c.seed_index,worldEqual(w,w1),"frozen generator delegate");
        add("P01 canonical truth "+c.seed_index,truthEqual(w.truth,w1.truth),"addendum 5.2: scenario/effect/gain fields");
        add("P01 canonical support event table "+c.seed_index,supportTablesEqual(w,w1,base),"addendum 5.2: per-pair analysis support events");
        add("P01 canonical subject matrix/p/counts "+c.seed_index,canonicalInferenceEqual(w1,base,c),"addendum 5.2: canonical-world Mpower subject matrix, exact sign-flip p, pair counts vs stored");
    end
    if ismember(string(c.cell_id),["P04" "P05" "P06"]),add("pre kernel declaration "+c.seed_index,string(inj.kernel_position)=="pre"&&abs(inj.component_amplitude_uV-double(paper.additive.A_uV))<eps,"pre RMS/skip checked by wrapper");end
    if ismember(string(c.cell_id),["P07" "P08" "P09"]),add("symmetric energy declaration "+c.seed_index,abs(inj.component_amplitude_uV-double(paper.additive.A_uV)/sqrt(2))<eps,"A/sqrt(2) each component");end
    add("pair diagnostic schema "+idseed,numel(c.pair_diagnostics)==numel(c.subject_metrics.n_pairs)||~isempty(c.pair_diagnostics),"event-level arrays intentionally absent");
    add("Mpower CRN main-stream identity "+idseed,crnIdentityFromWorlds(w,n),"wrapper metadata equals canonical null world");
    add("Mpower deterministic event/pair recomputation "+idseed,recomputeMpower(c,paper,w,winj,base),"rebuild wrapper and recompute robustfit diagnostics");
end
r=struct("names",names,"status",status,"detail",detail,"threads",maxNumCompThreads);
    function add(nm,ok,d),names(end+1,1)=string(nm);status(end+1,1)=ternary(ok,"PASS","FAIL");detail(end+1,1)=string(d);end
end
function r=validateOneMphase(path,paper,root,contractSha)
names=strings(0,1);status=strings(0,1);detail=strings(0,1);
q=load(path,"checkpoint");c=q.checkpoint;idseed=string(c.cell_id)+"/"+c.seed_index;
ok=string(c.status)=="complete"&&string(c.arm)=="Mphase"&&string(c.contract_sha256)==contractSha&&filenameIdentity(path,c);
add("Mphase identity "+idseed,ok,"atomic replay identity; filename cell/seed cross-checked against stored fields");
if ok
    expected=double(paper.seed.M01_to_M02.base_seed)+100000*double(paper.seed.M01_to_M02.scenario_index)+c.seed_index;cell=phaseCell(paper,string(c.cell_id));add("Mphase replay seed and guard-only fields "+idseed,c.world_seed==expected&&c.effect_index==0&&c.guard.edge_guard_lower_s==cell.edge_guard_lower_s&&c.guard.edge_guard_upper_s==cell.edge_guard_upper_s&&c.guard.bad_dilation_s==.75,"expected replay pairing");
    for run=["R1" "R2"],if ~isfield(c.mphase_shape,char(run)),continue,end;s=c.mphase_shape.(char(run));x=double(s.subject_mphase_vector(:));if isempty(x),okr=s.orbit_n==0&&isnan(s.orbit_q50)&&isnan(s.orbit_sample_sd)&&isnan(s.observed_standardized_score)&&isnan(s.p_two_sided);add("Mphase empty-run shape "+idseed+"/"+run,okr,"smoke has no R2 subjects");continue,end;sf=acc01_exact_signflip_family([zeros(numel(x),1),x]);o=sf.statistics(:,2);okr=s.orbit_n==numel(o)&&abs(s.orbit_q50-prctile(o,50,"Method","inclusive"))<1e-14&&abs(s.orbit_sample_sd-std(o,0))<1e-14&&abs(s.observed_standardized_score-sf.studentized_observed(2))<1e-14&&s.p_two_sided==sf.p_two_sided(2);add("Mphase independent orbit shape "+idseed+"/"+run,okr,"exact signflip from persisted subject vector");end
    if modeOf(c)=="synthetic",add("Mphase historical replay match "+idseed,historicalReplayEqual(c,paper,root),"M02 bitwise; M01 <=1e-12 tolerance (route b)");else,add("Mphase historical replay match "+idseed,true,"not applicable to 4-pair smoke");end
end
r=struct("names",names,"status",status,"detail",detail,"threads",maxNumCompThreads);
    function add(nm,ok,d),names(end+1,1)=string(nm);status(end+1,1)=ternary(ok,"PASS","FAIL");detail(end+1,1)=string(d);end
end
function ok=historicalReplayEqual(c,paper,root)
% M02 bitwise: reference pool nullpool_edgeguard_20260829@E0930B71 and the replay both use
% the 632ca70 optimized phase path. M01 tolerance 1e-12 (user-adjudicated route b,
% 2026-09-02): its reference pool nullpool@6656A33A predates 632ca70 and used the
% reference phase path, whose equivalence to the optimized path was gated at <1e-12;
% observed max|diff| on the profiling M01 world = 8.7e-19 (R1) / 4.2e-15 (R2) on
% ~1e-4-scale Mphase values.
try
    if string(c.cell_id)=="M01",ref=paper.historical_references.Mphase_guard_0p75;tol=1e-12;else,ref=paper.historical_references.Mphase_guard_1p5;tol=0;end
    h=load(fullfile(root,string(ref.path),sprintf("world_%04d.mat",c.seed_index)),"checkpoint");ok=true;
    for run=["R1" "R2"],s=c.mphase_shape.(char(run));if isfield(h.checkpoint,"subject_matrix")&&istable(h.checkpoint.subject_matrix),t=h.checkpoint.subject_matrix;z=t(t.run==run,:);[sub,ix]=sort(string(z.subject));hv=double(z.Mphase(ix));rv=double(s.subject_mphase_vector(:));ok=ok&&isequal(string(s.subjects_lexicographic(:)),sub(:))&&isequal(size(hv),size(rv))&&closeEnough(rv,hv,tol);else,ok=ok&&closeEnough(s.observed,h.checkpoint.inference.(char(run)).observed(2),tol)&&closeEnough(s.p_two_sided,h.checkpoint.inference.(char(run)).p_two_sided(2),tol);end,end
catch,ok=false;end
end
function ok=closeEnough(a,b,tol),if tol==0,ok=isequaln(a,b);else,ok=all((abs(a-b)<=tol)|(isnan(a)&isnan(b)),"all");end,end
function ok=canonicalInferenceEqual(w1,base,c)
% Addendum 5.2: recompute the Mpower subject matrix and exact sign-flip p from
% the directly invoked canonical generator world and compare against the stored
% (wrapper-produced) checkpoint values.
try
    rows=repmat(struct("subject","","run","","power",NaN,"phase",0),numel(w1.pairs),1);
    for i=1:numel(w1.pairs)
        [events,~,~]=acc01_build_support_and_events(w1.pairs(i),w1.covariates,base);pw=acc01_compute_power_metric(w1.pairs(i),events,base);rows(i)=struct("subject",w1.pairs(i).subject,"run",w1.pairs(i).run,"power",pw.slope,"phase",0);
    end
    agg=acc01_aggregate_subject_metrics(rows);
    ok=isequal(string(agg.table.subject),string(c.subject_metrics.subject))&&isequal(string(agg.table.run),string(c.subject_metrics.run))&&isequaln(agg.table.Mpower,c.subject_metrics.Mpower)&&isequaln(agg.table.n_pairs,c.subject_metrics.n_pairs);
    for run=["R1" "R2"],sf=acc01_exact_signflip_family(agg.matrix(agg.runs==run,:));ok=ok&&signflipEqual(sf,c.signflip.(char(run)));end
catch,ok=false;end
end
function ok=signflipEqual(sf,s),ok=isequaln(sf.observed,s.observed)&&isequaln(sf.p_two_sided,s.p_two_sided)&&isequaln(sf.p_maxstat,s.p_maxstat)&&isequaln(sf.studentization_sd,s.studentization_sd)&&isequaln(sf.studentized_observed,s.studentized_observed)&&isequaln(sf.max_observed,s.max_observed);end
function ok=truthEqual(a,b),ok=string(a.scenario)==string(b.scenario)&&isequaln(a.a_erp_uV,b.a_erp_uV)&&isequaln(a.lambda_a_uV,b.lambda_a_uV)&&isequaln(a.kappa,b.kappa)&&logical(a.provisional)==logical(b.provisional)&&isequaln(a.gain_intercept_a,b.gain_intercept_a)&&isequaln(a.gain_slope_b,b.gain_slope_b);end
function ok=supportTablesEqual(w0,w1,base)
ok=numel(w0.pairs)==numel(w1.pairs);
for i=1:numel(w0.pairs)
    [e0,~,~]=acc01_build_support_and_events(w0.pairs(i),w0.covariates,base);[e1,~,~]=acc01_build_support_and_events(w1.pairs(i),w1.covariates,base);ok=ok&&tableEqual(e0,e1);
end
end
function ok=tableEqual(a,b)
ok=isequal(string(a.Properties.VariableNames),string(b.Properties.VariableNames))&&height(a)==height(b);
for v=string(a.Properties.VariableNames)
    x=a.(char(v));y=b.(char(v));
    if isnumeric(x)||islogical(x),ok=ok&&isequaln(double(x),double(y));else,ok=ok&&isequal(string(x),string(y));end
end
end
function h=analysisEventHash(world,base)
% SHA-256 over the per-pair analysis (strict-support power) event tables in pair
% order. NaN payloads are normalized before hashing: parser NaN bits are
% read-nondeterministic across readtable calls (see isequaln note below).
md=java.security.MessageDigest.getInstance("SHA-256");
for i=1:numel(world.pairs)
    [events,~,~]=acc01_build_support_and_events(world.pairs(i),world.covariates,base);md.update(typecast(typecast(uint32(height(events)),"uint8"),"int8"));
    for v=string(events.Properties.VariableNames)
        x=events.(char(v));
        md.update(typecast(unicode2native(char(v),"UTF-8"),"int8"));md.update(typecast(typecast(uint32(numel(x)),"uint8"),"int8"));% column name + element count framing keeps the byte stream injective across schemas (round-5 review)
        if isnumeric(x)||islogical(x)
            x=double(x(:));x(isnan(x))=NaN;b=typecast(x,"uint8");
        else
            s=string(x(:));s(ismissing(s))="<missing>";b=unicode2native(char(strjoin(s,char(0))),"UTF-8");% NUL delimiter: covariate strings cannot contain NUL
        end
        md.update(typecast(b,"int8"));
    end
end
h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end
function ok=crnIdentityFromWorlds(w,n)
try
    ok=isequaln(w.covariates,n.covariates)&&numel(w.pairs)==numel(n.pairs);% isequaln: real covariates contain NaN (parser NaN payloads are read-nondeterministic; isequal is false even for two reads of the same file)
    for i=1:numel(w.pairs),a=w.pairs(i);b=n.pairs(i);ok=ok&&isequaln([a.eta_subject a.eta_session a.eta_pair a.polarity a.phi_u_rad a.T],[b.eta_subject b.eta_session b.eta_pair b.polarity b.phi_u_rad b.T])&&isequaln(a.bad_segments_s,b.bad_segments_s);end
catch,ok=false;end
end
function ok=recomputeMpower(c,paper,w,inj,base)
try
    ok=numel(w.pairs)==numel(c.pair_diagnostics);rows=repmat(struct("subject","","run","","power",NaN,"phase",0),numel(w.pairs),1);
    for i=1:numel(w.pairs)
        [events,~,~]=acc01_build_support_and_events(w.pairs(i),w.covariates,base);pw=acc01_compute_power_metric(w.pairs(i),events,base,"ReturnAperiodicDiagnostics",true);rows(i)=struct("subject",w.pairs(i).subject,"run",w.pairs(i).run,"power",pw.slope,"phase",0);q=independentPairDiagnostic(i,events,pw,inj,powerCell(paper,string(c.cell_id)),w.covariates);ok=ok&&structNumericEqual(q,c.pair_diagnostics(i));
    end
    agg=acc01_aggregate_subject_metrics(rows);ok=ok&&isequaln(agg.table.Mpower,c.subject_metrics.Mpower)&&isequaln(agg.table.n_pairs,c.subject_metrics.n_pairs);
catch
    ok=false;
end
end
function q=independentPairDiagnostic(i,events,pw,inj,pCell,covariates)
a=pw.aperiodic;X=[ones(height(events),1),double(events.z_e),double(events{:,endsWith(events.Properties.VariableNames,"_delta_z")|endsWith(events.Properties.VariableNames,"movie_time_linear_z")|endsWith(events.Properties.VariableNames,"movie_time_quadratic_z")})];q=struct("pair_row_index",i,"n_power_events",height(events),"offset_pre_mean",mean(a.offset_pre_log10,"omitnan"),"offset_post_mean",mean(a.offset_post_log10,"omitnan"),"offset_delta_mean",mean(a.delta_offset_log10,"omitnan"),"exponent_pre_mean",mean(a.exponent_pre,"omitnan"),"exponent_post_mean",mean(a.exponent_post,"omitnan"),"exponent_delta_mean",mean(a.delta_exponent,"omitnan"),"offset_pre_zbeta",olsZ(X,a.offset_pre_log10),"offset_post_zbeta",olsZ(X,a.offset_post_log10),"offset_delta_zbeta",olsZ(X,a.delta_offset_log10),"exponent_pre_zbeta",olsZ(X,a.exponent_pre),"exponent_post_zbeta",olsZ(X,a.exponent_post),"exponent_delta_zbeta",olsZ(X,a.delta_exponent),"raw_theta_pre_mean",mean(a.raw_theta_pre_log10,"omitnan"),"raw_theta_post_mean",mean(a.raw_theta_post_log10,"omitnan"),"raw_theta_delta_mean",mean(a.raw_theta_delta_log10,"omitnan"),"residual_theta_delta_mean",mean(a.residual_theta_delta_log10,"omitnan"),"residual_theta_delta_zbeta",olsZ(X,a.residual_theta_delta_log10),"injection_kernel_attempted_event_count",inj.per_pair_attempted_event_count(i),"injection_kernel_skipped_event_count",inj.per_pair_skipped_out_of_bounds_count(i),"predecessor_in_pre_overlap_fraction_post_component",inj.predecessor_in_pre_overlap_fraction);
q.predecessor_in_pre_overlap_fraction_post_component=perPairOverlap(events,inj,pCell,covariates);
end
function f=perPairOverlap(events,inj,pCell,covariates)
if string(pCell.kernel_position)=="pre",f=NaN;return,end;times=double(covariates.label_time_s(double(inj.sorted_selected_source_row_indices)));n=0;for e=1:height(events),t=double(events.label_time_s(e));n=n+any(times+.10<t-.10 & times+.50>t-.50);end;f=n/max(height(events),1);
end
function b=olsZ(X,y),keep=all(isfinite([X y]),2);if nnz(keep)<=size(X,2),b=NaN;else,q=X(keep,:)\y(keep);b=q(2);end,end
function ok=structNumericEqual(a,b),fa=string(fieldnames(a));ok=isequal(fa,string(fieldnames(b)));for f=fa',ok=ok&&isequaln(double(a.(f)),double(b.(f)));end,end
function [okPerSeed,seeds]=m01m02SemanticIdentity(outputRoot,paper)
% Addendum 5.7: for every seed present in both M01 and M02, all semantic
% fields/inputs/seeds must be identical; only guard edge lower/upper differs
% (pinned to the contract cell values here). Downstream guard-dependent
% numerics (domain fractions, surrogate counts, Mphase values, shape) are
% excluded by design, as are timestamps/timing.
seeds=union(seedList(outputRoot,"M01"),seedList(outputRoot,"M02"));okPerSeed=true(numel(seeds),1);% union, not intersect: a seed missing on either side is a FAIL, never silently skipped (round-5 review)
m01c=phaseCell(paper,"M01");m02c=phaseCell(paper,"M02");
for j=1:numel(seeds)
    fa=fullfile(outputRoot,"mphase","checkpoints","M01",sprintf("world_%03d.mat",seeds(j)));fb=fullfile(outputRoot,"mphase","checkpoints","M02",sprintf("world_%03d.mat",seeds(j)));
    if ~isfile(fa)||~isfile(fb),okPerSeed(j)=false;continue,end
    qa=load(fa,"checkpoint");a=qa.checkpoint;
    qb=load(fb,"checkpoint");b=qb.checkpoint;
    okPerSeed(j)=string(a.status)==string(b.status)&&string(a.analysis_id)==string(b.analysis_id)&&string(a.arm)==string(b.arm)&&a.world_seed==b.world_seed&&a.seed_index==b.seed_index&&a.scenario_index==b.scenario_index&&a.effect_index==b.effect_index&&string(a.contract_sha256)==string(b.contract_sha256)&&string(a.parent_p0_2_contract_sha256)==string(b.parent_p0_2_contract_sha256)&&string(a.mode)==string(b.mode)&&a.guard.bad_dilation_s==b.guard.bad_dilation_s&&a.guard.edge_guard_lower_s==double(m01c.edge_guard_lower_s)&&a.guard.edge_guard_upper_s==double(m01c.edge_guard_upper_s)&&b.guard.edge_guard_lower_s==double(m02c.edge_guard_lower_s)&&b.guard.edge_guard_upper_s==double(m02c.edge_guard_upper_s)&&string(a.world_summary.fpr_role)==string(b.world_summary.fpr_role)&&string(a.world_summary.note)==string(b.world_summary.note)&&inputHashesEqual(a.input_hashes,b.input_hashes)&&pairIdentityEqual(a.pair_diagnostics,b.pair_diagnostics)&&subjectMetaEqual(a.subject_metrics,b.subject_metrics);
end
end
function s=seedList(outputRoot,cellId)
f=dir(fullfile(outputRoot,"mphase","checkpoints",cellId,"world_*.mat"));s=sort(str2double(erase(string({f.name}),["world_" ".mat"])));
end
function ok=inputHashesEqual(a,b),fa=string(fieldnames(a));ok=isequal(fa,string(fieldnames(b)));for f=fa',ok=ok&&string(a.(f))==string(b.(f));end,end
function ok=pairIdentityEqual(a,b)
ok=numel(a)==numel(b);
for f=["subject" "run" "session_id" "pair_id"],ok=ok&&isequal(string({a.(char(f))})',string({b.(char(f))})');end
end
function ok=subjectMetaEqual(a,b),ok=isequal(string(a.subject),string(b.subject))&&isequal(string(a.run),string(b.run))&&isequaln(a.Mpower,b.Mpower)&&isequaln(a.n_pairs,b.n_pairs);end
function [ok,d]=runManifestCrossCheck(outputRoot,paper,contractSha,mpowerFiles,mphaseFiles,requireFull)
% Addendum 5.9: cross-check run_manifest.json against the contract and the
% actual checkpoint tree (bijection + per-file SHA-256 recompute).
path=fullfile(outputRoot,"run_manifest.json");
if ~isfile(path)
    if requireFull,ok=false;d="run_manifest.json missing under requireFull";else,ok=true;d="not applicable: subset/smoke output root has no run_manifest.json";end
    return
end
m=jsondecode(fileread(path));ids=[string({paper.cells.Mpower.id})';string({paper.cells.Mphase.id})'];worlds=[double([paper.cells.Mpower.worlds])';double([paper.cells.Mphase.worlds])'];% two arms carry different struct fields; never struct-concat (jsondecode column struct pitfall, cf. orchestrator 72d9971)
countsOK=numel(m.per_cell_completion)==numel(ids);
for i=1:min(numel(m.per_cell_completion),numel(ids)),e=m.per_cell_completion(i);countsOK=countsOK&&string(e.cell_id)==ids(i)&&double(e.expected)==worlds(i)&&double(e.complete)==worlds(i);end
entries=m.checkpoint_sha256;mp=string({entries.filename})';mh=string({entries.sha256})';
actual=[mpowerFiles;mphaseFiles];ap=string(arrayfun(@(f)fullfile(f.folder,f.name),actual,"UniformOutput",false));
hashOK=true;for i=1:numel(mp),if ~isfile(mp(i))||sha256File(mp(i))~=mh(i),hashOK=false;break,end,end
[entryCell,entrySeed,parseOK]=parseCheckpointPaths(mp);for i=1:numel(mp),parseOK=parseOK&&string(entries(i).cell_id)==entryCell(i);end
perCellOK=perCellCountsOK(outputRoot,paper);
for i=1:numel(ids),perCellOK=perCellOK&&nnz(entryCell==ids(i))==worlds(i)&&isequal(sort(entrySeed(entryCell==ids(i)))',1:worlds(i));end
ok=string(m.status)=="complete"&&string(m.analysis_id)==string(paper.analysis_id)&&double(m.worlds)==768&&string(m.contract_sha256)==contractSha&&countsOK&&perCellOK&&parseOK&&numel(mp)==numel(ap)&&isequal(sort(mp),sort(ap))&&hashOK;
d=sprintf("addendum 5.9: manifest entries=%d files=%d; per-cell actual/manifest counts and seed sets pinned to contract; entry cell_id cross-checked against path",numel(mp),numel(ap));
end
function ok=perCellCountsOK(outputRoot,paper)
% Independent ground truth: actual per-cell directory counts must equal the
% contract worlds for every cell (arm totals alone cannot catch a redistributed
% or single-cell-depleted tree; round-5 review P0).
ok=perCellCountsOneArm(outputRoot,paper.cells.Mpower,"mpower")&&perCellCountsOneArm(outputRoot,paper.cells.Mphase,"mphase");
end
function ok=perCellCountsOneArm(outputRoot,cells,arm)
ok=true;for i=1:numel(cells),ok=ok&&numel(dir(fullfile(outputRoot,arm,"checkpoints",string(cells(i).id),"world_*.mat")))==double(cells(i).worlds);end
end
function ok=filenameIdentity(path,c)
% Filename cell/seed must equal the checkpoint's stored cell_id/seed_index, so a
% duplicated checkpoint under a different world_NNN name cannot pass identity.
[folder,fn,~]=fileparts(path);[~,cellDir]=fileparts(folder);tok=regexp(fn,"^world_(\d+)$","tokens","once");
ok=~isempty(tok)&&string(c.cell_id)==string(cellDir)&&c.seed_index==str2double(tok{1});
end
function [cells,seeds,ok]=parseCheckpointPaths(paths)
cells=strings(numel(paths),1);seeds=zeros(numel(paths),1);ok=true;
for i=1:numel(paths)
    [folder,fn,~]=fileparts(paths(i));[~,cellDir]=fileparts(folder);tok=regexp(fn,"^world_(\d+)$","tokens","once");
    if isempty(tok),ok=false;continue,end
    cells(i)=string(cellDir);seeds(i)=str2double(tok{1});
end
end
function files=sortByPath(files)
if isempty(files),return,end
p=string(arrayfun(@(f)fullfile(f.folder,f.name),files,"UniformOutput",false));[~,ix]=sort(p);files=files(ix);
end
function [m,h]=independentMask(t,cell,paper,seed),source=find(logical(t.source_complete));N=numel(source);r=double(paper.event_mask.levels.(char(cell.retention)));K=round(r*N);rs=RandStream("Threefry","Seed",seed);rs.Substream=double(paper.event_mask.rng.substream);take=sort(randperm(rs,N,K));v=false(height(t),1);v(source(take))=true;m=struct("N",N,"K",K,"rows",source(take));h=sha256Bytes(uint8(v));end
function ok=allInputHashes(root,inputs),ok=true;for k=string(fieldnames(inputs))',ok=ok&&sha256File(fullfile(root,string(inputs.(k).path)))==string(inputs.(k).sha256);end,end
function ok=runP10(root,out),try,r=check_acc00_p0_3_p10_fft(root,out);ok=r.status=="PASS";catch,ok=false;end,end
function c=powerCell(p,id),z=p.cells.Mpower;c=z(find(string({z.id})==id,1));end
function c=phaseCell(p,id),z=p.cells.Mphase;c=z(find(string({z.id})==id,1));end
function m=modeOf(c),m=string(c.mode);end
function out=legacyMode(m),if m=="smoke",out="smoke";else,out="benchmark";end,end
function ok=worldEqual(a,b),ok=numel(a.pairs)==numel(b.pairs)&&isequaln(a.covariates,b.covariates);if ~ok,return,end;for i=1:numel(a.pairs),if ~isequal(a.pairs(i).signal,b.pairs(i).signal)||a.pairs(i).background_scale_m~=b.pairs(i).background_scale_m||a.pairs(i).polarity~=b.pairs(i).polarity,ok=false;return,end,end,end
function s=ternary(c,a,b),if c,s=a;else,s=b;end,end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));h=sha256Bytes(fread(fid,Inf,"*uint8"));end
function h=sha256Bytes(b),md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(uint8(b(:)),"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
