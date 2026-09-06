function outputs = run_acc00_nullpool(root, options)
%RUN_ACC00_NULLPOOL Authorized ACC00-SIM null seeds 1:400 only.
arguments
    root (1,1) string
    options.OutputDir (1,1) string = ""
    options.Workers (1,1) double {mustBeInteger,mustBePositive} = 12
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
addpath(fullfile(root,"scripts","ieeg","acc01"));
contractPath=fullfile(root,"scripts","ieeg","acc00_sim_contract.json");
contract=jsondecode(fileread(contractPath));contractHash=sha256File(contractPath);
assert(contract.world_plan.null_seeds==400,"Null pool requires exactly 400 seeds.");
workers=options.Workers;nWorlds=contract.world_plan.null_seeds;
seedIndex=(1:nWorlds)';seeds=contract.randomization.base_seed+100000+seedIndex;
outDir=options.OutputDir;
if strlength(outDir)==0,outDir=fullfile(root,"processed","subject","group","ieeg_acc00_sim","BangYoureDead","checkpoints","nullpool");end
if ~isfolder(outDir),mkdir(outDir);end
inputHashes=inputHashManifest(contract.inputs);reuseStatus=tryReuseBenchmark(root,outDir,contractHash,inputHashes,seeds(1));

todo=zeros(0,1);
for i=1:nWorlds
    final=worldPath(outDir,i);tmp=worldTmpPath(outDir,i);
    if isfile(final)
        assertValidCheckpoint(final,i,seeds(i),contractHash,inputHashes);
    else
        if isfile(tmp),delete(tmp);end
        todo(end+1,1)=i; %#ok<AGROW>
    end
end
fprintf("ACC00-SIM null pool: %d complete, %d missing; benchmark reuse: %s\n",nWorlds-numel(todo),numel(todo),reuseStatus);
pool=gcp("nocreate");
if isempty(pool)
    try
        pool=parpool("Processes",workers);
    catch firstPoolError
        c=parcluster('local');c.NumWorkers=workers;
        try
            pool=parpool(c,workers);
        catch secondPoolError
            throwAsCaller(addCause(secondPoolError,firstPoolError));
        end
    end
else
    assert(pool.NumWorkers==workers,"Existing pool has %d workers; expected %d.",pool.NumWorkers,workers);
end
fprintf("ACC00-SIM null pool using %d process workers.\n",pool.NumWorkers);baseSeed=contract.randomization.base_seed;stageClock=tic;
parfor k=1:numel(todo)
    i=todo(k);seed=baseSeed+100000+i;startedAt=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));worldClock=tic;
    world=acc00_sim_generate_world(root,contract,"benchmark",seed,"null",struct());
    result=acc01_run_metrics_optimized(world,contract);elapsed=toc(worldClock);
    subjectMatrix=result.subject_metrics.table;
    inference=reduceInference(result.signflip);
    pairRows=struct2table(result.pair_rows);
    surrogateSummary=pairRows(:,["subject" "run" "session_id" "pair_id" "n_phase_events" "domain_fraction" ...
        "surrogate_count" "refit_surrogate_count" "mphase_surrogate_reducer" "corr_series_diagnostic_count"]);
    timing=struct('elapsed_s',elapsed,'synthesis_s',world.timing.synthesis_s,'power_s',result.timing.power_s, ...
        'crossfit_s',result.timing.crossfit_s,'surrogate_ppc_s',result.timing.surrogate_ppc_s, ...
        'domain_s',result.timing.domain_s,'peak_memory_bytes',world.timing.peak_memory_bytes);
    diagnostics=result.diagnostics;completedAt=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
    checkpoint=struct('analysis_id',contract.analysis_id,'stage',"nullpool",'scenario',"null",'seed_index',i,'seed',seed, ...
        'status',"complete",'contract_sha256',contractHash,'input_hashes',inputHashes,'truth',world.truth, ...
        'started_at_utc',startedAt,'completed_at_utc',completedAt,'subject_matrix',subjectMatrix, ...
        'inference',inference,'surrogate_summary',surrogateSummary,'timing',timing,'diagnostics',diagnostics);
    tmp=worldTmpPath(outDir,i);final=worldPath(outDir,i);saveCheckpoint(tmp,checkpoint);movefile(tmp,final,"f");
    fprintf("ACC00-SIM null world %03d/%03d complete (seed %d, %.1f s, robustfit limits %d/%d)\n", ...
        i,nWorlds,seed,elapsed,diagnostics.robustfit_iteration_limit_count,diagnostics.robustfit_fit_count);
end
stageElapsed=toc(stageClock);
outputs=aggregateNullPool(outDir,nWorlds,seeds,contractHash,inputHashes,reuseStatus,workers,stageElapsed);
end

function status=tryReuseBenchmark(root,outDir,contractHash,inputHashes,expectedSeed)
target=worldPath(outDir,1);if isfile(target),status="existing_nullpool_checkpoint";return,end
source=fullfile(root,"processed","subject","group","ieeg_acc00_sim","BangYoureDead","checkpoints","benchmark","world_01.mat");
if ~isfile(source),status="not_reused_benchmark_missing";return,end
s=load(source,"checkpoint");
if ~isfield(s,"checkpoint"),status="not_reused_benchmark_schema";return,end
q=s.checkpoint;
if string(q.contract_sha256)~=contractHash,status="not_reused_contract_hash_mismatch";return,end
if q.seed~=expectedSeed||string(q.scenario)~="null",status="not_reused_seed_or_scenario_mismatch";return,end
if ~isfield(q,"input_hashes")||~isequaln(q.input_hashes,inputHashes),status="not_reused_input_hash_schema_or_mismatch";return,end
if ~isfield(q,"diagnostics")||~isfield(q.diagnostics,"robustfit_iteration_limit_count"),status="not_reused_missing_warning_diagnostic";return,end
status="not_reused_unrecognized_benchmark_schema";
end

function outputs=aggregateNullPool(outDir,nWorlds,seeds,contractHash,inputHashes,reuseStatus,workers,stageElapsed)
observedR1=nan(nWorlds,2);observedR2=nan(nWorlds,2);pR1=nan(nWorlds,2);pR2=nan(nWorlds,2);pMaxR1=nan(nWorlds,1);pMaxR2=nan(nWorlds,1);
elapsed=nan(nWorlds,1);synthesis=elapsed;power=elapsed;crossfit=elapsed;surrogatePpc=elapsed;domain=elapsed;peakMemory=elapsed;warningCount=zeros(nWorlds,1);fitCount=zeros(nWorlds,1);
started=strings(nWorlds,1);completed=strings(nWorlds,1);
for i=1:nWorlds
    s=load(worldPath(outDir,i),"checkpoint");q=s.checkpoint;assertValidCheckpointStruct(q,i,seeds(i),contractHash,inputHashes);
    observedR1(i,:)=q.inference.R1.observed;observedR2(i,:)=q.inference.R2.observed;pR1(i,:)=q.inference.R1.p_two_sided;pR2(i,:)=q.inference.R2.p_two_sided;
    pMaxR1(i)=q.inference.R1.p_maxstat;pMaxR2(i)=q.inference.R2.p_maxstat;elapsed(i)=q.timing.elapsed_s;synthesis(i)=q.timing.synthesis_s;power(i)=q.timing.power_s;
    crossfit(i)=q.timing.crossfit_s;surrogatePpc(i)=q.timing.surrogate_ppc_s;domain(i)=q.timing.domain_s;peakMemory(i)=q.timing.peak_memory_bytes;
    warningCount(i)=q.diagnostics.robustfit_iteration_limit_count;fitCount(i)=q.diagnostics.robustfit_fit_count;started(i)=q.started_at_utc;completed(i)=q.completed_at_utc;
end
worldTable=table((1:nWorlds)',seeds,started,completed,elapsed,synthesis,power,crossfit,surrogatePpc,domain,peakMemory,warningCount,fitCount,warningCount./fitCount, ...
    'VariableNames',["seed_index" "seed" "started_at_utc" "completed_at_utc" "elapsed_s" "synthesis_s" "power_s" "crossfit_s" ...
    "surrogate_ppc_s" "domain_s" "peak_memory_bytes" "robustfit_iteration_limit_count" "robustfit_fit_count" "robustfit_iteration_limit_fraction"]);
writetable(worldTable,fullfile(outDir,"nullpool_worlds.tsv"),"FileType","text","Delimiter","\t");

runs=["R1";"R1";"R2";"R2"];metrics=["Mpower";"Mphase";"Mpower";"Mphase"];values={observedR1(:,1);observedR1(:,2);observedR2(:,1);observedR2(:,2)};
n=repmat(nWorlds,4,1);meanValue=nan(4,1);stdValue=meanValue;medianValue=meanValue;q025=meanValue;q25=meanValue;q75=meanValue;q975=meanValue;
for i=1:4,v=values{i};meanValue(i)=mean(v);stdValue(i)=std(v,0);medianValue(i)=median(v);qq=prctile(v,[2.5 25 75 97.5],Method="inclusive");q025(i)=qq(1);q25(i)=qq(2);q75(i)=qq(3);q975(i)=qq(4);end
nullSummary=table(runs,metrics,n,meanValue,stdValue,medianValue,q025,q25,q75,q975, ...
    'VariableNames',["run" "metric" "n_worlds" "mean" "std_sigma_null" "median" "q2_5" "q25" "q75" "q97_5"]);
writetable(nullSummary,fullfile(outDir,"nullpool_distribution_summary.tsv"),"FileType","text","Delimiter","\t");

fprRun=["R1";"R1";"R1";"R2";"R2";"R2"];family=["Mpower";"Mphase";"maxstat_family";"Mpower";"Mphase";"maxstat_family"];
pValues={pR1(:,1);pR1(:,2);pMaxR1;pR2(:,1);pR2(:,2);pMaxR2};triggerCount=zeros(6,1);fpr=zeros(6,1);ciLow=zeros(6,1);ciHigh=zeros(6,1);
for i=1:6,triggerCount(i)=nnz(pValues{i}<=.05);[fpr(i),ci]=binofit(triggerCount(i),nWorlds,.05);ciLow(i)=ci(1);ciHigh(i)=ci(2);end
fprTable=table(fprRun,family,repmat(nWorlds,6,1),triggerCount,fpr,ciLow,ciHigh,repmat("Clopper-Pearson exact 95%",6,1), ...
    'VariableNames',["run" "test" "n_worlds" "trigger_count_alpha_0_05" "empirical_fpr" "ci95_low" "ci95_high" "ci_method"]);
writetable(fprTable,fullfile(outDir,"nullpool_fpr.tsv"),"FileType","text","Delimiter","\t");

qWarn=prctile(warningCount,[0 25 50 75 100],Method="inclusive");
diagnostics=table(nWorlds,sum(warningCount),mean(warningCount),qWarn(1),qWarn(2),qWarn(3),qWarn(4),qWarn(5),sum(fitCount),sum(warningCount)/sum(fitCount), ...
    'VariableNames',["n_worlds" "total_iteration_limit_count" "mean_per_world" "min_per_world" "q25_per_world" "median_per_world" "q75_per_world" "max_per_world" "total_robustfit_fits" "overall_fraction"]);
writetable(diagnostics,fullfile(outDir,"nullpool_robustfit_diagnostics.tsv"),"FileType","text","Delimiter","\t");

throughput=nWorlds/(stageElapsed/3600);conservativePeak=sum(maxk(peakMemory,min(workers,nWorlds)));
performance=table(nWorlds,workers,stageElapsed,throughput,max(peakMemory),conservativePeak,reuseStatus, ...
    'VariableNames',["n_worlds" "workers" "stage_elapsed_s" "worlds_per_hour" "max_single_world_peak_bytes" "sum_largest_worker_peaks_bytes" "benchmark_seed1_reuse_status"]);
writetable(performance,fullfile(outDir,"nullpool_performance.tsv"),"FileType","text","Delimiter","\t");
completion=struct('status',"complete",'n_worlds',nWorlds,'workers',workers,'contract_sha256',contractHash,'input_hashes',inputHashes, ...
    'completed_at_utc',string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")), ...
    'worlds_per_hour',throughput,'max_single_world_peak_bytes',max(peakMemory),'sum_12_largest_world_peaks_bytes',conservativePeak, ...
    'benchmark_seed1_reuse_status',reuseStatus);
writeJson(fullfile(outDir,"nullpool_complete.json"),completion);
outputs=struct('status',"complete",'directory',string(outDir),'distribution_summary',nullSummary,'fpr',fprTable,'diagnostics',diagnostics,'performance',performance,'contract_sha256',contractHash);
end

function inference=reduceInference(signflip)
inference=struct();
for run=["R1" "R2"]
    s=signflip.(run);inference.(run)=struct('subjects',s.subjects,'n_subjects',numel(s.subjects), ...
        'observed',s.observed,'p_two_sided',s.p_two_sided,'p_maxstat',s.p_maxstat, ...
        'studentization_sd',s.studentization_sd,'studentized_observed',s.studentized_observed, ...
        'max_observed',s.max_observed);
end
end

function hashes=inputHashManifest(inputs)
hashes=struct();roles=string(fieldnames(inputs));
for role=roles',hashes.(role)=upper(string(inputs.(role).sha256));end
end

function assertValidCheckpoint(path,index,seed,contractHash,inputHashes)
s=load(path,"checkpoint");assert(isfield(s,"checkpoint"),"Checkpoint schema missing: %s",path);assertValidCheckpointStruct(s.checkpoint,index,seed,contractHash,inputHashes);
end
function assertValidCheckpointStruct(q,index,seed,contractHash,inputHashes)
assert(string(q.status)=="complete"&&string(q.scenario)=="null"&&q.seed_index==index&&q.seed==seed,"Checkpoint identity mismatch at index %d.",index);
assert(string(q.contract_sha256)==contractHash,"Checkpoint contract hash mismatch at index %d.",index);
assert(isequaln(q.input_hashes,inputHashes),"Checkpoint input hash mismatch at index %d.",index);
end

function p=worldPath(outDir,index),p=fullfile(outDir,sprintf("world_%04d.mat",index));end
function p=worldTmpPath(outDir,index),p=fullfile(outDir,sprintf("world_%04d.tmp.mat",index));end
function saveCheckpoint(path,checkpoint),save(path,"checkpoint","-v7.3");end
function writeJson(path,value)
fid=fopen(path,"w");assert(fid>=0,"Cannot write %s",path);cleanup=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",jsonencode(value,"PrettyPrint",true));clear cleanup
end
function hash=sha256File(path)
fid=fopen(path,"rb");assert(fid>=0,"Cannot open %s",path);cleanup=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance("SHA-256");
while true,q=fread(fid,1048576,"*uint8");if isempty(q),break,end;md.update(typecast(q,"int8"));end
hash=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end
