% Release port of the analysis-workspace tool of the same name; see README in this folder.
function report = jrobust_smoke(root, workers)
%JROBUST_SMOKE Two full benchmark worlds; requires separate authorization.
arguments
    root (1,1) string = string(pwd) % release port: run from the repository root
        % (the workspace original resolved the root from the tools/ layout)
    workers (1,1) double {mustBeInteger,mustBePositive} = 2
end
root = string(char(java.io.File(char(root)).getCanonicalPath()));
toolDir = fileparts(mfilename("fullpath"));
c = jsondecode(fileread(fullfile(toolDir,"confirm_mphase_jrobust_contract.json")));
assert(logical(c.authorization.smoke_2_world),"Two-world smoke is not authorized by the child contract.");
assert(workers<=2,"Smoke uses at most two process workers.");
stamp = string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmssSSS'Z'"));
outputRoot = fullfile(toolDir,"output","smoke_"+stamp);
mkdir(outputRoot);
seedIndices = (1:2)';
pool = gcp("nocreate");
if isempty(pool), pool=parpool("Processes",workers); end
assert(pool.NumWorkers==workers,"Existing pool has %d workers; expected %d.",pool.NumWorkers,workers);
started = tic;
parfor i = 1:numel(seedIndices)
    jrobust_run_world(root,seedIndices(i),"smoke",outputRoot);
end
first = jrobust_run_world(root,1,"smoke",outputRoot);
assert(first.seed_index==1 && first.world_seed==20860828,"Smoke resume identity failed.");
tmpFiles = dir(fullfile(outputRoot,"checkpoints",string(c.world.arm_id),"*.tmp.mat"));
assert(isempty(tmpFiles),"Smoke left temporary checkpoints.");
bridge = jrobust_export_bridge(root,outputRoot,"smoke");
runPython(toolDir,"summarize_jrobust.py",outputRoot,"smoke","Smoke summary");
runPython(toolDir,"validate_jrobust.py",outputRoot,"smoke","Smoke validator");
decision = jsondecode(fileread(fullfile(outputRoot,"main_decision.json")));
validator = jsondecode(fileread(fullfile(outputRoot,"validator_report.json")));
worldSeconds = zeros(2,1);
checkpointBytes = zeros(2,1);
for i = 1:2
    path = fullfile(outputRoot,"checkpoints",string(c.world.arm_id),sprintf("world_%03d.mat",i));
    z=load(path,"checkpoint");
    worldSeconds(i)=double(z.checkpoint.timing.elapsed_s);
    info=dir(path);checkpointBytes(i)=double(info.bytes);
end
estimatedHours = median(worldSeconds)*120/(double(c.engineering.parallel_workers_formal)*3600);
report = struct("status","PASS","output",outputRoot,"worlds",2,"workers",workers, ...
    "elapsed_s",toc(started),"world_elapsed_s",worldSeconds,"checkpoint_bytes",checkpointBytes, ...
    "estimated_formal_hours",estimatedHours, ...
    "runtime_review_threshold_hours",double(c.engineering.smoke_runtime_review_threshold_hours), ...
    "requires_protocol_review",estimatedHours>double(c.engineering.smoke_runtime_review_threshold_hours), ...
    "bridge",bridge,"validator_status",string(validator.status), ...
    "decision_status",string(decision.status),"scientific_rules_evaluated",false);
writeJsonAtomic(fullfile(outputRoot,"smoke_report.json"),report);
disp(jsonencode(report,"PrettyPrint",true));
end

function runPython(toolDir,script,outputRoot,mode,label)
cmd=sprintf('python "%s" "%s" --mode %s',fullfile(toolDir,script),outputRoot,mode);
[status,text]=system(cmd,"-echo");
assert(status==0,"%s failed: %s",label,text);
end

function writeJsonAtomic(path,value)
tmp=path+".tmp";
fid=fopen(tmp,"w");assert(fid>=0,"Cannot write %s",tmp);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,"%s\n",jsonencode(value,"PrettyPrint",true));
clear cleanup
movefile(tmp,path,"f");
end
