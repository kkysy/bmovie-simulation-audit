% Release port of the analysis-workspace tool of the same name; see README in this folder.
function manifest = jrobust_main(root, seedIndices, outputRoot, workers)
%JROBUST_MAIN Formal 120-world orchestrator. Requires contract authorization.
arguments
    root (1,1) string = string(pwd) % release port: run from the repository root
        % (the workspace original resolved the root from the tools/ layout)
    seedIndices (:,1) double {mustBeInteger,mustBePositive} = (1:120)'
    outputRoot (1,1) string = fullfile(fileparts(mfilename("fullpath")),"output")
    workers (1,1) double {mustBeInteger,mustBePositive} = 8
end
toolDir = fileparts(mfilename("fullpath"));
c = jsondecode(fileread(fullfile(toolDir,"confirm_mphase_jrobust_contract.json")));
assert(logical(c.authorization.formal_120_world),"Formal 120-world run is not authorized by the child contract.");
assert(isequal(seedIndices(:),(1:120)'),"Formal J-robust run requires exactly seed_index 1:120.");
assert(workers==double(c.engineering.parallel_workers_formal),"Frozen formal worker count is %d.",double(c.engineering.parallel_workers_formal));
if ~isfolder(outputRoot), mkdir(outputRoot); end
started = tic;
manifestPath = fullfile(outputRoot,"run_manifest.json");
manifest = struct("status","running","analysis_id",string(c.analysis_id),"run_kind","formal", ...
    "seed_indices",seedIndices,"workers",workers,"started_at_utc",utcNow());
if isfile(manifestPath)
    previous = jsondecode(fileread(manifestPath));
    if isfield(previous,"first_started_at_utc"), manifest.first_started_at_utc=string(previous.first_started_at_utc); end
    if isfield(previous,"invocations"), manifest.invocations=previous.invocations; end
end
if ~isfield(manifest,"first_started_at_utc"), manifest.first_started_at_utc=manifest.started_at_utc; end
writeJsonAtomic(manifestPath,manifest);
pool = gcp("nocreate");
if isempty(pool), pool=parpool("Processes",workers); end
assert(pool.NumWorkers==workers,"Existing pool has %d workers; expected %d.",pool.NumWorkers,workers);
parfor i = 1:numel(seedIndices)
    jrobust_run_world(root,seedIndices(i),"formal",outputRoot);
end
jrobust_export_bridge(root,outputRoot,"formal");
runPython(toolDir,"summarize_jrobust.py",outputRoot,"formal","Main summary");
runPython(toolDir,"validate_jrobust.py",outputRoot,"formal","Independent validator");
manifest.status = "complete";
manifest.completed_at_utc = utcNow();
manifest.elapsed_s = toc(started);
invocation = struct("started_at_utc",manifest.started_at_utc, ...
    "completed_at_utc",manifest.completed_at_utc,"elapsed_s",manifest.elapsed_s);
if isfield(manifest,"invocations")
    prior = manifest.invocations;
    if isstruct(prior), prior=num2cell(prior(:)); end
    manifest.invocations = [prior;{invocation}];
else
    manifest.invocations = {invocation};
end
writeJsonAtomic(manifestPath,manifest);
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

function s=utcNow()
s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end
