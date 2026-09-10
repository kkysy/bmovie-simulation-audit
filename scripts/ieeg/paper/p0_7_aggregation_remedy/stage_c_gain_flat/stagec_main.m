% Release port of the analysis-workspace tool of the same name; see README in this folder.
function manifest = stagec_main(root,seedIndices,outputRoot,workers)
%STAGEC_MAIN Formal Stage C orchestrator. Do not run before ZCode approval.
arguments
    root (1,1) string = string(pwd) % release port: run from the repository root
        % (the workspace original resolved the root from the tools/ layout)
    seedIndices (:,1) double {mustBeInteger,mustBePositive} = (1:60)'
    outputRoot (1,1) string = fullfile(fileparts(mfilename("fullpath")),"output")
    workers (1,1) double {mustBeInteger,mustBePositive} = 12
end
assert(isequal(seedIndices(:),(1:60)'),"Formal Stage C requires exactly seed_index 1:60.");
assert(workers==12,"Formal Stage C requires 12 workers.");if ~isfolder(outputRoot),mkdir(outputRoot);end
started=tic;
manifest=struct("status","running","analysis_id","Bmovie-MPHASE-GAINFLAT-DIAGNOSTIC-STAGE-C-v1.0.0", ...
    "seed_indices",seedIndices,"worlds",60,"workers",workers,"started_at_utc",utcNow());
manifestPath=fullfile(outputRoot,"run_manifest.json");
% A resume rewrites this file; carry forward the first invocation's start time and
% the prior invocation log so the true first-run wall time is never overwritten.
if isfile(manifestPath)
    prev=jsondecode(fileread(manifestPath));
    if isfield(prev,"first_started_at_utc"),manifest.first_started_at_utc=string(prev.first_started_at_utc);end
    if isfield(prev,"invocations"),manifest.invocations=prev.invocations;end
end
if ~isfield(manifest,"first_started_at_utc"),manifest.first_started_at_utc=manifest.started_at_utc;end
writeJsonAtomic(manifestPath,manifest);
pool=gcp("nocreate");if isempty(pool),pool=parpool("Processes",workers);end
assert(pool.NumWorkers==12,"Existing pool has %d workers; expected 12.",pool.NumWorkers);
parfor i=1:numel(seedIndices),stagec_run_world(root,seedIndices(i),"benchmark",outputRoot);end
stagec_export_bridge(root,outputRoot,"benchmark");toolDir=fileparts(mfilename("fullpath"));
[s1,t1]=system(sprintf('python "%s" "%s" --mode formal',fullfile(toolDir,"stagec_summarize.py"),outputRoot),"-echo");
assert(s1==0,"Stage C summary failed: %s",t1);
[s2,t2]=system(sprintf('python "%s" "%s" --mode formal',fullfile(toolDir,"stagec_validate.py"),outputRoot),"-echo");
assert(s2==0,"Stage C validator failed: %s",t2);
manifest.status="complete";manifest.completed_at_utc=utcNow();manifest.elapsed_s=toc(started);
inv=struct("started_at_utc",manifest.started_at_utc,"completed_at_utc",manifest.completed_at_utc,"elapsed_s",manifest.elapsed_s);
if isfield(manifest,"invocations")
    prevInv=manifest.invocations;if isstruct(prevInv),prevInv=num2cell(prevInv(:));end
    manifest.invocations=[prevInv;{inv}];
else
    manifest.invocations={inv};
end
writeJsonAtomic(manifestPath,manifest);
end
function writeJsonAtomic(path,value),tmp=path+".tmp";fid=fopen(tmp,"w");assert(fid>=0);cl=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",jsonencode(value,"PrettyPrint",true));clear cl;movefile(tmp,path,"f");end
function s=utcNow(),s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));end
