% Release port of the analysis-workspace tool of the same name; see README in this folder.
function manifest = stageb_main(root,seedIndices,outputRoot,workers)
%STAGEB_MAIN Formal Stage B orchestrator. Do not run before ZCode approval.
arguments
    root (1,1) string = string(pwd) % release port: run from the repository root
        % (the workspace original resolved the root from the tools/ layout)
    seedIndices (:,1) double {mustBeInteger,mustBePositive} = (1:60)'
    outputRoot (1,1) string = fullfile(fileparts(mfilename("fullpath")),"output")
    workers (1,1) double {mustBeInteger,mustBePositive} = 12
end
assert(isequal(seedIndices(:),(1:60)'),"Formal Stage B requires exactly seed_index 1:60.");
assert(workers==12,"Formal Stage B requires 12 workers.");if ~isfolder(outputRoot),mkdir(outputRoot);end
started=tic;manifest=struct("status","running","analysis_id","Bmovie-MPHASE-SUBJECT-REDUCER-POWER-STAGE-B-v1.0.0", ...
    "seed_indices",seedIndices,"matched_worlds",60,"arm_worlds",120,"workers",workers,"started_at_utc",utcNow());
writeJsonAtomic(fullfile(outputRoot,"run_manifest.json"),manifest);
pool=gcp("nocreate");if isempty(pool),pool=parpool("Processes",workers);end
assert(pool.NumWorkers==12,"Existing pool has %d workers; expected 12.",pool.NumWorkers);
parfor i=1:numel(seedIndices),stageb_run_world(root,seedIndices(i),"benchmark",outputRoot);end
stageb_export_bridge(root,outputRoot,"benchmark");toolDir=fileparts(mfilename("fullpath"));
[s1,t1]=system(sprintf('python "%s" "%s" --mode formal',fullfile(toolDir,"stageb_summarize.py"),outputRoot),"-echo");
assert(s1==0,"Stage B summary failed: %s",t1);
[s2,t2]=system(sprintf('python "%s" "%s" --mode formal',fullfile(toolDir,"stageb_validate.py"),outputRoot),"-echo");
assert(s2==0,"Stage B validator failed: %s",t2);
manifest.status="complete";manifest.completed_at_utc=utcNow();manifest.elapsed_s=toc(started);
writeJsonAtomic(fullfile(outputRoot,"run_manifest.json"),manifest);
end
function writeJsonAtomic(path,value),tmp=path+".tmp";fid=fopen(tmp,"w");assert(fid>=0);cl=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",jsonencode(value,"PrettyPrint",true));clear cl;movefile(tmp,path,"f");end
function s=utcNow(),s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));end
