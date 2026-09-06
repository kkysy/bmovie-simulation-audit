function report = run_acc00_p0_3_ablation_profile(root, workers)
%RUN_ACC00_P0_3_ABLATION_PROFILE User-authorized full-size profiling, 6 worlds.
% One representative world per distinct mechanism path (seed_index=1) written
% to the separate profiling/ tree; the canonical checkpoint tree is untouched
% so the pre-authorization contract SHA never collides with the main run.
arguments
    root (1,1) string
    workers (1,1) double {mustBeInteger,mustBePositive} = 12
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(fullfile(root,"scripts","ieeg","paper"));
out=fullfile(root,"processed","subject","group","ieeg_methods_paper","p0_3_ablation","profiling");
if ~isfolder(out),mkdir(out),end
jobs=[ "P01" "mpower"; "P02" "mpower"; "P04" "mpower"; "P07" "mpower"; "M01" "mphase"; "M02" "mphase"];
pool=gcp("nocreate");if isempty(pool),parpool(workers);elseif pool.NumWorkers~=workers,delete(pool);parpool(workers);end
clock=tic;elapsed=nan(size(jobs,1),1);peak=elapsed;
parfor k=1:size(jobs,1)
    if jobs(k,2)=="mpower"
        c=run_acc00_p0_3_mpower_world(root,jobs(k,1),1,"synthetic",out);
    else
        c=run_acc00_p0_3_mphase_world(root,jobs(k,1),1,"synthetic",out);
    end
    elapsed(k)=c.timing.elapsed_s;peak(k)=c.timing.peak_memory_bytes;
end
report=struct("cells",jobs(:,1),"seed_index",1,"elapsed_s",elapsed,"peak_memory_bytes",peak,"wall_s",toc(clock),"workers",workers);
fid=fopen(fullfile(out,"profile_report.json"),"w");fprintf(fid,"%s\n",jsonencode(report,"PrettyPrint",true));fclose(fid);
fprintf("P0-3 profiling complete: wall %.1f s\n",report.wall_s);
end
