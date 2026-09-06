function export_acc00_p0_3_checkpoint_json(listPath, outputPath)
%EXPORT_ACC00_P0_3_CHECKPOINT_JSON Explicit MATLAB table-to-JSON bridge.
paths=readlines(listPath);records=cell(numel(paths),1);
for i=1:numel(paths)
    s=load(paths(i),"checkpoint");c=s.checkpoint;
    r=struct("filename",paths(i),"checkpoint",c);
    r.checkpoint.subject_metrics=table2struct(c.subject_metrics);
    if isfield(c,"pair_diagnostics"),r.checkpoint.pair_diagnostics=c.pair_diagnostics;end
    if isfield(c,"world_summary")&&isfield(c.world_summary,"subject_median_by_run")&&istable(c.world_summary.subject_median_by_run),r.checkpoint.world_summary.subject_median_by_run=table2struct(c.world_summary.subject_median_by_run);end
    records{i}=r;
end
fid=fopen(outputPath,"w");assert(fid>=0,"Cannot open JSON output.");clean=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",jsonencode(records));
end
