function export_acc00_p0_4_checkpoint_json(listPath, outputPath)
%EXPORT_ACC00_P0_4_CHECKPOINT_JSON MATLAB v7.3 checkpoint metadata bridge.
paths=readlines(listPath);records=cell(numel(paths),1);
for i=1:numel(paths)
    s=load(paths(i),"checkpoint");q=s.checkpoint;
    q.subject_metrics=table2struct(q.subject_metrics);
    records{i}=struct("filename",paths(i),"checkpoint",q);
end
fid=fopen(outputPath,"w");assert(fid>=0,"Cannot open JSON output.");cleanup=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",jsonencode(records));
end
