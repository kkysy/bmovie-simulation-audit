function run_acc00_p0_4_historical_smoke_driver(root, outputPath)
%RUN_ACC00_P0_4_HISTORICAL_SMOKE_DRIVER Durable result for detached smoke.
try
    report=check_acc00_p0_4_mref_historical(root,1);
catch ME
    report=struct("status","FAIL","message",string(ME.message));
end
fid=fopen(outputPath,"w");assert(fid>=0,"Cannot write historical smoke result.");cleanup=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",jsonencode(report));
assert(string(report.status)=="PASS","P0-4 historical M-REF smoke failed: %s",string(report.message));
end
