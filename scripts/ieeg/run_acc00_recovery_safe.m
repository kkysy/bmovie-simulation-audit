function run_acc00_recovery_safe(root)
%RUN_ACC00_RECOVERY_SAFE Persist an unexpected recovery-run exception.
outDir=fullfile(root,"processed","subject","group","ieeg_acc00_sim", ...
    "BangYoureDead","checkpoints","recovery_mpower_only_20260829");
try
    run_acc00_recovery(root,Workers=12);
catch ME
    path=fullfile(outDir,"recovery_exception.txt");fid=fopen(path,"w");
    if fid>=0
        cleanup=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",getReport(ME,"extended","hyperlinks","off"));clear cleanup
    end
    rethrow(ME)
end
end
