function run_p0_3_validator_driver(outputRoot, requireFull, poolSize)
%RUN_P0_3_VALIDATOR_DRIVER Detached-run wrapper for the P0-3 validator.
% Prints machine-readable start/finish lines for log-based monitoring.
arguments
    outputRoot (1,1) string
    requireFull (1,1) logical
    poolSize (1,1) double
end
root=string(fileparts(fileparts(fileparts(fileparts(mfilename("fullpath"))))));% <root>/scripts/ieeg/paper
cd(root);
fprintf("VALIDATOR_START outputRoot=%s requireFull=%d poolSize=%d pid=%d utc=%s\n",outputRoot,requireFull,poolSize,feature("getpid"),string(datetime("now","TimeZone","UTC")));
t0=tic;
try
    report=validate_acc00_p0_3_ablation(root,outputRoot,requireFull,poolSize);
    fprintf("VALIDATOR_DONE pass=%d fail=%d elapsed_s=%.1f\n",report.pass,report.fail,toc(t0));
catch ME
    fprintf("VALIDATOR_FAILED elapsed_s=%.1f message=%s\n",toc(t0),strrep(strrep(string(ME.message),newline," | "),sprintf("\r")," "));% machine-readable terminal state on failure; rethrow keeps nonzero exit
    rethrow(ME);
end
end
