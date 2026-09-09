function outputs = run_acc01_realdata_benchmark(varargin)
%RUN_ACC01_REALDATA_BENCHMARK Single-session ACC01 performance benchmark.

p = inputParser;
p.addParameter("ProjectRoot", "", @(x)ischar(x) || isstring(x));
p.addParameter("OutputDir", "", @(x)ischar(x) || isstring(x));
p.parse(varargin{:});
root = string(p.Results.ProjectRoot);
if strlength(root) == 0, root = string(fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))))); end
root = string(char(java.io.File(char(root)).getCanonicalPath()));
addpath(fullfile(root, "scripts", "ieeg", "acc01"));

contractPath = fullfile(root, "scripts", "ieeg", "acc01_contract.json");
contract = jsondecode(fileread(contractPath));
contractHash = sha256File(contractPath);
frozenHash = "D15F89687CCB1FEB4D0DFD9C902B12255876CA54551F8014229FD65E9FE32ABD";
assert(contractHash == frozenHash, "Bmovie:ACC01:ContractHash", "ACC01 contract SHA must remain frozen.");

entries = contract.inputs.raw_nwbs;
if iscell(entries), entries = [entries{:}]; end
match = string({entries.session_id})' == "P41CS_R1";
assert(nnz(match) == 1, "Bmovie:ACC01:BenchmarkSession", "Contract must contain exactly one P41CS_R1 raw-NWB entry.");
entry = entries(match);
assert(string(entry.subject) == "CS41" && string(entry.run) == "R1", ...
    "Bmovie:ACC01:BenchmarkSession", "P41CS_R1 manifest identity is not CS41/R1.");

outDir = string(p.Results.OutputDir);
if strlength(outDir) == 0
    outDir = fullfile(root, "processed", "subject", "group", "ieeg_acc01", "BangYoureDead", ...
        "checkpoints", "benchmark_cs41_p41cs_r1");
end
checkpointDir = fullfile(outDir, "sessions");
if ~isfolder(checkpointDir), mkdir(checkpointDir); end
checkpointPath = fullfile(checkpointDir, "session_P41CS_R1.mat");
assert(~isfile(checkpointPath), "Bmovie:ACC01:BenchmarkExists", "Refusing to overwrite benchmark checkpoint: %s", checkpointPath);

metadataInputs = requiredMetadataInputs(root);
covariates = readTsv(contract.inputs.label_covariates.path);
pid = feature("getpid");
snapshotStage = strings(0, 1);
snapshotMatlab = zeros(0, 1);
snapshotProcess = zeros(0, 1);

loadClock = tic;
[session, loadTiming] = acc01_load_real_session(root, contract, "P41CS_R1", "Mode", "full");
loadWall = toc(loadClock);
assert(string(session.status) == "complete", "Bmovie:ACC01:BenchmarkLoad", ...
    "Session load failed: %s", string(session.failure_reason));
session.metadata_inputs = metadataInputs;
assert(numel(session.analysis_pairs) == 14 && numel(unique(session.series_channel_indices)) == 16, ...
    "Bmovie:ACC01:BenchmarkShape", "Expected 14 pairs and 16 endpoints.");
observeStage("nwb_load");
observeStage("resample_filter");

[checkpoint, processTiming] = acc01_process_session(root, session, covariates, contract, contractHash, false, ...
    "CollectTiming", true, "StageObserver", @observeStage);
assert(string(checkpoint.status) == "complete", "Bmovie:ACC01:BenchmarkProcess", ...
    "Session processing returned a non-complete checkpoint.");

timing = struct("nwb_load_s", loadTiming.nwb_load_s, ...
    "resample_filter_s", loadTiming.resample_filter_s, ...
    "mpower_s", processTiming.mpower_s, "mphase_s", processTiming.mphase_s, ...
    "load_wall_s", loadWall, "support_s", processTiming.support_s);
checkpoint.benchmark = struct("timing", timing, "session_id", "P41CS_R1", ...
    "endpoint_count", numel(session.series_channel_indices), "pair_count", numel(session.analysis_pairs));

checkpointClock = tic;
acc01_write_session_checkpoint(checkpointPath, checkpoint);
checkpointSeconds = toc(checkpointClock);
observeStage("checkpoint_write");

stageNames = ["nwb_load"; "resample_filter"; "Mpower"; "Mphase"; "checkpoint_write"];
stageSeconds = [loadTiming.nwb_load_s; loadTiming.resample_filter_s; processTiming.mpower_s; ...
    processTiming.mphase_s; checkpointSeconds];
stageMem = [stageMemory("nwb_load"); stageMemory("resample_filter"); stageMemory("Mpower"); ...
    stageMemory("Mphase"); stageMemory("checkpoint_write")];
stagePeak = [stagePeakMemory("nwb_load"); stagePeakMemory("resample_filter"); stagePeakMemory("Mpower"); ...
    stagePeakMemory("Mphase"); stagePeakMemory("checkpoint_write")];
totalSeconds = sum(stageSeconds);
timestamp = repmat(utcNow(), numel(stageNames), 1);
benchmark = table(repmat("P41CS_R1", numel(stageNames), 1), stageNames, stageSeconds, stageMem, stagePeak, ...
    repmat(contractHash, numel(stageNames), 1), timestamp, ...
    'VariableNames', ["session_id" "stage" "seconds" "mem_matlab_bytes" "proc_peak_bytes" "contract_sha" "timestamp"]);
writetable(benchmark, fullfile(outDir, "benchmark.tsv"), "FileType", "text", "Delimiter", "\t");

manifest = struct("status", "complete", "session_id", "P41CS_R1", "subject", "CS41", "run", "R1", ...
    "contract_sha256", contractHash, "checkpoint", string(checkpointPath), ...
    "benchmark_tsv", string(fullfile(outDir, "benchmark.tsv")), "total_seconds", totalSeconds, ...
    "peak_memory_matlab_bytes", max(stageMem), "peak_working_set_bytes", max(stagePeak));
writeJson(fullfile(outDir, "benchmark_manifest.json"), manifest);
fprintf("ACC01 benchmark complete: %s (%.6f s total)\n", checkpointPath, totalSeconds);
fprintf("ACC01 benchmark peak memory: MATLAB %.0f bytes; process %.0f bytes\n", max(stageMem), max(stagePeak));
for i = 1:height(benchmark)
    fprintf("%s\t%.6f s\tMATLAB %.0f\tprocess %.0f\n", benchmark.stage(i), benchmark.seconds(i), ...
        benchmark.mem_matlab_bytes(i), benchmark.proc_peak_bytes(i));
end
outputs = struct("status", "complete", "checkpoint", string(checkpointPath), ...
    "benchmark_tsv", string(fullfile(outDir, "benchmark.tsv")), "contract_sha256", contractHash, ...
    "total_seconds", totalSeconds, "peak_memory_matlab_bytes", max(stageMem), ...
    "peak_working_set_bytes", max(stagePeak), "benchmark", benchmark);

    function observeStage(stage)
        [matlabBytes, processBytes] = memorySnapshot(pid);
        snapshotStage(end + 1, 1) = string(stage); %#ok<AGROW>
        snapshotMatlab(end + 1, 1) = matlabBytes; %#ok<AGROW>
        snapshotProcess(end + 1, 1) = processBytes; %#ok<AGROW>
    end
    function value = stageMemory(stage)
        x = snapshotMatlab(snapshotStage == stage);
        if isempty(x), value = NaN; else, value = max(x); end
    end
    function value = stagePeakMemory(stage)
        x = snapshotProcess(snapshotStage == stage);
        if isempty(x), value = NaN; else, value = max(x); end
    end
end

function metadata = requiredMetadataInputs(root)
i00Path = fullfile(root, "processed", "subject", "group", "ieeg_i00", "BangYoureDead", ...
    "task-bangyouredead_desc-i00-macro-local-car-groups.tsv");
anatomyPath = fullfile(root, "processed", "subject", "group", "ieeg_accg00_preflight", "BangYoureDead", ...
    "task-bangyouredead_desc-accg00preflight-acc-anatomy.tsv");
assert(isfile(i00Path) && isfile(anatomyPath), "Bmovie:ACC01:MissingMetadataInput", ...
    "Missing benchmark metadata input.");
metadata = struct("i00_local_car_groups_path", string(i00Path), ...
    "i00_local_car_groups_sha256", sha256File(i00Path), ...
    "accg00_acc_anatomy_path", string(anatomyPath), ...
    "accg00_acc_anatomy_sha256", sha256File(anatomyPath));
end

function [matlabBytes, processBytes] = memorySnapshot(pid)
matlabBytes = NaN; processBytes = NaN;
try
    m = memory;
    matlabBytes = double(m.MemUsedMATLAB);
catch
end
try
    [status, text] = system(sprintf('powershell -NoProfile -Command "(Get-Process -Id %d).PeakWorkingSet64"', pid));
    if status == 0, processBytes = str2double(strtrim(text)); end
catch
end
end

function out = utcNow()
out = string(datetime("now", "TimeZone", "UTC", "Format", "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end

function t = readTsv(path)
t = readtable(path, "FileType", "text", "Delimiter", "\t", "TextType", "string", ...
    "VariableNamingRule", "preserve");
end

function writeJson(path, value)
fid = fopen(path, "w"); assert(fid >= 0, "Bmovie:ACC01:BenchmarkWrite", "Cannot write %s", path);
cleanup = onCleanup(@() fclose(fid)); fprintf(fid, "%s\n", jsonencode(value, "PrettyPrint", true)); clear cleanup
end

function hash = sha256File(path)
fid = fopen(path, "rb"); assert(fid >= 0, "Bmovie:ACC01:BenchmarkHash", "Cannot open %s", path);
cleanup = onCleanup(@() fclose(fid)); md = java.security.MessageDigest.getInstance("SHA-256");
while true
    bytes = fread(fid, 1048576, "*uint8"); if isempty(bytes), break, end
    md.update(typecast(bytes, "int8"));
end
hash = upper(string(reshape(dec2hex(typecast(md.digest(), "uint8"), 2)', 1, [])));
end
