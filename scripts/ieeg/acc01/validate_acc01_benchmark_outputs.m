function report = validate_acc01_benchmark_outputs(varargin)
%VALIDATE_ACC01_BENCHMARK_OUTPUTS Independent single-session benchmark validator.
% This validator only reads artifacts and never calls ACC01 metric functions.

p = inputParser;
p.addParameter("ProjectRoot", "", @(x)ischar(x) || isstring(x));
p.addParameter("OutputDir", "", @(x)ischar(x) || isstring(x));
p.parse(varargin{:});
root = string(p.Results.ProjectRoot);
if strlength(root) == 0, root = string(fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))))); end
root = string(char(java.io.File(char(root)).getCanonicalPath()));
outDir = string(p.Results.OutputDir);
if strlength(outDir) == 0
    outDir = fullfile(root, "processed", "subject", "group", "ieeg_acc01", "BangYoureDead", ...
        "checkpoints", "benchmark_cs41_p41cs_r1");
end

names = strings(0, 1); states = strings(0, 1); details = strings(0, 1);
contractPath = fullfile(root, "scripts", "ieeg", "acc01_contract.json");
contract = jsondecode(fileread(contractPath));
contractHash = sha256File(contractPath);
frozenHash = "F2C874719CA9F83010666538017570D4980E0E7D622A0BC0A1E5F4501368D2C9";
add("contract SHA", contractHash == frozenHash, contractHash);

manifestPath = fullfile(outDir, "benchmark_manifest.json");
tsvPath = fullfile(outDir, "benchmark.tsv");
checkpointPath = fullfile(outDir, "sessions", "session_P41CS_R1.mat");
add("benchmark artifacts", isfile(manifestPath) && isfile(tsvPath) && isfile(checkpointPath), ...
    "manifest, TSV and checkpoint");
if ~isfile(manifestPath) || ~isfile(tsvPath) || ~isfile(checkpointPath)
    checks = table(names, states, details, 'VariableNames', ["check" "status" "detail"]);
    report = struct("checks", checks, "summary", struct("pass", nnz(states == "PASS"), ...
        "fail", nnz(states == "FAIL"), "contract_sha256", contractHash));
    assert(report.summary.fail == 0, "Bmovie:ACC01:BenchmarkValidatorFailed", "Required benchmark artifact is missing.");
    return
end

manifest = jsondecode(fileread(manifestPath));
benchmark = readTsv(tsvPath);
loaded = load(checkpointPath, "checkpoint");
checkpoint = loaded.checkpoint;
loadedAgain = load(checkpointPath, "checkpoint");
add("manifest identity", string(manifest.session_id) == "P41CS_R1" && string(manifest.subject) == "CS41" && ...
    string(manifest.run) == "R1" && upper(string(manifest.contract_sha256)) == contractHash, "CS41/P41CS_R1/R1");

requiredColumns = ["session_id" "stage" "seconds" "mem_matlab_bytes" "proc_peak_bytes" "contract_sha" "timestamp"];
columnOK = all(ismember(requiredColumns, string(benchmark.Properties.VariableNames)));
requiredStages = ["nwb_load" "resample_filter" "Mpower" "Mphase" "checkpoint_write"];
stageOK = columnOK && all(arrayfun(@(s)nnz(string(benchmark.stage) == s) == 1, requiredStages));
positiveOK = columnOK && all(isfinite(double(benchmark.seconds))) && all(double(benchmark.seconds) > 0) && ...
    all(isfinite(double(benchmark.mem_matlab_bytes))) && all(double(benchmark.mem_matlab_bytes) > 0) && ...
    all(isfinite(double(benchmark.proc_peak_bytes))) && all(double(benchmark.proc_peak_bytes) > 0);
identityOK = columnOK && all(string(benchmark.session_id) == "P41CS_R1") && ...
    all(upper(string(benchmark.contract_sha)) == contractHash);
add("five stage timing rows", stageOK, strjoin(requiredStages, ","));
add("positive timing and memory", positiveOK, "seconds, MATLAB memory and process peak");
add("TSV identity", identityOK, "session and contract columns");

tmpFiles = dir(fullfile(outDir, "sessions", "*.tmp.mat"));
add("atomic checkpoint", isempty(tmpFiles), sprintf("%d temporary files", numel(tmpFiles)));

fields = ["analysis_id" "stage" "status" "contract_sha256" "subject" "session_id" "run" ...
    "input_hashes" "analysis_pairs" "pair_rows" "mphase_rng" "mphase_provenance" ...
    "bad_boundary_distance" "metadata_inventory" "synthetic_only" "source" "benchmark"];
schemaOK = all(isfield(checkpoint, fields)) && string(checkpoint.status) == "complete" && ...
    string(checkpoint.session_id) == "P41CS_R1" && string(checkpoint.subject) == "CS41" && ...
    string(checkpoint.run) == "R1" && upper(string(checkpoint.contract_sha256)) == contractHash;
add("checkpoint frozen schema", schemaOK, "complete identity and provenance fields");

pairCount = numel(checkpoint.analysis_pairs) == 14 && numel(checkpoint.pair_rows) == 14;
pairIds = string({checkpoint.analysis_pairs.pair_id});
pairCount = pairCount && numel(unique(pairIds)) == 14;
endpointCount = isfield(checkpoint.source, "series_channel_indices") && ...
    numel(unique(double(checkpoint.source.series_channel_indices))) == 16;
add("14 pairs and 16 endpoints", pairCount && endpointCount, ...
    sprintf("pairs=%d endpoints=%d", numel(checkpoint.analysis_pairs), numel(unique(double(checkpoint.source.series_channel_indices)))));

finiteSignals = false;
if ~isempty(checkpoint.analysis_pairs)
    finiteSignals = all(arrayfun(@(q)~isempty(q.signal) && all(isfinite(double(q.signal))), checkpoint.analysis_pairs));
end
finitePower = ~isempty(checkpoint.pair_rows) && all(isfinite([checkpoint.pair_rows.power]));
finiteCore = finiteSignals && finitePower && isequaln(checkpoint.analysis_pairs, loadedAgain.checkpoint.analysis_pairs);
add("checkpoint reread and finite core arrays", finiteCore, "bipolar signals and Mpower rows");

inputOK = isfield(checkpoint.input_hashes, "raw_nwb_expected") && ...
    upper(string(checkpoint.input_hashes.raw_nwb_expected)) == rawNwbHash(contract, "P41CS_R1");
add("session NWB provenance", inputOK, "raw_nwb_expected matches frozen manifest");

benchmarkFieldOK = isfield(checkpoint.benchmark, "timing") && ...
    double(checkpoint.benchmark.endpoint_count) == 16 && double(checkpoint.benchmark.pair_count) == 14;
add("checkpoint benchmark metadata", benchmarkFieldOK, "pair/endpoint counts and timing struct");

checks = table(names, states, details, 'VariableNames', ["check" "status" "detail"]);
summary = struct("pass", nnz(states == "PASS"), "fail", nnz(states == "FAIL"), "contract_sha256", contractHash, ...
    "session_id", "P41CS_R1", "pair_count", 14, "endpoint_count", 16);
writetable(checks, fullfile(outDir, "benchmark_validation.tsv"), "FileType", "text", "Delimiter", "\t");
writeJson(fullfile(outDir, "benchmark_validation.json"), summary);
report = struct("checks", checks, "summary", summary);
fprintf("ACC01 benchmark validator: %d pass / %d fail\n", summary.pass, summary.fail);
assert(summary.fail == 0, "Bmovie:ACC01:BenchmarkValidatorFailed", "Benchmark validator failed %d checks.", summary.fail);

    function add(name, ok, detail)
        names(end + 1, 1) = string(name); %#ok<AGROW>
        states(end + 1, 1) = ternary(ok, "PASS", "FAIL"); %#ok<AGROW>
        details(end + 1, 1) = string(detail); %#ok<AGROW>
    end
end

function hash = rawNwbHash(contract, sessionId)
entries = contract.inputs.raw_nwbs;
if iscell(entries), entries = [entries{:}]; end
match = string({entries.session_id})' == string(sessionId);
assert(nnz(match) == 1, "Bmovie:ACC01:BenchmarkManifest", "Raw-NWB manifest identity is not unique.");
hash = string(entries(match).sha256);
end

function out = ternary(condition, yes, no)
if condition, out = yes; else, out = no; end
end

function t = readTsv(path)
t = readtable(path, "FileType", "text", "Delimiter", "\t", "TextType", "string", ...
    "VariableNamingRule", "preserve");
end

function writeJson(path, value)
fid = fopen(path, "w"); assert(fid >= 0, "Bmovie:ACC01:BenchmarkValidatorWrite", "Cannot write output.");
cleanup = onCleanup(@() fclose(fid)); fprintf(fid, "%s\n", jsonencode(value, "PrettyPrint", true)); clear cleanup
end

function hash = sha256File(path)
fid = fopen(path, "rb"); assert(fid >= 0, "Bmovie:ACC01:BenchmarkValidatorHash", "Cannot open input.");
cleanup = onCleanup(@() fclose(fid)); md = java.security.MessageDigest.getInstance("SHA-256");
while true
    bytes = fread(fid, 1048576, "*uint8"); if isempty(bytes), break, end
    md.update(typecast(bytes, "int8"));
end
hash = upper(string(reshape(dec2hex(typecast(md.digest(), "uint8"), 2)', 1, [])));
end
