function [checkpoint, timing] = acc01_process_session(root, session, covariates, contract, contractHash, syntheticOnly, varargin)
%ACC01_PROCESS_SESSION Shared per-session ACC01 checkpoint construction.
% Optional timing and observation are for the isolated benchmark wrapper only.

p = inputParser;
p.addParameter("CollectTiming", false, @(x)islogical(x) && isscalar(x));
p.addParameter("StageObserver", [], @(x)isempty(x) || isa(x, "function_handle"));
p.parse(varargin{:});
collectTiming = p.Results.CollectTiming;
observer = p.Results.StageObserver;
timing = emptyTiming();

started = utcNow();
if syntheticOnly
    checkpointHashes = inputHashes(contract.inputs);
    checkpointHashes.synthetic_world = "ACC00-SIM smoke fixture";
else
    checkpointHashes = sessionInputHashes(contract.inputs, session.session_id);
end
if string(session.status) ~= "complete"
    checkpoint = struct("analysis_id", contract.analysis_id, "stage", "realdata_session", ...
        "status", string(session.status), "contract_sha256", contractHash, ...
        "subject", string(session.subject), "session_id", string(session.session_id), "run", string(session.run), ...
        "started_at_utc", started, "completed_at_utc", utcNow(), "failure_reason", string(session.failure_reason), ...
        "input_hashes", checkpointHashes, ...
        "analysis_pairs", struct([]), "pair_rows", struct([]), "diagnostics", struct());
    return
end
if ~syntheticOnly && isfield(session, "metadata_inputs")
    checkpointHashes.i00_local_car_groups = session.metadata_inputs.i00_local_car_groups_sha256;
    checkpointHashes.accg00_acc_anatomy = session.metadata_inputs.accg00_acc_anatomy_sha256;
end
try
    [stream, rngInfo] = phaseStreamForSession(contractHash, session.session_id);
    rows = repmat(pairRowTemplate(), numel(session.analysis_pairs), 1);
    phaseProvenance = repmat(struct("pair_id", "", "refit_delta_s", [], "refit_surrogate_count", NaN), numel(session.analysis_pairs), 1);
    distanceRows = table();
    for i = 1:numel(session.analysis_pairs)
        pair = session.analysis_pairs(i);
        supportClock = tic;
        [powerEvents, phaseEvents, supportDiag] = acc01_build_support_and_events(pair, covariates, contract);
        supportSeconds = toc(supportClock);
        powerClock = tic;
        power = acc01_compute_power_metric(pair, powerEvents, contract);
        powerSeconds = toc(powerClock);
        if collectTiming
            timing.support_s = timing.support_s + supportSeconds;
            timing.mpower_s = timing.mpower_s + supportSeconds + powerSeconds;
        end
        observe(observer, "Mpower");

        phaseClock = tic;
        phase = acc01_compute_phase_metric(pair, phaseEvents, contract, stream);
        phaseSeconds = toc(phaseClock);
        if collectTiming
            timing.mphase_s = timing.mphase_s + phaseSeconds;
        end
        observe(observer, "Mphase");

        rows(i) = makePairRow(pair, power, phase, supportDiag, powerSeconds);
        phaseProvenance(i) = struct("pair_id", string(pair.pair_id), "refit_delta_s", phase.refit_delta_s, ...
            "refit_surrogate_count", phase.refit_surrogate_count);
        distanceRows = [distanceRows; badDistanceRows(pair, powerEvents)]; %#ok<AGROW>
    end
    metadata = metadataInventory(session, syntheticOnly);
    checkpoint = struct("analysis_id", contract.analysis_id, "stage", "realdata_session", ...
        "status", "complete", "contract_sha256", contractHash, "subject", string(session.subject), ...
        "session_id", string(session.session_id), "run", string(session.run), ...
        "started_at_utc", started, "completed_at_utc", utcNow(), ...
        "input_hashes", checkpointHashes, "source", stripSignals(session), ...
        "analysis_pairs", session.analysis_pairs, "pair_rows", rows, ...
        "mphase_rng", rngInfo, "mphase_provenance", phaseProvenance, ...
        "bad_boundary_distance", distanceRows, "metadata_inventory", metadata, ...
        "synthetic_only", syntheticOnly, "diagnostics", struct("n_pairs", numel(rows)));
catch ME
    checkpoint = struct("analysis_id", contract.analysis_id, "stage", "realdata_session", ...
        "status", "session_processing_failed", "contract_sha256", contractHash, ...
        "subject", string(session.subject), "session_id", string(session.session_id), "run", string(session.run), ...
        "started_at_utc", started, "completed_at_utc", utcNow(), "failure_reason", string(ME.identifier)+": "+string(ME.message), ...
        "input_hashes", checkpointHashes, "analysis_pairs", struct([]), "pair_rows", struct([]), "diagnostics", struct());
end
end

function timing = emptyTiming()
timing = struct("support_s", 0, "mpower_s", 0, "mphase_s", 0);
end

function observe(observer, stage)
if ~isempty(observer), observer(stage); end
end

function row = pairRowTemplate()
row = struct("subject", "", "run", "", "session_id", "", "pair_id", "", "power", NaN, "phase", NaN, ...
    "n_power_events", NaN, "n_phase_events", NaN, "domain_fraction", NaN, "surrogate_count", NaN, ...
    "refit_surrogate_count", NaN, "mphase_surrogate_reducer", "", "corr_series_diagnostic_count", NaN, ...
    "robustfit_iteration_limit_count", NaN, "robustfit_fit_count", NaN, "power_seconds", NaN, ...
    "crossfit_seconds", NaN, "surrogate_ppc_seconds", NaN, "domain_seconds", NaN);
end

function row = makePairRow(pair, power, phase, support, powerSeconds)
cheapCount = NaN;
if isfield(phase.diagnostics, "corr_series_cheap")
    cheapCount = phase.diagnostics.corr_series_cheap.count;
end
row = struct("subject", string(pair.subject), "run", string(pair.run), "session_id", string(pair.session_id), ...
    "pair_id", string(pair.pair_id), "power", power.slope, "phase", phase.difference, ...
    "n_power_events", support.n_power, "n_phase_events", support.n_phase, ...
    "domain_fraction", phase.surrogate.domain_fraction, "surrogate_count", phase.surrogate.sample_count, ...
    "refit_surrogate_count", phase.refit_surrogate_count, ...
    "mphase_surrogate_reducer", string(contractField(phase, "diagnostics", "reducer")), ...
    "corr_series_diagnostic_count", cheapCount, ...
    "robustfit_iteration_limit_count", power.robustfit_iteration_limit_count, "robustfit_fit_count", power.robustfit_fit_count, ...
    "power_seconds", powerSeconds, "crossfit_seconds", phase.timing.crossfit_s, ...
    "surrogate_ppc_seconds", phase.timing.surrogate_ppc_s, "domain_seconds", phase.timing.domain_s);
end

function value = contractField(s, first, second)
value = "";
if isfield(s, first) && isfield(s.(first), second), value = s.(first).(second); end
end

function rows = badDistanceRows(pair, events)
if height(events) == 0
    rows = table(); return
end
t = double(events.label_time_s(:)); d = inf(size(t));
for i = 1:size(pair.bad_segments_s, 1)
    b = double(pair.bad_segments_s(i, :));
    d = min(d, min(abs(t - b(1)), abs(t - b(2))));
end
rows = table(repmat(string(pair.subject), numel(t), 1), repmat(string(pair.session_id), numel(t), 1), ...
    repmat(string(pair.pair_id), numel(t), 1), double(events.label_index), t, d, ...
    'VariableNames', ["subject" "session_id" "pair_id" "label_index" "label_time_s" "nearest_bad_boundary_s"]);
end

function source = stripSignals(session)
source = session;
if isfield(source, "analysis_pairs")
    for i = 1:numel(source.analysis_pairs)
        source.analysis_pairs(i).signal = [];
    end
end
end

function metadata = metadataInventory(session, syntheticOnly)
ids = ["raw_PPC" "ERP_RMS" "uncorrected_ERP_time_frequency_power" "band_5_to_10_Hz"];
notComputed = table(ids', repmat("not_computed_pending_definition", numel(ids), 1), ...
    'VariableNames', ["diagnostic_id" "status"]);
pairTable = struct2table(session.analysis_pairs);
hemisphere = table(); subregion = table(); localCar = table();
if ~syntheticOnly && isfield(session, "pairs") && istable(session.pairs)
    hemisphere = groupsummary(session.pairs(:, ["hemisphere" "pair_name"]), "hemisphere");
    assert(isfield(session, "metadata_inputs"), "Bmovie:ACC01:MetadataInputs", "Missing preflight metadata input manifest.");
    car = readTsv(session.metadata_inputs.i00_local_car_groups_path);
    localCar = car(string(car.subject) == string(session.subject) & ...
        string(car.session_id) == string(session.session_id) & logical(double(car.qc_eligible)), :);
    anatomy = readTsv(session.metadata_inputs.accg00_acc_anatomy_path);
    indices = unique([double(session.pairs.channel_a_index); double(session.pairs.channel_b_index)]);
    subregion = anatomy(string(anatomy.subject) == string(session.subject) & ...
        string(anatomy.session_id) == string(session.session_id) & ...
        ismember(double(anatomy.series_channel_index), indices), ...
        ["subject" "session_id" "series_channel_index" "hemisphere" "region_family"]);
end
metadata = struct("not_computed", notComputed, "hemisphere", hemisphere, ...
    "frozen_subregion", subregion, "local_car_inventory", localCar, ...
    "n_analysis_pairs", height(pairTable));
end

function [stream, info] = phaseStreamForSession(contractHash, sessionId)
payload = char(string(contractHash) + "|" + string(sessionId));
bytes = unicode2native(payload, "UTF-8");
md = java.security.MessageDigest.getInstance("SHA-256");
md.update(typecast(uint8(bytes), "int8"));
digest = typecast(md.digest(), "uint8");
seed = uint32(digest(1)) + bitshift(uint32(digest(2)), 8) + ...
    bitshift(uint32(digest(3)), 16) + bitshift(uint32(digest(4)), 24);
stream = RandStream("Threefry", "Seed", double(seed));
info = struct("algorithm", "Threefry", "seed_uint32", double(seed), ...
    "seed_digest_sha256", upper(string(reshape(dec2hex(digest, 2)', 1, []))), ...
    "payload_utf8", string(payload));
end

function hashes = inputHashes(inputs)
hashes = struct();
for f = string(fieldnames(inputs))'
    if f == "raw_nwbs", continue, end
    if isfield(inputs.(f), "sha256"), hashes.(f) = string(inputs.(f).sha256); end
end
end

function hashes = sessionInputHashes(inputs, sessionId)
hashes = inputHashes(inputs);
entries = inputs.raw_nwbs;
if iscell(entries), entries = [entries{:}]; end
found = string({entries.session_id})' == string(sessionId);
assert(nnz(found) == 1, "Bmovie:ACC01:RawNwbManifest", "Session lacks one raw-NWB manifest entry.");
hashes.raw_nwb_expected = string(entries(found).sha256);
hashes.raw_nwb_relative_path = string(entries(found).relative_path);
end

function out = utcNow()
out = string(datetime("now", "TimeZone", "UTC", "Format", "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end

function t = readTsv(path)
t = readtable(path, "FileType", "text", "Delimiter", "\t", "TextType", "string", ...
    "VariableNamingRule", "preserve");
end
