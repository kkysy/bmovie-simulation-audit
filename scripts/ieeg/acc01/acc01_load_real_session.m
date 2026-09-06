function [out, stageTiming] = acc01_load_real_session(root, contract, sessionId, varargin)
%ACC01_LOAD_REAL_SESSION Assemble one real ACC01 session without metric math.
% The full mode implements only the frozen I00 -> endpoint -> 250-Hz bipolar
% data path.  Connectivity mode is deliberately limited to one channel and a
% short HDF5 slice; it does not return a signal or calculate a metric.

p = inputParser;
p.addParameter("Mode", "full", @(x)ischar(x) || isstring(x));
p.addParameter("ConnectivitySamples", 16, @(x)isnumeric(x) && isscalar(x) && x >= 1);
p.parse(varargin{:});
mode = lower(string(p.Results.Mode));
assert(ismember(mode, ["full" "connectivity_check"]), ...
    "Bmovie:ACC01:LoaderMode", "Mode must be full or connectivity_check.");
stageTiming = struct("nwb_load_s", NaN, "resample_filter_s", NaN, ...
    "raw_sample_count", NaN, "endpoint_count", NaN);

root = string(char(java.io.File(char(root)).getCanonicalPath()));
sessionId = string(sessionId);

try
    pairs = readTsv(contract.inputs.i00_acc_pairs.path);
    channels = readTsv(contract.inputs.i00_channel_qc.path);
    bad = readTsv(contract.inputs.i00_bad_segments.path);
    timing = readTsv(contract.inputs.i00_timing.path);
    eligible = logical(double(pairs.qc_eligible)) & contains(string(pairs.location), "ACC");
    pairs = pairs(eligible, :);
    sessionPairs = pairs(string(pairs.session_id) == sessionId, :);
    timingRow = timing(string(timing.session_id) == sessionId, :);
    if height(sessionPairs) == 0
        out = failed(sessionId, mode, "input_integrity_failed", "No eligible ACC pair for session.");
        return
    end
    if height(timingRow) ~= 1
        out = failed(sessionId, mode, "input_integrity_failed", "Timing join is not unique.");
        return
    end
    if ~all(string(sessionPairs.subject) == string(timingRow.subject(1))) || ...
            ~all(string(sessionPairs.run) == string(timingRow.run(1)))
        out = failed(sessionId, mode, "input_integrity_failed", "Pair/timing subject or run mismatch.");
        return
    end
    if double(timingRow.macro_rate_hz(1)) ~= 1000
        out = failed(sessionId, mode, "input_integrity_failed", "Frozen macro rate is not 1000 Hz.");
        return
    end

    first = round(double(timingRow.macro_movie_first_sample(1)));
    last = round(double(timingRow.macro_movie_last_sample(1)));
    n = last - first + 1;
    if first < 1 || last < first || n < 1
        out = failed(sessionId, mode, "input_integrity_failed", "Invalid movie sample bounds.");
        return
    end
    nwbFile = fullfile(root, string(timingRow.nwb_file(1)));
    if ~isfile(nwbFile)
        out = failed(sessionId, mode, "input_integrity_failed", "NWB file is missing.");
        return
    end

    endpoint = unique([double(sessionPairs.channel_a_index); double(sessionPairs.channel_b_index)], "stable");
    sessionChannels = channels(channels.series == "LFP_macro" & ...
        string(channels.session_id) == sessionId, :);
    found = ismember(endpoint, double(sessionChannels.series_channel_index));
    if ~all(found)
        out = failed(sessionId, mode, "input_integrity_failed", "Pair endpoint lacks channel-QC series index.");
        return
    end

    hdf5Path = string(contract.raw_lfp_binding.hdf5_path);
    if mode == "connectivity_check"
        count = min(double(p.Results.ConnectivitySamples), n);
        trace = double(h5read(nwbFile, hdf5Path, [endpoint(1), first], [1, count]));
        out = struct("status", "complete_io_connectivity_only", "mode", mode, ...
            "failure_reason", "", "subject", string(timingRow.subject(1)), ...
            "session_id", sessionId, "run", string(timingRow.run(1)), ...
            "nwb_file", string(nwbFile), "hdf5_path", hdf5Path, ...
            "orientation", string(contract.raw_lfp_binding.orientation), ...
            "series_channel_index", endpoint(1), "first_sample", first, ...
            "sample_count", count, "sampling_rate_hz", double(timingRow.macro_rate_hz(1)), ...
            "all_finite", all(isfinite(trace(:))), "pairs", struct([]), ...
            "analysis_pairs", struct([]), "diagnostics", struct());
        if ~out.all_finite
            out.status = "input_integrity_failed";
            out.failure_reason = "Connectivity slice contains non-finite values.";
        end
        return
    end

    stageTiming.endpoint_count = numel(endpoint);
    stageTiming.raw_sample_count = n;
    loadClock = tic;
    raw = nan(numel(endpoint), n);
    for i = 1:numel(endpoint)
        raw(i, :) = double(h5read(nwbFile, hdf5Path, [endpoint(i), first], [1, n]));
    end
    stageTiming.nwb_load_s = toc(loadClock);
    if any(~isfinite(raw), "all")
        out = failed(sessionId, mode, "input_integrity_failed", "Complete endpoint read contains non-finite values.");
        out.subject = string(timingRow.subject(1)); out.run = string(timingRow.run(1));
        out.nwb_file = string(nwbFile); out.series_channel_indices = endpoint;
        return
    end

    filterClock = tic;
    [b, a] = butter(3, 100 / 500);
    endpoint250 = nan(numel(endpoint), numel(1:4:n));
    for i = 1:numel(endpoint)
        filtered = filtfilt(b, a, raw(i, :)');
        endpoint250(i, :) = filtered(1:4:end)';
    end
    stageTiming.resample_filter_s = toc(filterClock);
    clear raw

    movieZeroCommon = double(timingRow.macro_common_start_s(1)) + ...
        (double(timingRow.macro_movie_first_sample(1)) - 1) / double(timingRow.macro_rate_hz(1));
    T = double(timingRow.usable_movie_stop_s(1));
    template = struct("subject", "", "session_id", "", "run", "", "pair_id", "", ...
        "contact_a", "", "contact_b", "", "T", NaN, "signal", [], ...
        "bad_segments_s", zeros(0, 2));
    analysisPairs = repmat(template, height(sessionPairs), 1);
    for i = 1:height(sessionPairs)
        [aFound, ia] = ismember(double(sessionPairs.channel_a_index(i)), endpoint);
        [bFound, ib] = ismember(double(sessionPairs.channel_b_index(i)), endpoint);
        if ~aFound || ~bFound
            out = failed(sessionId, mode, "input_integrity_failed", "Endpoint lookup failed after HDF5 read.");
            return
        end
        br = bad(string(bad.session_id) == sessionId & bad.series == "LFP_macro" & ...
            ismember(string(bad.origchannel_name), [string(sessionPairs.contact_a(i)), string(sessionPairs.contact_b(i))]), :);
        iv = [double(br.start_common_time_s), double(br.stop_common_time_s)] - movieZeroCommon;
        if isempty(iv), iv = zeros(0, 2); end
        if ~isempty(iv)
            iv(:, 1) = max(iv(:, 1), 0);
            iv(:, 2) = min(iv(:, 2), T);
            iv = mergeIntervals(iv(iv(:, 2) > iv(:, 1), :));
        end
        analysisPairs(i) = struct("subject", string(sessionPairs.subject(i)), ...
            "session_id", sessionId, "run", string(sessionPairs.run(i)), ...
            "pair_id", string(sessionPairs.pair_name(i)), ...
            "contact_a", string(sessionPairs.contact_a(i)), ...
            "contact_b", string(sessionPairs.contact_b(i)), "T", T, ...
            "signal", endpoint250(ia, :)' - endpoint250(ib, :)', "bad_segments_s", iv);
    end

    out = struct("status", "complete", "mode", mode, "failure_reason", "", ...
        "subject", string(timingRow.subject(1)), "session_id", sessionId, ...
        "run", string(timingRow.run(1)), "nwb_file", string(nwbFile), ...
        "hdf5_path", hdf5Path, "orientation", string(contract.raw_lfp_binding.orientation), ...
        "movie_first_sample", first, "movie_last_sample", last, "raw_sample_count", n, ...
        "analysis_sample_count", size(endpoint250, 2), "sampling_rate_hz", 1000, ...
        "analysis_rate_hz", 250, "series_channel_indices", endpoint, ...
        "pairs", sessionPairs(:, ["subject" "session_id" "run" "pair_name" "contact_a" "contact_b" "channel_a_index" "channel_b_index" "hemisphere" "location"]), ...
        "analysis_pairs", analysisPairs, "diagnostics", struct());
catch ME
    out = failed(sessionId, mode, "input_integrity_failed", string(ME.identifier) + ": " + string(ME.message));
end
end

function out = failed(sessionId, mode, status, reason)
out = struct("status", string(status), "mode", string(mode), "failure_reason", string(reason), ...
    "subject", "", "session_id", string(sessionId), "run", "", "nwb_file", "", ...
    "hdf5_path", "", "orientation", "", "pairs", struct([]), ...
    "analysis_pairs", struct([]), "diagnostics", struct());
end

function out = mergeIntervals(iv)
if isempty(iv), out = zeros(0, 2); return, end
iv = sortrows(iv, [1 2]); out = iv(1, :);
for i = 2:size(iv, 1)
    if iv(i, 1) <= out(end, 2) + 1e-9
        out(end, 2) = max(out(end, 2), iv(i, 2));
    else
        out(end + 1, :) = iv(i, :); %#ok<AGROW>
    end
end
end

function t = readTsv(path)
t = readtable(path, "FileType", "text", "Delimiter", "\t", "TextType", "string", ...
    "VariableNamingRule", "preserve");
end
