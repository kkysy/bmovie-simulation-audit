function out = acc01_compute_phase_metric_optimized(pair, events, contract, stream, options)
%ACC01_COMPUTE_PHASE_METRIC_OPTIMIZED Exact-preserving batched Mphase path.
% The estimator, random stream, J=32 selection, mean reducer, and cheap
% diagnostic are unchanged.  Only schedule evaluation is batched.
%
% CrossfitCoefficientMode is a backwards-compatible P0-4 diagnostic hook.
% Its default ("refit") is the frozen behaviour used by every pre-existing
% caller.  "observed_fixed" is deliberately opt-in and is only for the
% P0-4 cross-fit-contamination perturbation.
arguments
    pair
    events
    contract
    stream
    options.CrossfitCoefficientMode (1,1) string {mustBeMember(options.CrossfitCoefficientMode,["refit" "observed_fixed"])} = "refit"
end

domainClock = tic;
s = acc01_build_empirical_surrogates(pair, events, contract, stream);
domainSeconds = toc(domainClock);
emptyTiming = struct("crossfit_s", 0, "surrogate_ppc_s", 0, ...
    "domain_s", domainSeconds);
if isempty(s.delta_s)
    out = struct("difference", NaN, "observed_ppc", NaN, ...
        "surrogate_ppc", NaN, "surrogate", s, "refit_delta_s", [], ...
        "refit_ppc", [], "refit_surrogate_count", 0, ...
        "diagnostics", struct(), "timing", emptyTiming);
    return
end

events = carryFoldIdentity(events, contract);
cacheClock = tic;
cache = acc01_crossfit_operator_cache(pair, events, contract);
cacheBuildSeconds = toc(cacheClock);

B = numel(s.delta_s);
J = min(double(contract.phase.refit_surrogates.J), B);
if B > J
    refitIndex = randperm(stream, B, J);
else
    refitIndex = 1:B;
end
refitDelta = double(s.delta_s(refitIndex));

fs = double(contract.sampling.analysis_hz);
n = height(events);
nSchedules = J + 1;
centers = nan(n, nSchedules);
centers(:, 1) = round(double(events.label_time_s) * fs) + 1;
for j = 1:J
    shiftedLabels = mod(double(events.label_time_s) + refitDelta(j), double(pair.T));
    centers(:, j + 1) = round(shiftedLabels * fs) + 1;
end

[allPpc, allDiag, batchTiming] = ...
    acc01_phase_ppc_for_schedules_optimized(pair, events, cache, centers, contract, options.CrossfitCoefficientMode);
obs = allPpc(1);
refitPpc = allPpc(2:end);
obsDiag = allDiag;

% Keep the pre-existing 1024-point circular-correlation family as a cheap
% diagnostic only; it is not used to estimate the Mphase surrogate reducer.
[cheapObserved, cheapSurrogate] = acc01_corr_series_diagnostic( ...
    obsDiag.phase_trace, pair, events, contract, s.delta_s);

reducer = string(contract.phase.mphase_surrogate_reducer);
assert(reducer == "mean", "ACC00-SIM Mphase reducer must be frozen to mean.");
surrogatePpc = mean(refitPpc);
out = struct("difference", obs - surrogatePpc, "observed_ppc", obs, ...
    "surrogate_ppc", surrogatePpc, "surrogate", s, ...
    "refit_delta_s", refitDelta, "refit_ppc", refitPpc, ...
    "refit_surrogate_count", J, ...
    "diagnostics", struct("reducer", reducer, ...
        "corr_series_cheap", struct("observed_ppc", cheapObserved, ...
            "surrogate_ppc", cheapSurrogate, "count", numel(cheapSurrogate), ...
            "family_size", 1024, "role", "diagnostic_only"), ...
        "refit_zero_is_explicitly_available", true, ...
        "observed_path", obsDiag), ...
    "timing", struct("crossfit_s", cacheBuildSeconds + batchTiming.crossfit_s, ...
        "surrogate_ppc_s", batchTiming.filter_ppc_s, ...
        "domain_s", domainSeconds, ...
        "operator_build_s", cacheBuildSeconds, ...
        "operator_apply_s", batchTiming.operator_apply_s));
end

function events = carryFoldIdentity(events, contract)
if ~ismember("crossfit_fold", string(events.Properties.VariableNames))
    n = height(events);
    events.crossfit_fold = min(contract.phase.erp_crossfit_folds, ...
        ceil((1:n)' * contract.phase.erp_crossfit_folds / max(n, 1)));
end
end
