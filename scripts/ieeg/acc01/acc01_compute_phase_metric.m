function out = acc01_compute_phase_metric(pair, events, contract, stream)
%ACC01_COMPUTE_PHASE_METRIC Mphase with full-path refit empirical null.
domainClock = tic;
s = acc01_build_empirical_surrogates(pair, events, contract, stream);
domainSeconds = toc(domainClock);
emptyTiming = struct('crossfit_s', 0, 'surrogate_ppc_s', 0, 'domain_s', domainSeconds);
if isempty(s.delta_s)
    out = struct('difference', NaN, 'observed_ppc', NaN, 'surrogate_ppc', NaN, ...
        'surrogate', s, 'refit_delta_s', [], 'refit_ppc', [], ...
        'refit_surrogate_count', 0, 'diagnostics', struct(), 'timing', emptyTiming);
    return
end

events = carryFoldIdentity(events, contract);
obsClock = tic;
[obs, obsDiag] = acc01_phase_ppc_for_schedule(pair, events, contract);
obsSeconds = toc(obsClock);

% Keep the pre-existing 1024-point circular-correlation family as a cheap
% diagnostic only. It is not used to estimate the Mphase surrogate reducer.
[cheapObserved, cheapSurrogate] = corrSeriesDiagnostic(obsDiag.phase_trace, pair, events, contract, s.delta_s);

B = numel(s.delta_s);
J = min(double(contract.phase.refit_surrogates.J), B);
if B > J
    refitIndex = randperm(stream, B, J);
else
    refitIndex = 1:B;
end
refitDelta = double(s.delta_s(refitIndex));
refitPpc = nan(J, 1);
refitCrossfitSeconds = 0;
refitTotalSeconds = 0;
for j = 1:J
    shifted = events;
    shifted.label_time_s = mod(double(events.label_time_s) + refitDelta(j), double(pair.T));
    refitClock = tic;
    [refitPpc(j), refitDiag] = acc01_phase_ppc_for_schedule(pair, shifted, contract);
    elapsed = toc(refitClock);
    refitTotalSeconds = refitTotalSeconds + elapsed;
    refitCrossfitSeconds = refitCrossfitSeconds + refitDiag.crossfit.operator_seconds;
end

reducer = string(contract.phase.mphase_surrogate_reducer);
assert(reducer == "mean", "ACC00-SIM Mphase reducer must be frozen to mean.");
surrogatePpc = mean(refitPpc);
out = struct('difference', obs - surrogatePpc, 'observed_ppc', obs, ...
    'surrogate_ppc', surrogatePpc, 'surrogate', s, ...
    'refit_delta_s', refitDelta, 'refit_ppc', refitPpc, ...
    'refit_surrogate_count', J, ...
    'diagnostics', struct('reducer', reducer, ...
        'corr_series_cheap', struct('observed_ppc', cheapObserved, ...
            'surrogate_ppc', cheapSurrogate, 'count', numel(cheapSurrogate), ...
            'family_size', 1024, 'role', "diagnostic_only"), ...
        'refit_zero_is_explicitly_available', true, ...
        'observed_path', obsDiag), ...
    'timing', struct('crossfit_s', obsDiag.crossfit.operator_seconds + refitCrossfitSeconds, ...
        'surrogate_ppc_s', max(0, refitTotalSeconds - refitCrossfitSeconds), ...
        'domain_s', domainSeconds + obsSeconds));
end

function events = carryFoldIdentity(events, contract)
if ~ismember("crossfit_fold", string(events.Properties.VariableNames))
    n = height(events);
    events.crossfit_fold = min(contract.phase.erp_crossfit_folds, ...
        ceil((1:n)' * contract.phase.erp_crossfit_folds / max(n, 1)));
end
end

function [observed, surrogate] = corrSeriesDiagnostic(phase, pair, events, contract, deltas)
fs = double(contract.sampling.analysis_hz);
nGrid = min(round(double(pair.T) * fs), numel(phase));
n = height(events);
z = exp(1i * phase(1:nGrid));
impulse = zeros(nGrid, 1);
center = mod(round(double(events.label_time_s) * fs), nGrid) + 1;
for e = 1:numel(center)
    impulse(center(e)) = impulse(center(e)) + 1;
end
corrSeries = ifft(conj(fft(impulse)) .* fft(z));
ppcSeries = (abs(corrSeries).^2 - n) / (n * (n - 1));
responseOffsets = round(double(contract.windows_s.phase_response(1)) * fs): ...
    round(double(contract.windows_s.phase_response(2)) * fs);
observed = mean(ppcSeries(mod(responseOffsets, nGrid) + 1));
deltaIndex = round(double(deltas(:)) * fs);
surrogate = nan(numel(deltaIndex), 1);
for j = 1:numel(deltaIndex)
    surrogate(j) = mean(ppcSeries(mod(deltaIndex(j) + responseOffsets, nGrid) + 1));
end
end
