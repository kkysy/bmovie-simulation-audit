function [observed, surrogate] = acc01_corr_series_diagnostic( ...
    phase, pair, events, contract, deltas)
%ACC01_CORR_SERIES_DIAGNOSTIC Preserve the cheap 1024-point diagnostic.

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
