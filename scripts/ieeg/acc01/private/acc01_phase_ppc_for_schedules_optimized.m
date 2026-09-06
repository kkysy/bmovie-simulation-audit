function [ppc, diag, timing] = acc01_phase_ppc_for_schedules_optimized( ...
    pair, events, cache, centers, contract, coefficientMode)
%ACC01_PHASE_PPC_FOR_SCHEDULES_OPTIMIZED Evaluate observed plus refit paths.

if nargin < 6
    coefficientMode = "refit";
end

q = tic;
[residuals, crossfitDiag] = acc01_crossfit_apply_operator_cache(pair, cache, centers, coefficientMode);
crossfitSeconds = toc(q);

q = tic;
[filtered, filterDiag] = acc01_filter_clean_blocks_matrix(residuals, pair, contract);
phase = angle(hilbert(filtered));

fs = double(contract.sampling.analysis_hz);
nGrid = min(round(double(pair.T) * fs), size(phase, 1));
nSamples = size(phase, 1);
n = height(events);
responseOffsets = round(double(contract.windows_s.phase_response(1)) * fs): ...
    round(double(contract.windows_s.phase_response(2)) * fs);
nSchedules = size(centers, 2);
if n < 2 || nGrid < 1
    ppc = NaN(nSchedules, 1);
    ppcByOffset = NaN(nSchedules, numel(responseOffsets));
else
    % Put events in the first (contiguous) dimension.  This preserves the
    % scalar sum order for each schedule/offset while eliminating the 101 x
    % 32 point loops.
    ppcCenters = mod(centers(1:n, :) - 1, nGrid) + 1;
    positions = mod(ppcCenters - 1 + ...
        reshape(responseOffsets, 1, 1, []), nGrid) + 1;
    % Column strides use the actual matrix height.  nGrid can be one sample
    % shorter than the stored trace when the generator uses ceil(T*fs).
    linearIndex = positions + reshape((0:nSchedules-1) * nSamples, 1, [], 1);
    angles = phase(linearIndex);
    ppcByOffset = reshape((abs(sum(exp(1i * angles), 1)).^2 - n) ./ ...
        (n * (n - 1)), nSchedules, numel(responseOffsets));
    ppc = mean(ppcByOffset, 2);
end
filterPpcSeconds = toc(q);

% Return the observed-path diagnostics in the same schema as the scalar
% helper.  The extra schedule columns stay internal to this optimized path.
obsDiag = struct("crossfit", crossfitDiag, "filter", filterDiag, ...
    "phase_trace", phase(:, 1), "ppc_by_response_sample", ppcByOffset(1, :), ...
    "response_offsets_samples", responseOffsets);
diag = obsDiag;
timing = struct("crossfit_s", crossfitSeconds, ...
    "filter_ppc_s", filterPpcSeconds, ...
    "operator_apply_s", crossfitDiag.operator_seconds, ...
    "operator_build_s", cache.operator_build_seconds, ...
    "schedule_count", nSchedules);
end
