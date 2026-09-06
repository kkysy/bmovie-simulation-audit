function [ppc, diag] = acc01_phase_ppc_for_schedule(pair, events, contract)
%ACC01_PHASE_PPC_FOR_SCHEDULE Full observed/refit phase path for one schedule.
% The event rows, including crossfit_fold when present, define the schedule.
[residual, crossfitDiag] = acc01_crossfit_erp_subtract(pair, events, contract);
[filtered, filterDiag] = acc01_filter_clean_blocks(residual, pair, contract);
phase = angle(hilbert(filtered));

fs = double(contract.sampling.analysis_hz);
nGrid = min(round(double(pair.T) * fs), numel(phase));
n = height(events);
responseOffsets = round(double(contract.windows_s.phase_response(1)) * fs): ...
    round(double(contract.windows_s.phase_response(2)) * fs);
if n < 2 || nGrid < 1
    ppc = NaN;
    ppcByOffset = NaN(size(responseOffsets));
else
    center = mod(round(double(events.label_time_s) * fs), nGrid) + 1;
    ppcByOffset = nan(size(responseOffsets));
    for k = 1:numel(responseOffsets)
        ix = mod(center - 1 + responseOffsets(k), nGrid) + 1;
        angles = phase(ix);
        ppcByOffset(k) = (abs(sum(exp(1i * angles)))^2 - n) / (n * (n - 1));
    end
    ppc = mean(ppcByOffset);
end

diag = struct('crossfit', crossfitDiag, 'filter', filterDiag, ...
    'phase_trace', phase, 'ppc_by_response_sample', ppcByOffset, ...
    'response_offsets_samples', responseOffsets);
end
