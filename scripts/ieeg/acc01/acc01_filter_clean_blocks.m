function [filtered, diag] = acc01_filter_clean_blocks(signal, pair, contract)
%ACC01_FILTER_CLEAN_BLOCKS Zero-phase filter each continuous clean block.
% Bad samples remain zero in the returned trace. All event/support logic
% keeps them out of the PPC windows; the important invariant here is that a
% filtfilt call never receives samples from two clean blocks.
signal = double(signal(:));
fs = double(contract.sampling.analysis_hz);
n = numel(signal);
badMask = false(n, 1);
for k = 1:size(pair.bad_segments_s, 1)
    first = max(1, floor(double(pair.bad_segments_s(k, 1)) * fs) + 1);
    last = min(n, ceil(double(pair.bad_segments_s(k, 2)) * fs));
    if last >= first
        badMask(first:last) = true;
    end
end

cleanMask = ~badMask;
blockStart = find(cleanMask & [true; ~cleanMask(1:end-1)]);
blockStop = find(cleanMask & [~cleanMask(2:end); true]);
[b, a] = butter(double(contract.phase.filter_order), ...
    double(contract.phase.band_hz) / (fs / 2));
filtered = zeros(n, 1);
for k = 1:numel(blockStart)
    ix = blockStart(k):blockStop(k);
    % Let filtfilt report an invalid short clean block rather than silently
    % changing the frozen filter or borrowing samples across a bad segment.
    filtered(ix) = filtfilt(b, a, signal(ix));
end

diag = struct('n_samples', n, 'n_bad_samples', nnz(badMask), ...
    'n_blocks', numel(blockStart), 'block_start_index', blockStart, ...
    'block_stop_index', blockStop, 'bad_segments_s', pair.bad_segments_s, ...
    'filter_order', contract.phase.filter_order, ...
    'band_hz', contract.phase.band_hz, ...
    'scope', contract.phase.filter_scope);
end
