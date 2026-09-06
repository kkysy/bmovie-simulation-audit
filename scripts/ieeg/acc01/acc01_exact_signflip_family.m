function out = acc01_exact_signflip_family(x)
n = size(x, 1);
if n == 0
    out = struct('signs', [], 'statistics', [], 'observed', [NaN NaN], ...
        'p_two_sided', [NaN NaN], 'p_maxstat', NaN, ...
        'studentization_sd', [NaN NaN], 'studentized_statistics', [], ...
        'studentized_observed', [NaN NaN], 'max_statistics', [], ...
        'max_observed', NaN);
    return
end
bits = dec2bin(0:2^n-1) - '0';
signs = 1 - 2 * bits;
% dec2bin(0,...) is the all-plus observed row; retain the explicit
% assignment as a guard against any future enumeration-order change.
signs(1, :) = 1;
stats = (signs * x) / n;
obs = stats(1, :);
twoSided = mean(abs(stats) >= abs(obs), 1);

studentizationSd = std(stats, 0, 1);
studentizationScale = studentizationSd;
% A constant zero metric has no scale; preserving it as zero is the only
% finite edge behavior and does not affect non-degenerate family tests.
studentizationScale(studentizationScale == 0) = 1;
studentizedStats = stats ./ studentizationScale;
studentizedObs = studentizedStats(1, :);
maxStats = max(abs(studentizedStats), [], 2);
maxObs = max(abs(studentizedObs));
out = struct('signs', signs, 'statistics', stats, 'observed', obs, ...
    'p_two_sided', twoSided, 'p_maxstat', mean(maxStats >= maxObs), ...
    'studentization_sd', studentizationSd, ...
    'studentized_statistics', studentizedStats, ...
    'studentized_observed', studentizedObs, ...
    'max_statistics', maxStats, 'max_observed', maxObs);
end
