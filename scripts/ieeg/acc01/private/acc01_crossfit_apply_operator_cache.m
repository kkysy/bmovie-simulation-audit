function [residuals, diag] = acc01_crossfit_apply_operator_cache(pair, cache, centers, coefficientMode)
%ACC01_CROSSFIT_APPLY_OPERATOR_CACHE Refit and reassemble many schedules.
% Each column of centers is one schedule.  Event overlap subtraction remains
% sequential within each schedule, matching the frozen observed path.
% coefficientMode="refit" is the frozen default. "observed_fixed" is an
% explicit P0-4 diagnostic perturbation: B_obs,k is fitted only from the
% observed schedule and then applied to every shifted schedule.

if nargin < 4
    coefficientMode = "refit";
end
coefficientMode = string(coefficientMode);
assert(ismember(coefficientMode,["refit" "observed_fixed"]),"Unknown coefficient mode.");

centers = double(centers);
nSchedules = size(centers, 2);
n = cache.n_events;
residuals = repmat(pair.signal(:), 1, nSchedules);
applySeconds = 0;
observedB = cell(cache.n_folds, 1);
observedSolveCount = 0;
shiftSolveCount = 0;
if coefficientMode == "observed_fixed"
    observedIndices = centers(:, 1) + cache.off;
    observedEpochs = pair.signal(observedIndices);
    for k = 1:cache.n_folds
        tr = cache.train_mask{k};
        if nnz(tr) < size(cache.X, 2)
            continue
        end
        observedB{k} = mldivide(cache.operators{k}, observedEpochs(tr, :));
        observedSolveCount = observedSolveCount + 1;
    end
end

for j = 1:nSchedules
    indices = centers(:, j) + cache.off;
    E = pair.signal(indices);
    prediction = zeros(size(E));
    q = tic;
    for k = 1:cache.n_folds
        tr = cache.train_mask{k};
        te = cache.test_index{k};
        if nnz(tr) < size(cache.X, 2)
            continue
        end
        if coefficientMode == "observed_fixed"
            B = observedB{k};
        else
            B = mldivide(cache.operators{k}, E(tr, :));
            if j == 1
                observedSolveCount = observedSolveCount + 1;
            else
                shiftSolveCount = shiftSolveCount + 1;
            end
        end
        prediction(te, :) = cache.X(te, :) * B;
    end
    applySeconds = applySeconds + toc(q);

    r = pair.signal(:);
    for i = 1:n
        r(indices(i, :)) = r(indices(i, :)) - prediction(i, :)';
    end
    residuals(:, j) = r;
end

diag = struct("n_events", n, "n_offsets", numel(cache.off), ...
    "operator_seconds", applySeconds, "reference_relative_error", NaN, ...
    "operator_build_seconds", cache.operator_build_seconds, ...
    "schedule_count", nSchedules, "coefficient_mode", coefficientMode, ...
    "observed_coefficient_solves", observedSolveCount, ...
    "shift_epoch_coefficient_solves", shiftSolveCount);
end
