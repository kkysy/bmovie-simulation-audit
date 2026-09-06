function cache = acc01_crossfit_operator_cache(~, events, contract)
%ACC01_CROSSFIT_OPERATOR_CACHE Cache the fixed cross-fit design operators.
% The event attributes and fold identity are unchanged by a whole-schedule
% shift.  Only the sampled epoch matrix E changes, so the QR decomposition
% of each training design matrix can be reused for all schedules.

fs = double(contract.sampling.analysis_hz);
off = round(double(contract.windows_s.filter_support(1)) * fs): ...
    round(double(contract.windows_s.filter_support(2)) * fs);
n = height(events);
if ismember("crossfit_fold", string(events.Properties.VariableNames))
    fold = double(events.crossfit_fold(:));
    assert(numel(fold) == n && all(isfinite(fold)) && ...
        all(fold >= 1 & fold <= contract.phase.erp_crossfit_folds), ...
        "crossfit_fold must carry a valid frozen fold identity.");
else
    fold = min(contract.phase.erp_crossfit_folds, ...
        ceil((1:n)' * contract.phase.erp_crossfit_folds / max(n, 1)));
end

cols = endsWith(events.Properties.VariableNames, "_delta_z") | ...
    endsWith(events.Properties.VariableNames, "movie_time_linear_z") | ...
    endsWith(events.Properties.VariableNames, "movie_time_quadratic_z");
X = [ones(n, 1), double(events.z_e), double(events{:, cols})];
center = round(double(events.label_time_s) * fs) + 1;
nFolds = max(fold);
operators = cell(nFolds, 1);
trainMask = cell(nFolds, 1);
testIndex = cell(nFolds, 1);
buildSeconds = 0;
for k = 1:nFolds
    tr = fold ~= k;
    te = find(fold == k);
    trainMask{k} = tr;
    testIndex{k} = te;
    if nnz(tr) < size(X, 2)
        continue
    end
    q = tic;
    % decomposition(...,'qr') and mldivide reproduce X(tr,:)\E(tr,:)
    % while retaining the factorization for every shifted schedule.
    operators{k} = decomposition(X(tr, :), "qr");
    buildSeconds = buildSeconds + toc(q);
end

cache = struct("X", X, "fold", fold, "center", center, "off", off, ...
    "n_events", n, "n_folds", nFolds, "operators", {operators}, ...
    "train_mask", {trainMask}, "test_index", {testIndex}, ...
    "operator_build_seconds", buildSeconds);
end
