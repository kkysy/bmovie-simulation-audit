% Corner-only rerun after the 2026-09-10 G06 amplitude fix.
% run_acc00_p0_6_main must NOT be used: its resume check compares checkpoint
% contract_sha256 against the current contract file, and the 360 existing
% null-arm checkpoints record the pre-relativization hash A6222A61..., so main
% would discard and recompute them. This driver recomputes only the 50
% corner_probe_G06 worlds, then revalidates and re-aggregates.
% run from the repository root (release port: no absolute paths)
root = string(char(java.io.File(char(pwd)).getCanonicalPath()));
assert(isfile(fullfile(root,'MANIFEST.md')), 'run from the repository root');
addpath(fullfile(root,'scripts','ieeg','paper'));
pool = gcp('nocreate');
if isempty(pool), parpool('local', 12); end
pctRunOnAll(sprintf("addpath('%s');addpath('%s');addpath('%s');maxNumCompThreads(1);", ...
    char(fullfile(root,"scripts","ieeg")), char(fullfile(root,"scripts","ieeg","acc01")), ...
    char(fullfile(root,"scripts","ieeg","paper"))));
clock = tic;
parfor si = 1:50
    run_acc00_p0_6_world(root, "corner_probe_G06", si); %#ok<PFBNS>
end
fprintf("corner rerun wallclock %.1f s\n", toc(clock));
r = validate_acc00_p0_6_perturbation(root);
fprintf("validator %s (%d/%d)\n", r.status, r.summary(1), r.summary(2));
out = summarize_acc00_p0_6(root, true);
disp(out.verdict);
fprintf("CORNER RERUN COMPLETE\n");
