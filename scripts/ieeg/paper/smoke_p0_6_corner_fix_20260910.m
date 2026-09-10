% One-off smoke: verify the fixed corner probe injects the true G06 amplitude.
% run from the repository root (release port: no absolute paths)
root = string(char(java.io.File(char(pwd)).getCanonicalPath()));
assert(isfile(fullfile(root,'MANIFEST.md')), 'run from the repository root');
addpath(fullfile(root,'scripts','ieeg','paper'));
q = run_acc00_p0_6_world(root, "corner_probe_G06", 1);
fprintf("SMOKE grid_id=%s a_erp_uV=%.16g seed=%d elapsed=%.1f\n", ...
    q.corner.grid_id, q.corner.a_erp_uV, q.world_seed, q.timing.elapsed_s);
assert(q.corner.a_erp_uV == 4.61799281479342, "smoke: amplitude is not the G06 cell");
fprintf("SMOKE PASS\n");
