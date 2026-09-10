% Finish step for the 2026-09-10 corner rerun: worlds already computed; validate + aggregate.
% run from the repository root (release port: no absolute paths)
root = string(char(java.io.File(char(pwd)).getCanonicalPath()));
assert(isfile(fullfile(root,'MANIFEST.md')), 'run from the repository root');
addpath(fullfile(root,'scripts','ieeg','paper'));
r = validate_acc00_p0_6_perturbation(root);
fprintf("validator %s (%d/%d)\n", r.status, r.summary(1), r.summary(2));
out = summarize_acc00_p0_6(root, true);
disp(out.verdict);
fprintf("CORNER RERUN FINISH COMPLETE\n");
