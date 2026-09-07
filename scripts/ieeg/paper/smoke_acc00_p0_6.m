function report=smoke_acc00_p0_6(root)
%SMOKE_ACC00_P0_6 One-shot smoke: dry-run 2 worlds per arm + assertions.
arguments
    root (1,1) string = pwd
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
fprintf("[smoke] bitwise background probe (CRN)...");
bitwiseBackgroundProbe(root);
fprintf(" PASS\n");
fprintf("[smoke] dry-run orchestration (2 worlds/arm)...\n");
run_acc00_p0_6_main(root,2,"dry-run");
fprintf("[smoke] summarize...\n");
summarize_acc00_p0_6(root,false);
fprintf("[smoke] validator...\n");
report=validate_acc00_p0_6_perturbation(root);
assert(string(report.status)=="PASS","P0-6 smoke validator FAIL.");
fprintf("[smoke] PASS (%d checks)\n",report.summary(2));
end

function bitwiseBackgroundProbe(root)
% Same seed => same RNG consumption for background synthesis (schedule-independent),
% so synthesized signals must be bitwise identical across schedules before any
% event-derived processing. Asserted on the first 8 pairs.
paperPath=fullfile(root,"scripts","ieeg","paper");addpath(fullfile(root,"scripts","ieeg"));addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(paperPath);
base=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
c=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_6_perturbation_contract.json")));
seed=20260827+100000+1;
sch15=fullfile(root,char(string(c.density_levels.("x1p5x").schedule_tsv)));sch05=fullfile(root,char(string(c.density_levels.("x0p5x").schedule_tsv)));
s1=base;s1.inputs.label_covariates.path=sch15;w1=acc00_sim_generate_world(root,s1,"benchmark",seed,"null",struct());
s0=base;s0.inputs.label_covariates.path=sch05;w0=acc00_sim_generate_world(root,s0,"benchmark",seed,"null",struct());
assert(numel(w1.pairs)==numel(w0.pairs),"Pair count mismatch.");
for i=1:8
    assert(isequal(w1.pairs(i).signal,w0.pairs(i).signal),"Background bitwise identity failed at pair %d.",i);
end
end
