% Release port of the analysis-workspace tool of the same name; see README in this folder.
function report = confirm_smoke(root)
%CONFIRM_SMOKE Minimal Stage A end-to-end smoke; no formal-world inference.
arguments
    root (1,1) string = string(pwd) % release port: run from the repository root
        % (the workspace original resolved the root from the tools/ layout)
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
toolDir=fileparts(mfilename("fullpath"));
stamp=string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmssSSS'Z'"));
out=fullfile(toolDir,"output","smoke_"+stamp);
mkdir(out);
c=jsondecode(fileread(fullfile(toolDir,"confirm_mphase_reducer_contract.json")));
seeds=double(c.independence.base)+double(c.independence.offset)+(double(c.independence.seed_index(1)):double(c.independence.seed_index(2)));
excluded=c.independence.excluded_seed_ranges;
for i=1:numel(excluded)
    assert(~any(seeds>=double(excluded(i).range(1)) & seeds<=double(excluded(i).range(2))), ...
        "Frozen seed range overlaps %s.",string(excluded(i).name));
end
toy=[0 0 9];assert(median(toy)==0 && mean(toy)==3,"Reducer distinction toy failed.");
q1=confirm_run_world(root,1,"smoke",out);
q2=confirm_run_world(root,1,"smoke",out);
assert(q1.world_seed==20860828 && q2.world_seed==q1.world_seed,"Smoke resume identity failed.");
assert(height(q1.pair_rows)==4 && all(q1.pair_rows.mphase_surrogate_reducer=="mean"), ...
    "Smoke production Mphase route failed.");
tmpFiles=dir(fullfile(out,"checkpoints",string(c.arm.id),"*.tmp.mat"));
assert(isempty(tmpFiles),"Smoke left tmp checkpoints.");
bridge=confirm_export_bridge(root,out,"smoke");
cmd1=sprintf('python "%s" "%s" --mode smoke',fullfile(toolDir,"summarize_confirm.py"),out);
[s1,t1]=system(cmd1,"-echo");assert(s1==0,"Smoke summary failed: %s",t1);
cmd2=sprintf('python "%s" "%s" --mode smoke',fullfile(toolDir,"validate_confirm.py"),out);
[s2,t2]=system(cmd2,"-echo");assert(s2==0,"Smoke validator failed: %s",t2);
validator=jsondecode(fileread(fullfile(out,"validator_report.json")));
decision=jsondecode(fileread(fullfile(out,"main_decision.json")));
report=struct("status","PASS","output",out,"world_seed",q1.world_seed, ...
    "pair_rows",height(q1.pair_rows),"subject_rows",height(q1.subject_rows), ...
    "bridge",bridge,"validator_status",string(validator.status), ...
    "decision_status",string(decision.status),"formal_worlds_run",0);
fid=fopen(fullfile(out,"smoke_report.json"),"w");assert(fid>=0);cl=onCleanup(@()fclose(fid));
fprintf(fid,"%s\n",jsonencode(report,"PrettyPrint",true));disp(jsonencode(report,"PrettyPrint",true));
end
