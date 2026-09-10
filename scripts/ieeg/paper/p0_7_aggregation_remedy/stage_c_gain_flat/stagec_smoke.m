% Release port of the analysis-workspace tool of the same name; see README in this folder.
function report = stagec_smoke(root)
%STAGEC_SMOKE Two gain-flattened worlds; implementation only.
% Do not run before ZCode approval. Observed separation must never tune design.
arguments
    root (1,1) string = string(pwd) % release port: run from the repository root
        % (the workspace original resolved the root from the tools/ layout)
end
toolDir=fileparts(mfilename("fullpath"));stamp=string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmssSSS'Z'"));
out=fullfile(toolDir,"output","smoke_"+stamp);mkdir(out);
q1=stagec_run_world(root,1,"smoke",out);q2=stagec_run_world(root,2,"smoke",out);q1r=stagec_run_world(root,1,"smoke",out);
assert(q1.world_seed==22360828&&q2.world_seed==22360829&&q1r.world_seed==q1.world_seed,"Smoke seed/resume failure.");
assert(numel(q1.arm_results)==1&&string(q1.arm_results.status)=="complete","Smoke arm incomplete.");
assert(isempty(dir(fullfile(out,"checkpoints","*.tmp.mat"))),"Smoke left tmp checkpoints.");
bridge=stagec_export_bridge(root,out,"smoke");
[s1,t1]=system(sprintf('python "%s" "%s" --mode smoke',fullfile(toolDir,"stagec_summarize.py"),out),"-echo");assert(s1==0,"Smoke summary failed: %s",t1);
[s2,t2]=system(sprintf('python "%s" "%s" --mode smoke',fullfile(toolDir,"stagec_validate.py"),out),"-echo");assert(s2==0,"Smoke validator failed: %s",t2);
report=struct("status","PASS","output",out,"worlds",2,"bridge",bridge, ...
    "scientific_rules_evaluated",false,"tuning_prohibited",true,"formal_worlds_run",0);
fid=fopen(fullfile(out,"smoke_report.json"),"w");assert(fid>=0);cl=onCleanup(@()fclose(fid));
fprintf(fid,"%s\n",jsonencode(report,"PrettyPrint",true));clear cl;disp(jsonencode(report,"PrettyPrint",true));
end
