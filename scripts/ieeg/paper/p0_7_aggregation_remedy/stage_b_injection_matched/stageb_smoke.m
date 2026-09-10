% Release port of the analysis-workspace tool of the same name; see README in this folder.
function report = stageb_smoke(root)
%STAGEB_SMOKE Two matched worlds (2 P0 + 2 P4); implementation only.
% Do not run before ZCode approval. Observed separation must never tune design.
arguments
    root (1,1) string = string(pwd) % release port: run from the repository root
        % (the workspace original resolved the root from the tools/ layout)
end
toolDir=fileparts(mfilename("fullpath"));stamp=string(datetime("now","TimeZone","UTC","Format","yyyyMMdd'T'HHmmssSSS'Z'"));
out=fullfile(toolDir,"output","smoke_"+stamp);mkdir(out);
q1=stageb_run_world(root,1,"smoke",out);q2=stageb_run_world(root,2,"smoke",out);q1r=stageb_run_world(root,1,"smoke",out);
assert(q1.world_seed==22360828&&q2.world_seed==22360829&&q1r.world_seed==q1.world_seed,"Smoke seed/resume failure.");
assert(numel(q1.arm_results)==2&&all(string({q1.arm_results.status})=="complete"),"Smoke matched arms incomplete.");
assert(isempty(dir(fullfile(out,"checkpoints","*.tmp.mat"))),"Smoke left tmp checkpoints.");
bridge=stageb_export_bridge(root,out,"smoke");
[s1,t1]=system(sprintf('python "%s" "%s" --mode smoke',fullfile(toolDir,"stageb_summarize.py"),out),"-echo");assert(s1==0,"Smoke summary failed: %s",t1);
[s2,t2]=system(sprintf('python "%s" "%s" --mode smoke',fullfile(toolDir,"stageb_validate.py"),out),"-echo");assert(s2==0,"Smoke validator failed: %s",t2);
report=struct("status","PASS","output",out,"matched_worlds",2,"arm_worlds",4,"bridge",bridge, ...
    "scientific_rules_evaluated",false,"tuning_prohibited",true,"formal_worlds_run",0);
writeJson(fullfile(out,"smoke_report.json"),report);disp(jsonencode(report,"PrettyPrint",true));
end
function writeJson(path,value),fid=fopen(path,"w");assert(fid>=0);cl=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",jsonencode(value,"PrettyPrint",true));end
