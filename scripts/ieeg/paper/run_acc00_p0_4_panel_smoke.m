function report = run_acc00_p0_4_panel_smoke(root, smokeLabel)
%RUN_ACC00_P0_4_PANEL_SMOKE Synthetic P0-4 smoke; never launches 480 worlds.
arguments
    root (1,1) string
    smokeLabel (1,1) string = "smoke"
end
assert(~contains(smokeLabel,[filesep "/" "\\"]),"smokeLabel must be one directory name.");root=string(char(java.io.File(char(root)).getCanonicalPath()));addpath(fullfile(root,"scripts","ieeg","paper"));maxNumCompThreads(1);% pin single-thread from process start so M02 historical replay reproduces bitwise (do NOT run multithreaded worlds before it)
c=jsondecode(fileread(fullfile(root,"scripts","ieeg","paper","acc00_sim_p0_4_panel_contract.json")));out=fullfile(root,string(c.outputs.root),smokeLabel);if ~isfolder(out),mkdir(out),end;clock=tic;
p=run_acc00_p0_4_panel_world(root,"PanelMember-01",1,"smoke",out);m=run_acc00_p0_4_panel_world(root,"PanelMember-02",1,"smoke",out);d2=run_acc00_p0_4_panel_world(root,"PanelMember-03",1,"smoke",out);d1=run_acc00_p0_4_panel_world(root,"PanelMember-04",1,"smoke",out);
final=fullfile(out,"checkpoints","PanelMember-01","world_001.mat");tmp=final+".tmp.mat";movefile(final,tmp,"f");pResume=run_acc00_p0_4_panel_world(root,"PanelMember-01",1,"smoke",out);assert(p.world_seed==pResume.world_seed,"tmp resume identity failed.");d1Toy=check_acc00_p0_4_d1_apply_toy(root);history=check_acc00_p0_4_mref_historical(root,1);validator=validate_acc00_p0_4_panel(root,out,false);
report=struct("status","PASS","validator",validator,"d1_apply_toy",d1Toy,"historical_mref",history,"elapsed_s",toc(clock),"member_elapsed_s",[p.timing.elapsed_s m.timing.elapsed_s d2.timing.elapsed_s d1.timing.elapsed_s],"max_peak_memory_bytes",max([p.timing.peak_memory_bytes m.timing.peak_memory_bytes d2.timing.peak_memory_bytes d1.timing.peak_memory_bytes]));fid=fopen(fullfile(out,"smoke_report.json"),"w");fprintf(fid,"%s\n",jsonencode(report,"PrettyPrint",true));fclose(fid);fprintf("P0-4 smoke PASS: %.3f s\n",report.elapsed_s);
end
