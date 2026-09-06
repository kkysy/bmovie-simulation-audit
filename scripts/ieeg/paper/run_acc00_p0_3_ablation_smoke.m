function report = run_acc00_p0_3_ablation_smoke(root, smokeLabel)
%RUN_ACC00_P0_3_ABLATION_SMOKE Frozen five-step clean-process smoke entrypoint.
arguments
    root (1,1) string
    smokeLabel (1,1) string = "smoke"
end
assert(~contains(smokeLabel,[filesep "/" "\\"]),"smokeLabel must be a single directory name.");root=string(char(java.io.File(char(root)).getCanonicalPath()));cd(root);addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(fullfile(root,"scripts","ieeg","paper"));paper=jsondecode(fileread(fullfile(root,"scripts","ieeg","paper","acc00_sim_p0_3_ablation_contract.json")));out=fullfile(root,string(paper.outputs.root),smokeLabel);if ~isfolder(out),mkdir(out),end;clock=tic;
fft=check_acc00_p0_3_p10_fft(root,out);
% P01 is deliberately run before the altered cells, and is resumed once to
% exercise atomic identity-verified resume without manufacturing a tmp file.
p01=run_acc00_p0_3_mpower_world(root,"P01",1,"smoke",out);p01resume=run_acc00_p0_3_mpower_world(root,"P01",1,"smoke",out);assert(p01.world_seed==p01resume.world_seed,"P01 resume identity failed.");
cells=["P02" "P03" "P04" "P07" "P10"];mpower=cell(numel(cells),1);for i=1:numel(cells),mpower{i}=run_acc00_p0_3_mpower_world(root,cells(i),1,"smoke",out);end
m01=run_acc00_p0_3_mphase_world(root,"M01",1,"smoke",out);m02=run_acc00_p0_3_mphase_world(root,"M02",1,"smoke",out);
validator=validate_acc00_p0_3_ablation(root,out,false);mpowerElapsed=cellfun(@(x)x.timing.elapsed_s,mpower);mpowerPeak=cellfun(@(x)x.timing.peak_memory_bytes,mpower);report=struct("status","PASS","fft",fft,"validator",validator,"elapsed_s",toc(clock),"mpower_world_elapsed_s",mpowerElapsed(:)',"mphase_world_elapsed_s",[m01.timing.elapsed_s m02.timing.elapsed_s],"max_peak_memory_bytes",max([p01.timing.peak_memory_bytes,mpowerPeak(:)',m01.timing.peak_memory_bytes,m02.timing.peak_memory_bytes]));fid=fopen(fullfile(out,"smoke_report.json"),"w");fprintf(fid,"%s\n",jsonencode(report,"PrettyPrint",true));fclose(fid);fprintf("P0-3 smoke PASS: %.3f s\n",report.elapsed_s);
end
