function report = validate_acc00_powercurve_profile_v2(root)
% Clean-process validator for v1/v2 profiling equivalence and checkpoint schema.
arguments
    root (1,1) string
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
paper=jsondecode(fileread(fullfile(root,"scripts","ieeg","paper","acc00_sim_powercurve_contract.json")));
v1=fullfile(root,"processed","subject","group","ieeg_methods_paper","p0_2_design","profiling");
v2=fullfile(root,"processed","subject","group","ieeg_methods_paper","p0_2_sim_powercurve","profiling","v2");
rows=cell(3,1); gids=[1 6 7];
for k=1:3
  gi=gids(k); a=load(fullfile(v1,sprintf("world_g%02d_n001.mat",gi))).checkpoint; b=load(fullfile(v2,sprintf("world_g%02d_n001.mat",gi))).checkpoint;
  assert(isequaln(a.subject_matrix,b.subject_matrix),"subject_matrix mismatch G%02d",gi);
  assert(isequaln(a.inference.R1.p_two_sided,b.inference.R1.p_two_sided)&&isequaln(a.inference.R2.p_two_sided,b.inference.R2.p_two_sided),"inference p mismatch G%02d",gi);
  required=["parent_historical_contract_sha256","input_hashes","beta_gen","expected_slope","counts","started_at_utc","completed_at_utc"];
  assert(all(isfield(b,cellstr(required))),"v2 schema missing G%02d",gi);
  assert(string(b.grid_id)==string(paper.additive.grid_id(gi))&&b.seed==91200000+1000*gi+1,"identity mismatch G%02d",gi);
  assert(isfinite(b.timing.peak_memory_bytes)&&b.timing.peak_memory_bytes>0,"memory instrumentation missing G%02d",gi);
  rows{k}=table(string(b.grid_id),b.seed,b.timing.elapsed_s,b.timing.peak_memory_bytes,b.beta_gen,b.expected_slope,"PASS",'VariableNames',["grid_id" "seed" "elapsed_s" "peak_memory_bytes" "beta_gen" "expected_slope" "status"]);
end
report=vertcat(rows{:}); out=fullfile(root,"processed","subject","group","ieeg_methods_paper","p0_2_sim_powercurve","profiling","v2_validator.tsv"); writetable(report,out,"FileType","text","Delimiter","\t"); disp(report);
end
