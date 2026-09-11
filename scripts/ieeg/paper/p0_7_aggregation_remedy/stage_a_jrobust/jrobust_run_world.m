function checkpoint = jrobust_run_world(root, seedIndex, runKind, outputRoot)
%JROBUST_RUN_WORLD Run one matched J32/J256 full null world.
arguments
    root (1,1) string
    seedIndex (1,1) double {mustBeInteger,mustBePositive}
    runKind (1,1) string {mustBeMember(runKind,["smoke" "formal"])}
    outputRoot (1,1) string = fullfile(fileparts(mfilename("fullpath")),"output")
end
maxNumCompThreads(1);
root = canonicalPath(root);
toolDir = fileparts(mfilename("fullpath"));
addpath(fullfile(root,"scripts","ieeg","acc01"));
contractPath = fullfile(toolDir,"confirm_mphase_jrobust_contract.json");
parentPath = fullfile(root,"scripts","ieeg","acc00_sim_contract.json");
c = jsondecode(fileread(contractPath));
base = jsondecode(fileread(parentPath));
validateFrozenContract(c,base,root,contractPath,parentPath,seedIndex,runKind);

worldSeed = double(c.seed_namespace.base)+double(c.seed_namespace.offset)+seedIndex;
outDir = fullfile(outputRoot,"checkpoints",string(c.world.arm_id));
if ~isfolder(outDir), mkdir(outDir); end
final = fullfile(outDir,sprintf("world_%03d.mat",seedIndex));
tmp = final+".tmp.mat";
if isfile(final)
    assert(~isfile(tmp),"Complete checkpoint has stale tmp residue: %s",tmp);
    z = load(final,"checkpoint");
    assertIdentity(z.checkpoint,c,root,contractPath,parentPath,seedIndex,worldSeed,runKind);
    checkpoint = z.checkpoint;
    return
end
if isfile(tmp)
    z = load(tmp,"checkpoint");
    assertIdentity(z.checkpoint,c,root,contractPath,parentPath,seedIndex,worldSeed,runKind);
    movefile(tmp,final,"f");
    checkpoint = z.checkpoint;
    return
end

started = utcNow();
worldClock = tic;
world = acc00_sim_generate_world(root,base,string(c.world.generator_mode),worldSeed,"null",struct());
nPairs = numel(world.pairs);
eventsByPair = cell(nPairs,1);
initialStates = cell(nPairs,1);
j32Rows = repmat(pairTemplate(),nPairs,1);
for i = 1:nPairs
    pair = world.pairs(i);
    [~,events] = acc01_build_support_and_events(pair,world.covariates,base);
    eventsByPair{i} = events;
    initialStates{i} = world.stream.State;
    pairClock = tic;
    ph = acc01_compute_phase_metric_optimized(pair,events,base,world.stream);
    j32Rows(i) = makePairRow(pair,events,ph,"J32_reference",32,initialStates{i},world.stream.State,toc(pairClock));
end
[j32Subjects,j32Inference,j32Diagnostics] = reduceBoth(j32Rows);
j32 = makeArmResult("J32_reference",32,j32Rows,j32Subjects,j32Inference,j32Diagnostics,base);
assertStageAReference(root,c,seedIndex,worldSeed,j32);

high = base;
high.phase.refit_surrogates.J = double(c.arms.J256_high.J);
assertOnlyJDiff(base,high,32,double(c.arms.J256_high.J));
j256Rows = repmat(pairTemplate(),nPairs,1);
for i = 1:nPairs
    pair = world.pairs(i);
    metricStream = RandStream("Threefry","Seed",worldSeed);
    metricStream.State = initialStates{i};
    pairClock = tic;
    ph = acc01_compute_phase_metric_optimized(pair,eventsByPair{i},high,metricStream);
    j256Rows(i) = makePairRow(pair,eventsByPair{i},ph,"J256_high",256,initialStates{i},metricStream.State,toc(pairClock));
end
[j256Subjects,j256Inference,j256Diagnostics] = reduceBoth(j256Rows);
j256 = makeArmResult("J256_high",256,j256Rows,j256Subjects,j256Inference,j256Diagnostics,high);
assertMatchedPairCRN(j32,j256);

armResults = [j32;j256];
checkpoint = struct("status","complete","analysis_id",string(c.analysis_id), ...
    "run_kind",runKind,"generator_mode",string(c.world.generator_mode), ...
    "arm_id",string(c.world.arm_id),"seed_index",seedIndex,"world_seed",worldSeed, ...
    "contract_sha256",sha256File(contractPath),"parent_contract_sha256",sha256File(parentPath), ...
    "stage_a_contract_sha256",sha256File(fullfile(root,string(c.stage_a_reference.contract.path))), ...
    "schedule_sha256",sha256File(fullfile(root,string(c.world.schedule_path))), ...
    "input_hashes",inputHashes(base.inputs),"production_code_hashes",fileHashRows(root,c.production_code), ...
    "local_code_hashes",fileHashRows(root,c.local_code),"background_world_sha256",worldFingerprint(j32Rows), ...
    "effective_contract_diff",string(c.arms.allowed_effective_contract_diff), ...
    "arm_order",string(c.arms.order),"arm_results",armResults, ...
    "started_at_utc",started,"completed_at_utc",utcNow(), ...
    "worker_threads",maxNumCompThreads,"matlab_version",string(version), ...
    "timing",struct("elapsed_s",toc(worldClock),"synthesis_s",world.timing.synthesis_s, ...
        "peak_memory_bytes",world.timing.peak_memory_bytes));
save(tmp,"checkpoint","-v7.3");
movefile(tmp,final,"f");
end

function row = makePairRow(pair,events,ph,arm,J,initialState,finalState,elapsed)
B = double(ph.surrogate.sample_count);
if B<=J, selectionMode="all_candidates"; else, selectionMode="uniform_without_replacement"; end
row = struct("arm",arm,"J_contract",J,"subject",string(pair.subject), ...
    "run",string(pair.run),"session_id",string(pair.session_id),"pair_id",string(pair.pair_id), ...
    "corrected_mphase",double(ph.difference),"observed_ppc",double(ph.observed_ppc), ...
    "mean_refit_surrogate_ppc",double(ph.surrogate_ppc),"n_phase_events",height(events), ...
    "quantized_feasible_count",double(ph.surrogate.quantized_feasible_count), ...
    "domain_fraction",double(ph.surrogate.domain_fraction),"surrogate_count",B, ...
    "surrogate_family_mode",string(ph.surrogate.mode),"refit_selection_mode",selectionMode, ...
    "refit_surrogate_count",double(ph.refit_surrogate_count), ...
    "mphase_surrogate_reducer",string(ph.diagnostics.reducer), ...
    "candidate_delta_sha256",numericHash(double(ph.surrogate.delta_s)), ...
    "refit_delta_sha256",numericHash(double(ph.refit_delta_s)), ...
    "pair_stream_initial_state_sha256",numericHash(initialState), ...
    "pair_stream_final_state_sha256",numericHash(finalState), ...
    "pair_signal_sha256",numericHash(double(pair.signal)), ...
    "pair_elapsed_s",elapsed,"surrogate_ppc_s",double(ph.timing.surrogate_ppc_s));
end

function arm = makeArmResult(id,J,pairRows,subjectRows,inference,diagnostics,effective)
arm = struct("status","complete","arm",id,"J_contract",J, ...
    "effective_contract_sha256",sha256Text(jsonencode(effective)), ...
    "pair_rows",struct2table(pairRows),"subject_rows",subjectRows, ...
    "inference",inference,"run_diagnostics",diagnostics);
end

function [t,inference,diagTable] = reduceBoth(rows)
subject = string({rows.subject})';
run = string({rows.run})';
value = double([rows.corrected_mphase])';
keys = unique(table(subject,run),'rows','sorted');
reducers = ["median" "mean"];
out = cell(height(keys)*numel(reducers),9);
k = 0;
for i = 1:height(keys)
    keep = subject==keys.subject(i) & run==keys.run(i);
    v = value(keep);
    finite = isfinite(v);
    vf = v(finite);
    assert(~isempty(vf),"Subject %s %s has no finite Mphase pairs.",keys.subject(i),keys.run(i));
    pairSkew = finiteSkewness(vf);
    for reducer = reducers
        k = k+1;
        if reducer=="median", statistic=median(vf); else, statistic=mean(vf); end
        out(k,:) = {keys.subject(i),keys.run(i),reducer,statistic,numel(v),numel(vf), ...
            pairSkew,mean(vf>0),nnz(vf==0)};
    end
end
t = cell2table(out,'VariableNames',["subject" "run" "reducer" "Mphase" ...
    "n_pairs_total" "n_pairs_finite" "pair_corrected_skewness" ...
    "pair_positive_fraction" "pair_zero_count"]);
inference = struct();
diagRows = cell(numel(reducers)*2,10);
diagIndex = 0;
for reducer = reducers
    for r = ["R1" "R2"]
        z = sortrows(t(t.reducer==reducer & t.run==r,:),"subject");
        x = double(z.Mphase);
        sf = acc01_exact_signflip_family([zeros(numel(x),1),x]);
        sf.subjects = string(z.subject);
        inference.(char(reducer)).(char(r)) = reduceSignflip(sf);
        diagIndex = diagIndex+1;
        diagRows(diagIndex,:) = {reducer,r,numel(x),finiteSkewness(x),nnz(x>0), ...
            nnz(x<0),nnz(x==0),mean(x>0),mean(x),std(x,0)};
    end
end
diagTable = cell2table(diagRows,'VariableNames',["reducer" "run" "n_subjects" ...
    "subject_statistic_skewness" "n_positive" "n_negative" "n_zero" ...
    "positive_fraction" "subject_mean" "subject_sample_sd"]);
end

function s = reduceSignflip(q)
s = struct("subjects",q.subjects,"n_subjects",numel(q.subjects), ...
    "observed",q.observed,"p_two_sided",q.p_two_sided, ...
    "p_maxstat",q.p_maxstat,"studentization_sd",q.studentization_sd, ...
    "studentized_observed",q.studentized_observed,"max_observed",q.max_observed);
end

function assertStageAReference(root,c,seedIndex,worldSeed,arm)
path = fullfile(root,sprintf(string(c.stage_a_reference.checkpoint_pattern),seedIndex));
assert(isfile(path),"Missing frozen Stage A checkpoint: %s",path);
z = load(path,"checkpoint");
q = z.checkpoint;
assert(string(q.status)=="complete" && string(q.analysis_id)==string(c.stage_a_reference.analysis_id) ...
    && q.seed_index==seedIndex && q.world_seed==worldSeed && string(q.mode)=="benchmark", ...
    "Stage A checkpoint identity mismatch.");
old = q.pair_rows;
new = arm.pair_rows;
assert(height(old)==height(new),"Stage A pair count mismatch.");
textFields = ["subject" "run" "session_id" "pair_id" "mphase_surrogate_reducer"];
for f = textFields, assert(all(string(old.(f))==string(new.(f))),"Stage A pair identity mismatch: %s",f); end
numericFields = ["corrected_mphase" "observed_ppc" "mean_refit_surrogate_ppc" ...
    "n_phase_events" "domain_fraction" "surrogate_count" "refit_surrogate_count"];
for f = numericFields
    assert(allClose(double(old.(f)),double(new.(f)),double(c.tolerances.stage_a_continuous_abs)), ...
        "Stage A pair value mismatch: %s",f);
end
assertTableMatch(q.subject_rows,arm.subject_rows,double(c.tolerances.stage_a_continuous_abs));
for reducer = ["median" "mean"]
    for run = ["R1" "R2"]
        a = q.inference.(char(reducer)).(char(run));
        b = arm.inference.(char(reducer)).(char(run));
        assert(allClose(double(a.observed),double(b.observed),double(c.tolerances.stage_a_continuous_abs)) ...
            && allClose(double(a.p_two_sided),double(b.p_two_sided),double(c.tolerances.stage_a_continuous_abs)) ...
            && allClose(double(a.studentization_sd),double(b.studentization_sd),double(c.tolerances.stage_a_continuous_abs)), ...
            "Stage A inference mismatch for %s/%s.",reducer,run);
    end
end
end

function assertTableMatch(old,new,tol)
assert(height(old)==height(new),"Stage A subject count mismatch.");
for f = ["subject" "run" "reducer"], assert(all(string(old.(f))==string(new.(f))),"Stage A subject identity mismatch."); end
for f = ["Mphase" "n_pairs_total" "n_pairs_finite" "pair_corrected_skewness" "pair_positive_fraction" "pair_zero_count"]
    assert(allClose(double(old.(f)),double(new.(f)),tol),"Stage A subject value mismatch: %s",f);
end
end

function assertMatchedPairCRN(j32,j256)
a = j32.pair_rows;
b = j256.pair_rows;
assert(height(a)==height(b),"Matched arm pair count mismatch.");
for f = ["subject" "run" "session_id" "pair_id" "candidate_delta_sha256" ...
        "pair_stream_initial_state_sha256" "pair_signal_sha256"]
    assert(all(string(a.(f))==string(b.(f))),"Matched pair CRN mismatch: %s",f);
end
for f = ["n_phase_events" "quantized_feasible_count" "domain_fraction" "surrogate_count"]
    assert(allClose(double(a.(f)),double(b.(f)),1e-12),"Matched candidate-family mismatch: %s",f);
end
assert(all(double(a.refit_surrogate_count)==min(32,double(a.surrogate_count))),"J32 count rule failed.");
assert(all(double(b.refit_surrogate_count)==min(256,double(b.surrogate_count))),"J256 count rule failed.");
end

function validateFrozenContract(c,base,root,contractPath,parentPath,seedIndex,runKind)
assert(logical(c.authorization.contract_and_code) && logical(c.authorization.static_review), ...
    "Contract/code authorization is missing.");
if runKind=="smoke"
    assert(logical(c.authorization.smoke_2_world),"Two-world smoke is not authorized by the child contract.");
else
    assert(logical(c.authorization.formal_120_world),"Formal 120-world run is not authorized by the child contract.");
end
assert(seedIndex>=double(c.seed_namespace.seed_index(1)) && seedIndex<=double(c.seed_namespace.seed_index(2)), ...
    "seed_index outside frozen range.");
assert(string(c.parent_contract.sha256)==sha256File(parentPath) ...
    && string(base.analysis_id)==string(c.parent_contract.required_analysis_id),"Parent identity mismatch.");
assert(double(base.phase.refit_surrogates.J)==double(c.parent_contract.required_J) ...
    && string(base.phase.mphase_surrogate_reducer)==string(c.parent_contract.required_surrogate_reducer), ...
    "Parent Mphase definition mismatch.");
assert(string(c.arms.order(1))=="J32_reference" && string(c.arms.order(2))=="J256_high", ...
    "Arm order mismatch.");
assert(double(c.arms.J32_reference.J)==32 && double(c.arms.J256_high.J)==256, ...
    "Frozen J values changed.");
assert(isequal(string(c.arms.allowed_effective_contract_diff),"phase.refit_surrogates.J"), ...
    "Effective-contract allowlist changed.");
assert(sha256File(fullfile(root,string(c.world.schedule_path)))==string(c.world.schedule_sha256), ...
    "Schedule identity mismatch.");
verifyPathHash(root,c.stage_a_reference.contract);
verifyPathHash(root,c.stage_a_reference.runner);
verifyPathHash(root,c.stage_a_reference.validator);
verifyPathHash(root,c.stage_a_reference.pair_rows);
verifyPathHash(root,c.stage_a_reference.subject_rows);
verifyPathHash(root,c.stage_a_reference.world_rows);
verifyFileRows(root,c.production_code);
verifyFileRows(root,c.local_code);
assert(isfile(contractPath),"Missing child contract.");
end

function assertOnlyJDiff(parent,derived,parentJ,derivedJ)
assert(double(parent.phase.refit_surrogates.J)==parentJ && double(derived.phase.refit_surrogates.J)==derivedJ, ...
    "Unexpected J values in effective contracts.");
restored = derived;
restored.phase.refit_surrogates.J = parentJ;
assert(isequaln(restored,parent),"Effective contract differs from parent outside phase.refit_surrogates.J.");
end

function assertIdentity(q,c,root,contractPath,parentPath,seedIndex,worldSeed,runKind)
assert(string(q.status)=="complete" && string(q.analysis_id)==string(c.analysis_id) ...
    && string(q.run_kind)==runKind && string(q.generator_mode)==string(c.world.generator_mode) ...
    && q.seed_index==seedIndex && q.world_seed==worldSeed,"Checkpoint identity mismatch.");
assert(string(q.contract_sha256)==sha256File(contractPath) ...
    && string(q.parent_contract_sha256)==sha256File(parentPath),"Checkpoint contract hash mismatch.");
assert(string(q.schedule_sha256)==sha256File(fullfile(root,string(c.world.schedule_path))), ...
    "Checkpoint schedule hash mismatch.");
assert(isequal(string(q.arm_order),string(c.arms.order)) && numel(q.arm_results)==2 ...
    && all(string({q.arm_results.status})=="complete"),"Checkpoint matched arms incomplete.");
assertMatchedPairCRN(q.arm_results(1),q.arm_results(2));
end

function verifyPathHash(root,item)
path = fullfile(root,string(item.path));
assert(isfile(path) && sha256File(path)==string(item.sha256),"Frozen file identity mismatch: %s",path);
end

function verifyFileRows(root,rows)
for i = 1:numel(rows), verifyPathHash(root,rows(i)); end
end

function t = fileHashRows(root,rows)
paths = strings(numel(rows),1);
hashes = strings(numel(rows),1);
for i = 1:numel(rows)
    paths(i) = string(rows(i).path);
    hashes(i) = sha256File(fullfile(root,paths(i)));
end
t = table(paths,hashes,'VariableNames',["path" "sha256"]);
end

function tf = allClose(a,b,tol)
sameFinite = isfinite(a) & isfinite(b) & abs(a-b)<=tol;
sameNaN = isnan(a) & isnan(b);
sameInf = isinf(a) & isinf(b) & a==b;
tf = all((sameFinite | sameNaN | sameInf),"all");
end

function x = finiteSkewness(v)
if numel(v)<3 || all(v==v(1)), x=NaN; else, x=skewness(v,0); end
end

function q = pairTemplate()
q = struct("arm","","J_contract",0,"subject","","run","","session_id","","pair_id","", ...
    "corrected_mphase",NaN,"observed_ppc",NaN,"mean_refit_surrogate_ppc",NaN, ...
    "n_phase_events",0,"quantized_feasible_count",0,"domain_fraction",NaN,"surrogate_count",0, ...
    "surrogate_family_mode","","refit_selection_mode","","refit_surrogate_count",0, ...
    "mphase_surrogate_reducer","","candidate_delta_sha256","","refit_delta_sha256","", ...
    "pair_stream_initial_state_sha256","","pair_stream_final_state_sha256","", ...
    "pair_signal_sha256","","pair_elapsed_s",NaN,"surrogate_ppc_s",NaN);
end

function h = worldFingerprint(rows)
text = join(string({rows.pair_id})'+"="+string({rows.pair_signal_sha256})',"|");
h = sha256Text(text);
end

function h = inputHashes(inputs)
h = struct();
for name = string(fieldnames(inputs))', h.(name)=upper(string(inputs.(name).sha256)); end
end

function h = numericHash(value)
md = java.security.MessageDigest.getInstance("SHA-256");
meta = uint8(char(string(class(value))+"|"+join(string(size(value)),"x")+"|"));
md.update(typecast(meta,"int8"));
if islogical(value), bytes=uint8(value(:)); else, bytes=typecast(value(:),"uint8"); end
if ~isempty(bytes), md.update(typecast(bytes,"int8")); end
h = upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end

function h = sha256Text(text)
md = java.security.MessageDigest.getInstance("SHA-256");
bytes = unicode2native(char(text),"UTF-8");
md.update(typecast(uint8(bytes),"int8"));
h = upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end

function h = sha256File(path)
fid=fopen(path,"rb");assert(fid>=0,"Cannot open %s",path);cleanup=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance("SHA-256");
while true, bytes=fread(fid,1048576,"*uint8");if isempty(bytes),break,end;md.update(typecast(bytes,"int8"));end
h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end

function s = utcNow()
s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end

function p = canonicalPath(p)
p=string(char(java.io.File(char(p)).getCanonicalPath()));
end
