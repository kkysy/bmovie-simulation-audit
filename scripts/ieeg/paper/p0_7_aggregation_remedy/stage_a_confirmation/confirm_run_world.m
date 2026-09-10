function checkpoint = confirm_run_world(root, seedIndex, mode, outputRoot)
%CONFIRM_RUN_WORLD Run one matched median/mean Stage A null world.
arguments
    root (1,1) string
    seedIndex (1,1) double {mustBeInteger,mustBePositive}
    mode (1,1) string {mustBeMember(mode,["smoke" "benchmark"])} = "benchmark"
    outputRoot (1,1) string = fullfile(fileparts(mfilename("fullpath")),"output")
end
maxNumCompThreads(1);
root = canonicalPath(root);
toolDir = fileparts(mfilename("fullpath"));
addpath(fullfile(root,"scripts","ieeg","acc01"));
contractPath = fullfile(toolDir,"confirm_mphase_reducer_contract.json");
parentPath = fullfile(root,"scripts","ieeg","acc00_sim_contract.json");
c = jsondecode(fileread(contractPath));
base = jsondecode(fileread(parentPath));
validateFrozenContract(c,base,parentPath,root,seedIndex);
seed = double(c.independence.base) + double(c.independence.offset) + seedIndex;
arm = string(c.arm.id);
outDir = fullfile(outputRoot,"checkpoints",arm);
if ~isfolder(outDir), mkdir(outDir); end
final = fullfile(outDir,sprintf("world_%03d.mat",seedIndex));
tmp = final + ".tmp.mat";
if isfile(final)
    assert(~isfile(tmp),"Final checkpoint has stale tmp residue: %s",tmp);
    z = load(final,"checkpoint");
    assertIdentity(z.checkpoint,c,contractPath,parentPath,root,seedIndex,seed,mode);
    checkpoint = z.checkpoint;
    return
end
if isfile(tmp)
    z = load(tmp,"checkpoint");
    assertIdentity(z.checkpoint,c,contractPath,parentPath,root,seedIndex,seed,mode);
    movefile(tmp,final,"f");
    checkpoint = z.checkpoint;
    return
end

started = utcNow();
clock = tic;
world = acc00_sim_generate_world(root,base,mode,seed,"null",struct());
pairRows = repmat(pairTemplate(),numel(world.pairs),1);
for i = 1:numel(world.pairs)
    pair = world.pairs(i);
    [~,events] = acc01_build_support_and_events(pair,world.covariates,base);
    ph = acc01_compute_phase_metric_optimized(pair,events,base,world.stream);
    pairRows(i) = struct("subject",string(pair.subject),"run",string(pair.run), ...
        "session_id",string(pair.session_id),"pair_id",string(pair.pair_id), ...
        "corrected_mphase",double(ph.difference),"observed_ppc",double(ph.observed_ppc), ...
        "mean_refit_surrogate_ppc",double(ph.surrogate_ppc), ...
        "n_phase_events",height(events),"domain_fraction",double(ph.surrogate.domain_fraction), ...
        "surrogate_count",double(ph.surrogate.sample_count), ...
        "refit_surrogate_count",double(ph.refit_surrogate_count), ...
        "mphase_surrogate_reducer",string(ph.diagnostics.reducer));
end
[subjectRows,inference,runDiagnostics] = reduceBoth(pairRows);
checkpoint = struct("status","complete","analysis_id",string(c.analysis_id), ...
    "arm",arm,"density","1p0x","scenario","null","effect_injection","none", ...
    "mode",mode,"seed_index",seedIndex,"world_seed",seed, ...
    "contract_sha256",sha256File(contractPath), ...
    "parent_contract_sha256",sha256File(parentPath), ...
    "schedule_sha256",sha256File(fullfile(root,string(c.arm.schedule_path))), ...
    "input_hashes",inputHashes(base.inputs),"pair_rows",struct2table(pairRows), ...
    "subject_rows",subjectRows,"inference",inference,"run_diagnostics",runDiagnostics, ...
    "started_at_utc",started,"completed_at_utc",utcNow(), ...
    "worker_threads",maxNumCompThreads,"matlab_version",string(version), ...
    "timing",struct("elapsed_s",toc(clock),"synthesis_s",world.timing.synthesis_s, ...
        "peak_memory_bytes",world.timing.peak_memory_bytes));
save(tmp,"checkpoint","-v7.3");
movefile(tmp,final,"f");
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
        k = k + 1;
        if reducer=="median", statistic=median(vf); else, statistic=mean(vf); end
        out(k,:) = {keys.subject(i),keys.run(i),reducer,statistic,numel(v),numel(vf), ...
            pairSkew,mean(vf>0),nnz(vf==0)};
    end
end
t = cell2table(out,'VariableNames',["subject" "run" "reducer" "Mphase" ...
    "n_pairs_total" "n_pairs_finite" "pair_corrected_skewness" ...
    "pair_positive_fraction" "pair_zero_count"]);
inference = struct();
diagRows = cell(0,10);
for reducer = reducers
    for r = ["R1" "R2"]
        z = t(t.reducer==reducer & t.run==r,:);
        z = sortrows(z,"subject");
        x = double(z.Mphase);
        sf = acc01_exact_signflip_family([zeros(numel(x),1),x]);
        sf.subjects = string(z.subject);
        inference.(char(reducer)).(char(r)) = reduceSignflip(sf);
        diagRows(end+1,:) = {reducer,r,numel(x),finiteSkewness(x),nnz(x>0), ... %#ok<AGROW>
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

function x = finiteSkewness(v)
if numel(v)<3 || all(v==v(1)), x=NaN; else, x=skewness(v,0); end
end

function q = pairTemplate()
q = struct("subject","","run","","session_id","","pair_id","", ...
    "corrected_mphase",NaN,"observed_ppc",NaN,"mean_refit_surrogate_ppc",NaN, ...
    "n_phase_events",0,"domain_fraction",NaN,"surrogate_count",0, ...
    "refit_surrogate_count",0,"mphase_surrogate_reducer","");
end

function validateFrozenContract(c,base,parentPath,root,seedIndex)
assert(seedIndex>=double(c.independence.seed_index(1)) && seedIndex<=double(c.independence.seed_index(2)), ...
    "seed_index is outside the frozen 1:120 range.");
assert(string(c.parent_contract.sha256)==sha256File(parentPath),"Parent contract hash mismatch.");
assert(string(base.analysis_id)==string(c.parent_contract.required_analysis_id),"Parent analysis_id mismatch.");
assert(base.surrogate.domain_guard_edge_lower_s==1.5 && base.surrogate.domain_guard_edge_upper_s==1.5 ...
    && base.surrogate.domain_guard_bad_dilation_s==0.75,"Frozen support guards changed.");
assert(base.phase.refit_surrogates.J==32 && string(base.phase.mphase_surrogate_reducer)=="mean", ...
    "Frozen production Mphase surrogate definition changed.");
schedule = fullfile(root,string(c.arm.schedule_path));
assert(isfile(schedule) && sha256File(schedule)==string(c.arm.schedule_sha256),"Reference schedule identity mismatch.");
seed = double(c.independence.base)+double(c.independence.offset)+seedIndex;
excluded = c.independence.excluded_seed_ranges;
for i=1:numel(excluded)
    assert(seed<double(excluded(i).range(1)) || seed>double(excluded(i).range(2)), ...
        "World seed overlaps excluded range %s.",string(excluded(i).name));
end
end

function assertIdentity(q,c,contractPath,parentPath,root,seedIndex,seed,mode)
assert(string(q.status)=="complete" && string(q.analysis_id)==string(c.analysis_id) ...
    && string(q.arm)==string(c.arm.id) && string(q.mode)==mode ...
    && q.seed_index==seedIndex && q.world_seed==seed,"Checkpoint identity mismatch.");
assert(string(q.contract_sha256)==sha256File(contractPath),"Checkpoint contract hash mismatch.");
assert(string(q.parent_contract_sha256)==sha256File(parentPath),"Checkpoint parent hash mismatch.");
assert(string(q.schedule_sha256)==sha256File(fullfile(root,string(c.arm.schedule_path))), ...
    "Checkpoint schedule hash mismatch.");
end

function h = inputHashes(inputs)
h = struct();
for name = string(fieldnames(inputs))', h.(name)=upper(string(inputs.(name).sha256)); end
end

function s = utcNow(), s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")); end
function p = canonicalPath(p), p=string(char(java.io.File(char(p)).getCanonicalPath())); end
function h = sha256File(path)
fid=fopen(path,"rb");assert(fid>=0,"Cannot open %s",path);cl=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance("SHA-256");
while true, b=fread(fid,1048576,"*uint8");if isempty(b),break,end;md.update(typecast(b,"int8"));end
h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end
