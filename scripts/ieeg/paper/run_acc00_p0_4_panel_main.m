function manifest = run_acc00_p0_4_panel_main(root, workers, authorize, outputRoot, dryRun)
%RUN_ACC00_P0_4_PANEL_MAIN Gated P0-4 panel orchestration and manifest writer.
% Formal mode can only run after a future contract/user authorization.  The
% explicit dry-run mode is limited to seeds 1:2 and can never target the
% production output root.
arguments
    root (1,1) string
    workers (1,1) double {mustBeInteger,mustBePositive} = 12
    authorize (1,1) logical = false
    outputRoot (1,1) string = ""
    dryRun (1,1) logical = false
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));paperPath=fullfile(root,"scripts","ieeg","paper");addpath(paperPath);c=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_4_panel_contract.json")));productionRoot=canonical(fullfile(root,string(c.outputs.root)));
assert(workers==double(c.execution.parallel_workers),"P0-4 frozen worker count is 12.");
if strlength(outputRoot)==0
    if dryRun,outputRoot=fullfile(productionRoot,"main_dryrun");else,outputRoot=productionRoot;end
end
outputRoot=canonical(outputRoot);
if dryRun
    assert(outputRoot~=productionRoot,"P0-4 dry-run must use an alternate output root.");
    seedIndices=1:2;status="dry_run";
else
    assert(authorize&&logical(c.execution.main_run_authorized),"P0-4 480-world main execution is not authorized.");
    assert(outputRoot==productionRoot,"Formal P0-4 run must target the production root.");
    seedIndices=1:double(c.panel.null_worlds);status="complete";
end
members=string(c.panel.opaque_members(:));
started=utcNow();maxNumCompThreads(1);pool=gcp("nocreate");if isempty(pool),pool=parpool("local",workers);elseif pool.NumWorkers~=workers,delete(pool);pool=parpool("local",workers);end
workerSetup=sprintf("addpath('%s'); addpath('%s'); maxNumCompThreads(1);",char(fullfile(root,"scripts","ieeg","acc01")),char(paperPath));pctRunOnAll(workerSetup);
spmd
    workerThreadCount=maxNumCompThreads;
    workerLabIndex=spmdIndex;
end
workerThreads=zeros(pool.NumWorkers,1);workerLabs=zeros(pool.NumWorkers,1);
for k=1:pool.NumWorkers,workerThreads(k)=workerThreadCount{k};workerLabs(k)=workerLabIndex{k};end
assert(isequal(sort(workerLabs(:))',(1:pool.NumWorkers))&&all(workerThreads==1),"Every P0-4 worker must be explicitly single-threaded.");
for mi=1:numel(members) % contract execution.worker_rationale: run members sequentially; do not co-schedule member classes
    memberSeeds=double(seedIndices(:));
    parfor si=1:numel(memberSeeds)
        run_acc00_p0_4_panel_world(root,members(mi),memberSeeds(si),"synthetic",outputRoot); %#ok<PFBNS> % member broadcast; seed is the parallel axis
    end
end
[counts,hashes]=checkpointManifest(outputRoot,members,seedIndices);assert(all([counts.complete]==numel(seedIndices)),"P0-4 orchestration did not complete its requested plan.");
if ~isfolder(outputRoot),mkdir(outputRoot),end
manifest=struct("status",status,"analysis_id",string(c.analysis_id),"contract_sha256",sha256File(fullfile(paperPath,"acc00_sim_p0_4_panel_contract.json")),"parent_p0_3_contract_sha256",string(c.parent_p0_3.sha256),"parallel_workers",workers,"null_worlds",double(c.panel.null_worlds),"opaque_members",members,"per_member_counts",counts,"total_checkpoints",4*double(c.panel.null_worlds),"completed_checkpoints",numel(hashes),"checkpoint_sha256",hashes,"started_at_utc",started,"completed_at_utc",utcNow(),"worker_max_num_threads_evidence",struct("requested",1,"labindex",workerLabs,"observed",workerThreads));
if ~dryRun
    assert(numel(hashes)==480&&all([counts.complete]==120),"Formal P0-4 manifest requires exactly 480 checkpoints.");
end
tmp=fullfile(outputRoot,"run_manifest.json.tmp");final=fullfile(outputRoot,"run_manifest.json");fid=fopen(tmp,"w");assert(fid>=0,"Cannot write P0-4 manifest tmp.");cleanup=onCleanup(@()fclose(fid));fprintf(fid,"%s\n",jsonencode(manifest,"PrettyPrint",true));clear cleanup;movefile(tmp,final,"f");
end
function [counts,hashes]=checkpointManifest(outputRoot,members,seedIndices)
counts=repmat(struct("opaque_id","","expected",0,"complete",0),numel(members),1);hashes=repmat(struct("opaque_id","","seed_index",0,"world_seed",0,"filename","","sha256",""),0,1);
for i=1:numel(members)
    id=members(i);files=dir(fullfile(outputRoot,"checkpoints",id,"world_*.mat"));files=sortByPath(files);counts(i)=struct("opaque_id",id,"expected",numel(seedIndices),"complete",numel(files));
    for j=1:numel(files)
        tok=regexp(files(j).name,"^world_(\d+)\.mat$","tokens","once");assert(~isempty(tok),"Unexpected P0-4 checkpoint filename.");si=str2double(tok{1});assert(ismember(si,seedIndices),"Unexpected checkpoint seed in requested output root.");path=fullfile(files(j).folder,files(j).name);q=load(path,"checkpoint");expected=91000000+100000*5+si;assert(string(q.checkpoint.opaque_id)==id&&q.checkpoint.seed_index==si&&q.checkpoint.world_seed==expected,"Checkpoint CRN identity mismatch.");hashes(end+1)=struct("opaque_id",id,"seed_index",si,"world_seed",q.checkpoint.world_seed,"filename",string(path),"sha256",sha256File(path)); %#ok<AGROW>
    end
end
end
function p=canonical(p),p=string(char(java.io.File(char(p)).getCanonicalPath()));end
function f=sortByPath(f),if isempty(f),return,end;p=string(arrayfun(@(x)fullfile(x.folder,x.name),f,"UniformOutput",false));[~,ix]=sort(p);f=f(ix);end
function s=utcNow(),s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cleanup=onCleanup(@()fclose(fid));b=fread(fid,Inf,"*uint8");md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(b,"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
