function manifest=run_acc00_p0_6_main(root,workers,runMode)
%RUN_ACC00_P0_6_MAIN P0-6 orchestration: serial arms, parfor over seeds.
% Contract execution.orchestration forbids co-scheduling arms; the seed axis
% is the only parallel one. "dry-run" is limited to seed_index 1:2 per arm.
arguments
    root (1,1) string = pwd
    workers (1,1) double = 12
    runMode (1,1) string {mustBeMember(runMode,["dry-run" "complete"])} = "dry-run"
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));paperPath=fullfile(root,"scripts","ieeg","paper");addpath(paperPath);
c=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_6_perturbation_contract.json")));
sdir=fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation","schedules");
assert(isfile(fullfile(sdir,"schedule_0p5x.tsv"))&&isfile(fullfile(sdir,"schedule_1p5x.tsv")),"Frozen schedules are missing; run generate_acc00_p0_6_schedules first.");
arms=["Mphase_0p5" "Mphase_1p5" "Mpower_1p5" "corner_probe_G06"];
started=utcNow();
pool=gcp("nocreate");if isempty(pool),parpool("local",workers);elseif pool.NumWorkers~=workers,delete(pool);parpool("local",workers);end
pctRunOnAll(sprintf("addpath('%s');addpath('%s');addpath('%s');maxNumCompThreads(1);",char(fullfile(root,"scripts","ieeg")),char(fullfile(root,"scripts","ieeg","acc01")),char(paperPath)));
for a=arms
    na=120;if a=="corner_probe_G06",na=50;end;if runMode=="dry-run",na=2;end
    todo=[];
    for si=1:na
        f=fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation","checkpoints",a,sprintf("world_%03d.mat",si));
        if isfile(f)
            q=load(f,"checkpoint");
            if string(q.checkpoint.status)=="complete"&&q.checkpoint.seed_index==si&&string(q.checkpoint.arm)==a&&string(q.checkpoint.contract_sha256)==sha256File(fullfile(paperPath,"acc00_sim_p0_6_perturbation_contract.json"))
                continue;
            end
            delete(f);% identity mismatch (e.g. contract re-freeze): discard and recompute
        end
        t=fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation","checkpoints",a,sprintf("world_%03d.mat.tmp.mat",si));
        if isfile(t),delete(t),end
        todo(end+1)=si; %#ok<AGROW>
    end
    parfor k=1:numel(todo)
        run_acc00_p0_6_world(root,a,todo(k));% arm broadcast; seed is the parallel axis; %#ok<PFBNS>
    end
    fprintf("P0-6 arm %s done (%d worlds).\n",a,na);
end
expectedPerArm=[120 120 120 50];if runMode=="dry-run",expectedPerArm(:)=2;end
[counts,hashes]=checkpointManifest(root,arms,expectedPerArm);
assert(all([counts.complete]==[counts.expected]),"P0-6 orchestration did not complete its requested plan.");
od=fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation");if ~isfolder(od),mkdir(od),end
manifest=struct("status",runMode,"analysis_id",string(c.analysis_id),"contract_sha256",sha256File(fullfile(paperPath,"acc00_sim_p0_6_perturbation_contract.json")),"parallel_workers",workers, ...
    "arms",arms,"per_arm_counts",counts,"total_checkpoints",sum([counts.expected]),"completed_checkpoints",numel(hashes),"checkpoint_sha256",hashes, ...
    "started_at_utc",started,"completed_at_utc",utcNow());
tmp=fullfile(od,"run_manifest.json.tmp");fid=fopen(tmp,"w");fprintf(fid,"%s\n",jsonencode(manifest,"PrettyPrint",true));fclose(fid);movefile(tmp,fullfile(od,"run_manifest.json"),"f");
end
function [counts,hashes]=checkpointManifest(root,arms,expected)
counts=repmat(struct("arm","","expected",0,"complete",0),numel(arms),1);hashes=repmat(struct("arm","","seed_index",0,"world_seed",0,"filename","","sha256",""),0,1);
for i=1:numel(arms)
    d=dir(fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation","checkpoints",arms(i),"world_*.mat"));d=sortFiles(d);
    counts(i)=struct("arm",arms(i),"expected",expected(i),"complete",numel(d));
    for j=1:numel(d)
        tok=regexp(d(j).name,"^world_(\d+)\.mat$","tokens","once");si=str2double(tok{1});
        q=load(fullfile(d(j).folder,d(j).name),"checkpoint");
        assert(string(q.checkpoint.arm)==arms(i)&&q.checkpoint.seed_index==si,"P0-6 checkpoint identity mismatch in %s.",d(j).name);
        hashes(end+1)=struct("arm",arms(i),"seed_index",si,"world_seed",q.checkpoint.world_seed,"filename",string(fullfile(d(j).folder,d(j).name)),"sha256",sha256File(fullfile(d(j).folder,d(j).name))); %#ok<AGROW>
    end
end
end
function d=sortFiles(d),[~,ix]=sort({d.name});d=d(ix);end
function s=utcNow(),s=string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));b=fread(fid,Inf,"*uint8");md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(b,"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
