function manifest = run_acc00_p0_3_ablation_main(root, workers, authorize)
%RUN_ACC00_P0_3_ABLATION_MAIN Frozen 768-world orchestration, deliberately gated.
arguments
    root (1,1) string
    workers (1,1) double {mustBeInteger,mustBePositive} = 1
    authorize (1,1) logical = false
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));paperPath=fullfile(root,"scripts","ieeg","paper");paper=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_3_ablation_contract.json")));
assert(authorize&&logical(paper.execution.main_run_authorized),"P0-3 full 768-world execution is not authorized by the immutable contract.");
assert(workers==double(paper.execution.parallel_workers),"Frozen worker count is 12.");mpowerPlan=buildPlan(paper.cells.Mpower);mphasePlan=buildPlan(paper.cells.Mphase);assert(size(mpowerPlan,1)==640&&size(mphasePlan,1)==128,"Frozen plan must be 640+128 worlds.");pool=gcp("nocreate");if isempty(pool),parpool(workers);elseif pool.NumWorkers~=workers,delete(pool);parpool(workers);end
% Arms are intentionally serial; each arm uses the same 12-worker pool.
parfor k=1:size(mpowerPlan,1),run_acc00_p0_3_mpower_world(root,mpowerPlan(k,1),str2double(mpowerPlan(k,2)),"synthetic");end
parfor k=1:size(mphasePlan,1),run_acc00_p0_3_mphase_world(root,mphasePlan(k,1),str2double(mphasePlan(k,2)),"synthetic");end
out=fullfile(root,string(paper.outputs.root));if ~isfolder(out),mkdir(out),end;[counts,hashes]=checkpointManifest(out,paper);assert(all([counts.complete]==[paper.cells.Mpower.worlds paper.cells.Mphase.worlds]),"Incomplete P0-3 cells before manifest.");manifest=struct("status","complete","analysis_id",paper.analysis_id,"worlds",768,"workers",workers,"contract_sha256",sha256File(fullfile(paperPath,"acc00_sim_p0_3_ablation_contract.json")),"per_cell_completion",counts,"checkpoint_sha256",hashes,"completed_at_utc",string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")));fid=fopen(fullfile(out,"run_manifest.json"),"w");fprintf(fid,"%s\n",jsonencode(manifest,"PrettyPrint",true));fclose(fid);
end
function plan=buildPlan(cells),plan=strings(0,2);for ci=1:numel(cells),c=cells(ci);for si=1:double(c.worlds),plan(end+1,:)=[string(c.id),string(si)];end;end,end %#ok<AGROW> % jsondecode struct arrays are Nx1; "for c=cells" would iterate columns (once), so index explicitly
function [counts,hashes]=checkpointManifest(out,paper),allCells=[num2cell(paper.cells.Mpower);num2cell(paper.cells.Mphase)];counts=repmat(struct("cell_id","","expected",0,"complete",0),numel(allCells),1);hashes=repmat(struct("cell_id","","filename","","sha256",""),0,1);for i=1:numel(allCells),c=allCells{i};arm="mpower";if startsWith(string(c.id),"M"),arm="mphase";end;files=dir(fullfile(out,arm,"checkpoints",string(c.id),"world_*.mat"));counts(i)=struct("cell_id",string(c.id),"expected",double(c.worlds),"complete",numel(files));for j=1:numel(files),hashes(end+1)=struct("cell_id",string(c.id),"filename",string(fullfile(files(j).folder,files(j).name)),"sha256",sha256File(fullfile(files(j).folder,files(j).name)));end,end,end %#ok<AGROW> % num2cell of jsondecode column structs yields columns; vertical concat
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));b=fread(fid,Inf,"*uint8");md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(b,"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
