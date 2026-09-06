function manifest = run_acc00_powercurve_main(root, workers, authorize)
% Future 986-world Mpower-only runner. Main execution is contract-gated.
arguments
    root (1,1) string
    workers (1,1) double {mustBeInteger,mustBePositive} = 1
    authorize (1,1) logical = false
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
paper=jsondecode(fileread(fullfile(root,"scripts","ieeg","paper","acc00_sim_powercurve_contract.json")));
assert(authorize && logical(paper.execution.main_run_authorized),"Main run is not authorized by paper-only contract; no world started.");
assert(workers>0,"workers must be positive");
plan=[]; for gi=1:numel(paper.additive.grid_worlds), for si=1:double(paper.additive.grid_worlds(gi)), plan=[plan; gi si]; end, end %#ok<AGROW>
outDir=fullfile(root,"processed","subject","group","ieeg_methods_paper","p0_2_sim_powercurve"); if ~isfolder(outDir),mkdir(outDir);end
pool=gcp("nocreate");
if ~isempty(pool)
    if pool.NumWorkers~=workers, delete(pool); parpool(workers); end
else
    parpool(workers);
end
parfor k=1:size(plan,1), run_acc00_powercurve_world(root,plan(k,1),plan(k,2),"main"); end
for gi=1:numel(paper.additive.grid_worlds)
    cellDir=fullfile(outDir,"profiling","main"); files=dir(fullfile(cellDir,sprintf("world_g%02d_n*.mat",gi)));
    assert(numel(files)==double(paper.additive.grid_worlds(gi)),"Incomplete grid cell G%02d before aggregation",gi);
    for j=1:numel(files), q=load(fullfile(files(j).folder,files(j).name),"checkpoint"); assertIdentity(q.checkpoint,paper,root,gi,q.checkpoint.seed_index,q.checkpoint.seed); end
end
manifest=struct("status","complete","workers",workers,"worlds",size(plan,1),"contract_sha256",upper(string(sha256File(fullfile(root,"scripts","ieeg","paper","acc00_sim_powercurve_contract.json")))),"completed_at_utc",string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z")));
fid=fopen(fullfile(outDir,"run_manifest.json"),"w"); fprintf(fid,"%s\n",jsonencode(manifest,"PrettyPrint",true)); fclose(fid);
end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance("SHA-256");while true,q=fread(fid,1048576,"*uint8");if isempty(q),break,end;md.update(typecast(q,"int8"));end;h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
function assertIdentity(q,paper,root,gi,si,seed), assert(q.grid_index==gi&&q.seed_index==si&&q.seed==seed&&string(q.contract_sha256)==sha256File(fullfile(root,"scripts","ieeg","paper","acc00_sim_powercurve_contract.json"))&&string(q.grid_id)==string(paper.additive.grid_id(gi)),"Main checkpoint identity mismatch"); end
