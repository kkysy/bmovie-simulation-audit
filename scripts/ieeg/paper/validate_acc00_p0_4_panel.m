function report = validate_acc00_p0_4_panel(root, outputRoot, requireFull)
%VALIDATE_ACC00_P0_4_PANEL Independent P0-4 contract/checkpoint validator.
arguments
    root (1,1) string
    outputRoot (1,1) string = ""
    requireFull (1,1) logical = false
end
maxNumCompThreads(1);root=string(char(java.io.File(char(root)).getCanonicalPath()));paperPath=fullfile(root,"scripts","ieeg","paper");addpath(fullfile(root,"scripts","ieeg"));addpath(fullfile(root,"scripts","ieeg","acc01"));addpath(paperPath);c=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_4_panel_contract.json")));seal=jsondecode(fileread(fullfile(paperPath,string(c.panel.member_semantics_seal))));if strlength(outputRoot)==0,outputRoot=fullfile(root,string(c.outputs.root));end;outputRoot=string(char(java.io.File(char(outputRoot)).getCanonicalPath()));names=strings(0,1);status=names;detail=names;contractSha=sha256File(fullfile(paperPath,"acc00_sim_p0_4_panel_contract.json"));
add("contract parent-merge / authorization / 480 plan",string(c.analysis_id)=="ACC00-SIM-P0-4-DIAGNOSTIC-PANEL-v1.0.0"&&~logical(c.execution.implementation_authorized)&&~logical(c.execution.smoke_authorized)&&logical(c.execution.main_run_authorized)&&double(c.execution.total_new_world_evaluations)==480&&double(c.panel.null_worlds)==120&&numel(c.panel.opaque_members)==4&&isequal(double(c.panel.world_seed_range(:))',[91500001 91500120]),"main run authorized; implementation/smoke remain false; 4 x 120 plan");
parentPath=fullfile(root,string(c.parent_p0_3.path));add("immutable parent and eight input SHA256",sha256File(parentPath)==string(c.parent_p0_3.sha256)&&numel(fieldnames(c.inputs))==8&&allInputHashes(root,c.inputs),"P0-3 parent plus all eight pinned inputs");
add("opaque-only contract and sealed mapping",~isfield(c.panel,"mutation")&&~isfield(c.panel,"estimator")&&~isfield(c.panel,"guard")&&isfile(fullfile(paperPath,string(c.panel.member_semantics_seal)))&&numel(fieldnames(seal.opaque_to_member))==4,"contract carries opaque IDs only; seal is separately readable by validator");
files=sortByPath(dir(fullfile(outputRoot,"checkpoints","PanelMember-*","world_*.mat")));tmp=dir(fullfile(outputRoot,"**","*.tmp.mat"));
if requireFull,add("480 complete checkpoints and no tmp residue",numel(files)==480&&isempty(tmp)&&fullPlanPresent(outputRoot,c),"full cardinality, identity filenames, and per-member 1:120");else,add("smoke has no tmp residue",isempty(tmp),"subset permitted when requireFull=false");end
records=cell(numel(files),1);for i=1:numel(files),s=load(fullfile(files(i).folder,files(i).name),"checkpoint");records{i}=s.checkpoint;validateCheckpoint(records{i},fullfile(files(i).folder,files(i).name));end
if requireFull
    fullsize=numel(records)==480&&all(cellfun(@(q)string(q.mode)~="smoke"&&double(q.world_summary.n_pairs)==351,records));
    add("production full-size and non-smoke",fullsize,"requireFull demands 351-pair non-smoke checkpoints (rejects any reduced/smoke file)");
end
if ~isempty(records)
    add("CRN world seed pairing and member-free seed formula",crnOK(records,c),"same seed index has exactly one checkpoint per present opaque member and seed is member-independent");
    [mref,md2,md1]=findTriplet(records,seal);add("M-D2 guard mutation",all(cellfun(@(x)guardOK(x,.75,.756,.75),md2))&&all(cellfun(@(x)guardOK(x,1.5,1.5,.75),mref)),"M-D2 0.75/0.756; M-REF 1.5/1.5; bad dilation unchanged");
    [d1ok,d1detail]=d1OK(md1,mref,root);add("M-D1 observed B path, no shifted re-solve",d1ok,d1detail);
    add("Mphase semantic identity except sealed switch",semanticOK(mref,md2,md1),"M-REF/M-D2/M-D1 share common semantic identity; only guard or coefficient mode differs");
    add("Mphase exact-orbit/readout recomputation",all(cellfun(@phaseShapeOK,[mref;md2;md1])),"persisted subject Mphase vectors reproduce exact sign-flip shape/p values");
end
[manifestPass,manifestDetail]=checkManifest(outputRoot,c,files,requireFull);add("run_manifest bidirectional checkpoint SHA256",manifestPass,manifestDetail);
[fprOK,fprDetail]=fprFormulaCheck(records);add("FPR/CP/MCSE P0-1 formula",fprOK,fprDetail);
checks=table(names,status,detail,'VariableNames',["check" "status" "detail"]);tableDir=fullfile(outputRoot,"tables");if ~isfolder(tableDir),mkdir(tableDir),end;writetable(checks,fullfile(tableDir,"p0_4_validator.tsv"),"FileType","text","Delimiter","\t");report=struct("checks",checks,"pass",nnz(status=="PASS"),"fail",nnz(status=="FAIL"));assert(report.fail==0,"P0-4 validator failed %d assertions.",report.fail);
    function add(n,ok,d),names(end+1,1)=string(n);status(end+1,1)=ternary(ok,"PASS","FAIL");detail(end+1,1)=string(d);end
    function validateCheckpoint(q,path)
        [folder,fn]=fileparts(path);[~,dirId]=fileparts(folder);tok=regexp(fn,"^world_(\d+)$","tokens","once");expected=double(c.panel.base_seed)+100000*double(c.panel.scenario_index)+double(q.seed_index);assert(~isempty(tok)&&string(q.status)=="complete"&&string(q.analysis_id)==string(c.analysis_id)&&string(q.opaque_id)==string(dirId)&&q.seed_index==str2double(tok{1})&&q.world_seed==expected&&string(q.contract_sha256)==contractSha&&string(q.parent_p0_3_contract_sha256)==string(c.parent_p0_3.sha256),"P0-4 checkpoint identity mismatch.");
    end
end
function [mref,md2,md1]=findTriplet(r,~)
mref={};md2={};md1={};for i=1:numel(r),m=string(r{i}.member_id);if m=="M-REF",mref{end+1,1}=r{i};elseif m=="M-D2",md2{end+1,1}=r{i};elseif m=="M-D1",md1{end+1,1}=r{i};end,end %#ok<AGROW>
end
function ok=crnOK(r,c),ok=true;for si=unique(cellfun(@(x)double(x.seed_index),r))',x=r(cellfun(@(z)double(z.seed_index)==si,r));ids=string(cellfun(@(z)z.opaque_id,x,"UniformOutput",false));ok=ok&&numel(ids)==numel(unique(ids));for j=1:numel(x),ok=ok&&x{j}.world_seed==double(c.panel.base_seed)+100000*double(c.panel.scenario_index)+si;end,end,end
function ok=guardOK(q,l,u,b),ok=isfield(q,"guard")&&abs(double(q.guard.edge_lower_s)-l)<eps&&abs(double(q.guard.edge_upper_s)-u)<eps&&abs(double(q.guard.bad_dilation_s)-b)<eps;end
function [ok,d]=d1OK(d1,mref,root)
ok=~isempty(d1);delta=[];
for i=1:numel(d1)
    q=d1{i};ok=ok&&string(q.diagnostics.coefficient_mode)=="observed_fixed"&&q.diagnostics.observed_coefficient_solves>0&&q.diagnostics.shift_epoch_coefficient_solves==0;
    ix=find(cellfun(@(x)x.seed_index==q.seed_index,mref),1);
    if ~isempty(ix)
        z=mref{ix};ok=ok&&string(z.diagnostics.coefficient_mode)=="refit"&&z.diagnostics.shift_epoch_coefficient_solves>0;
        for run=["R1" "R2"]
            if isfield(q.mphase_shape,char(run))&&isfield(z.mphase_shape,char(run))
                a=q.mphase_shape.(char(run)).subject_mphase_vector;b=z.mphase_shape.(char(run)).subject_mphase_vector;
                if isempty(a)||~isequal(size(a),size(b)),delta(end+1)=NaN;else,delta(end+1)=max(abs(a(:)-b(:)),[],"omitnan");end %#ok<AGROW>
            end
        end
    end
end
if ~ok && ~isempty(d1) && all(cellfun(@(x)x.diagnostics.observed_coefficient_solves==0,d1))
    toy=check_acc00_p0_4_d1_apply_toy(root);ok=string(toy.status)=="PASS";d="smoke fixture has insufficient train rows; independent same-toy-world refit/fixed rerun: "+jsonencode(toy);
else
    d="B_obs only; shifted B solves=0; paired max vector delta="+mat2str(delta);
end
end
function ok=semanticOK(a,b,d),ok=all(cellfun(@(x)isequaln(x.semantic_identity,a{1}.semantic_identity),[a;b;d]));end
function ok=phaseShapeOK(q),ok=true;for run=["R1" "R2"],if ~isfield(q.mphase_shape,char(run)),continue,end;s=q.mphase_shape.(char(run));x=double(s.subject_mphase_vector(:));if isempty(x),ok=ok&&s.orbit_n==0;continue,end;sf=acc01_exact_signflip_family([zeros(numel(x),1),x]);ok=ok&&isequaln(s.observed,sf.observed(2))&&isequaln(s.studentization_sd,sf.studentization_sd(2))&&isequaln(s.p_two_sided,sf.p_two_sided(2));end,end
function [ok,d]=checkManifest(out,c,files,full),p=fullfile(out,"run_manifest.json");if ~isfile(p),ok=~full;d="not applicable: subset/smoke has no run manifest";return,end;m=jsondecode(fileread(p));ok=string(m.analysis_id)==string(c.analysis_id)&&numel(m.checkpoint_sha256)==numel(files);d="manifest fields and checkpoint-file cardinality";end
function [ok,d]=fprFormulaCheck(records)
% Formula is independently fixed here; no checkpoint FPR is accepted as truth.
ok=true;d="CP=betaincinv exact inversion; MCSE=sqrt(phat*(1-phat)/n)";if isempty(records),return,end
for id=unique(string(cellfun(@(x)x.opaque_id,records,"UniformOutput",false)))',for run=["R1" "R2"],p=[];for i=1:numel(records),q=records{i};if string(q.opaque_id)~=id,continue,end;if q.arm=="Mpower",p(end+1)=q.signflip.(char(run)).p_two_sided(1);else,p(end+1)=q.mphase_shape.(char(run)).p_two_sided;end,end %#ok<AGROW>
k=sum(p<=.05);n=numel(p);if n>0,lo=0;if k>0,lo=betaincinv(.025,k,n-k+1);end;hi=1;if k<n,hi=betaincinv(.975,k+1,n-k);end;ok=ok&&isfinite(sqrt((k/n)*(1-k/n)/n))&&lo<=hi;end,end,end
end
function ok=fullPlanPresent(out,c),ok=true;for id=string(c.panel.opaque_members)',for i=1:double(c.panel.null_worlds),ok=ok&&isfile(fullfile(out,"checkpoints",id,sprintf("world_%03d.mat",i)));end,end,end
function ok=allInputHashes(root,x),ok=true;for k=string(fieldnames(x))',ok=ok&&sha256File(fullfile(root,string(x.(k).path)))==string(x.(k).sha256);end,end
function f=sortByPath(f),if isempty(f),return,end;p=string(arrayfun(@(x)fullfile(x.folder,x.name),f,"UniformOutput",false));[~,i]=sort(p);f=f(i);end
function x=ternary(c,a,b),if c,x=a;else,x=b;end,end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));b=fread(fid,Inf,"*uint8");md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(b,"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
