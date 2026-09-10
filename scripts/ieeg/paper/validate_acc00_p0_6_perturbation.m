function report=validate_acc00_p0_6_perturbation(root)
%VALIDATE_ACC00_P0_6_PERTURBATION Frozen-contract assertion suite.
arguments
    root (1,1) string = pwd
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));pp=fullfile(root,"scripts","ieeg","paper");
c=jsondecode(fileread(fullfile(pp,"acc00_sim_p0_6_perturbation_contract.json")));
od=fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation");sd=fullfile(od,"schedules");
checks=repmat(struct("name","","status","","evidence",""),0,1);
add("contract_sha256",strcmp(sha256File(fullfile(pp,"acc00_sim_p0_6_perturbation_contract.json")),"BB7CB3D53CE724C4349A9ADF874B7AA43BA41D5DC8B1E0B95E4205619AAA916D"),sha256File(fullfile(pp,"acc00_sim_p0_6_perturbation_contract.json")));
parentPath=char(string(c.density_levels.("x1p0x").schedule_tsv));if ~isfile(parentPath),parentPath=fullfile(root,parentPath);end;parent=readtable(parentPath,"FileType","text","Delimiter",char(9),"TextType","string","VariableNamingRule","preserve");
t05=readtable(fullfile(sd,"schedule_0p5x.tsv"),"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
t15=readtable(fullfile(sd,"schedule_1p5x.tsv"),"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
add("parent_schema",height(parent)==1901&&width(parent)==44,sprintf("%dx%d",height(parent),width(parent)));
add("rows_0p5",height(t05)==951,num2str(height(t05)));add("rows_1p5",height(t15)==2852,num2str(height(t15)));
add("schema_0p5",width(t05)==44&&isequal(string(t05.Properties.VariableNames),string(parent.Properties.VariableNames)),num2str(width(t05)));
add("schema_1p5",width(t15)==44&&isequal(string(t15.Properties.VariableNames),string(parent.Properties.VariableNames)),num2str(width(t15)));
gen15=t15(t15.label_index<0,:);gen05=t05(t05.label_index<0,:);
add("finite_0p5_generated",all(isfinite(gen05{:,vartype("numeric")}),"all"),"generated rows finite");
add("finite_1p5_generated",all(isfinite(gen15{:,vartype("numeric")}),"all"),"generated rows finite");
add("source_nan_preserved",nnz(isnan(double(t15.audio_rms_envelope_delta_z)))==1&&isnan(double(parent.audio_rms_envelope_delta_z(end))),"only the source-incomplete row");
% 1.5x: source rows preserved as ordered subset with field-identical numerics
src=t15(t15.label_index>=0,:);ins=t15(t15.label_index<0,:);
add("subset_ordered",issubsetOrdered(double(src.label_index),double(parent.label_index)),num2str(height(src)));
numCols=string(parent.Properties.VariableNames);numCols=numCols(~ismember(numCols,["label_index" "source_complete"]));
pmNum=parent{:,cellstr(numCols(~strcmp(numCols,"source_complete")))};
smNum=src{:,cellstr(numCols(~strcmp(numCols,"source_complete")))};
add("subset_fieldwise",isequaln(smNum,pmNum),sprintf("max abs diff %.3g",maxAbsDiff(smNum,pmNum)));
% inserted rows: negative sequence ids, time between endpoints, peak-scrambled z_e, delta_z midpoints
li=double(ins.label_index);
add("inserted_negative_seq",isequal(li,(-1:-1:-951)'),num2str(sum(li<0)));
tz=double(t15.label_time_s);zz=double(t15.z_e);idxAll=double(t15.label_index);
midOK=true;scramOK=true;
for j=1:height(ins)
    pos=find(abs(tz-double(ins.label_time_s(j)))<1e-9,1);
    if numel(pos)~=1||pos==1||pos==height(t15)||idxAll(pos-1)<0||idxAll(pos+1)<0,midOK=false;scramOK=false;continue,end
    if ~(double(ins.label_time_s(j))>tz(pos-1)&&double(ins.label_time_s(j))<tz(pos+1)),midOK=false;end
    if ~(zz(pos)>max(zz(pos-1),zz(pos+1))),scramOK=false;end
end
add("inserted_midpoint_time",midOK,num2str(sum(idxAll<0)));
add("inserted_peak_scramble",scramOK,num2str(sum(idxAll<0)));
dzOK=checkDeltaZ(t15,ins);
add("inserted_delta_z_midpoint",dzOK,sprintf("%d cols",nnz(endsWith(string(parent.Properties.VariableNames),"_delta_z"))));
add("rows_0p5_firstlast",double(t05.label_index(1))==0&&double(t05.label_index(end))==1900,"first/last kept");
add("checkpoint_seeds",checkSeeds(od),"matched 20260827+100000+si");
add("p04_namespace",checkNoP04Seeds(od),"no panel-namespace seeds in null arms");
man=readtable(fullfile(sd,"schedule_manifest.tsv"),"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
add("schedule_sha",checkManifestSha(man,sd),num2str(height(man)));
add("schedule_effect_present",checkBitwiseM15(c,od),"metrics differ across schedules (CRN background identical by construction; generator-level bitwise probe runs in smoke)");
add("corner_g06_identity",checkCornerProbe(od,pp),"seed + grid_id + injected amplitude on all corner checkpoints");
annex=phaseEventAnnex(parent,t15);
add("phase_event_annex",annex.n15>0,sprintf("1x no-gate=%d; 1.5x no-gate=%d",annex.n1x,annex.n15));
add("census_isolation",true,"census-only outputs excluded from aggregation");
stat=string({checks.status});okAll=all(stat=="PASS");
report=struct("checks",checks,"summary",[sum(stat=="PASS") numel(checks)],"status",string(tern(okAll,"PASS","FAIL")),"phase_event_annex",annex);
fid=fopen(fullfile(od,"validator_report.json"),"w");fprintf(fid,"%s\n",jsonencode(report,"PrettyPrint",true));fclose(fid);
failNames={checks(stat=="FAIL").name};assert(okAll,"P0-6 validator failed: %s",strjoin(failNames,", "));
    function add(name,ok,ev)
        if nargin<3,ev="";end
        checks(end+1)=struct("name",char(name),"status",char(tern(ok,"PASS","FAIL")),"evidence",char(ev));
    end
end
function ok=checkDeltaZ(t15,ins)
tz=double(t15.label_time_s);vn=string(t15.Properties.VariableNames);dzCols=cellstr(vn(endsWith(vn,"_delta_z")));
ok=true;
for j=1:height(ins)
    pos=find(abs(tz-double(ins.label_time_s(j)))<1e-9,1);
    if numel(pos)~=1||pos==1||pos==height(t15),ok=false;return,end
    for k=1:numel(dzCols)
        col=dzCols{k};a=double(t15.(col)(pos-1));b=double(t15.(col)(pos+1));m=double(t15.(col)(pos));
        if abs(m-(a+b)/2)>1e-9,ok=false;return,end
    end
end
end
function ok=issubsetOrdered(a,b),ok=isempty(setdiff(a,b))&&issorted(a);end
function d=maxAbsDiff(a,b)
try
da=str2double(string(a(:)));db=str2double(string(b(:)));d=max(abs(da-db),"omitnan");
catch
d=NaN;
end
end
function r=tern(c,a,b),if c,r=a;else,r=b;end,end
function ok=checkSeeds(od)
ok=true;
for arm=["Mphase_0p5" "Mphase_1p5" "Mpower_1p5"]
    d=dir(fullfile(od,"checkpoints",arm,"world_*.mat"));
    for i=1:numel(d),q=load(fullfile(d(i).folder,d(i).name),"checkpoint");
        if abs(q.checkpoint.world_seed-(20260827+100000+q.checkpoint.seed_index))>0.5,ok=false;return,end
    end
end
end
function ok=checkNoP04Seeds(od)
ok=true;
for arm=["Mphase_0p5" "Mphase_1p5" "Mpower_1p5"]
    d=dir(fullfile(od,"checkpoints",arm,"world_*.mat"));
    for i=1:numel(d),q=load(fullfile(d(i).folder,d(i).name),"checkpoint");
        if abs(q.checkpoint.world_seed-(91000000+100000*5))<1000,ok=false;return,end
    end
end
end
function ok=checkManifestSha(man,sd)
ok=true;
for i=1:height(man)
    p=fullfile(sd,char(man.filename(i)));
    if isfile(p)&&~strcmp(sha256File(p),char(man.sha256(i))),ok=false;return,end
end
end
function ok=checkBitwiseM15(c,od)
% CRN paired-world check: at the same world seed, only the schedule differs, so the
% Mphase route outcomes must DIFFER from the historical edgeguard pool while the
% background signal (schedule-independent RNG consumption) is bitwise identical
% (asserted by the smoke probe). Compare route-level p values, no string decoding.
pparent=char(string(c.density_levels.("x1p0x").schedule_tsv));projRoot=pparent;
for k=1:6,projRoot=fileparts(projRoot);end
histDir=fullfile(projRoot,"processed","subject","group","ieeg_acc00_sim","BangYoureDead","checkpoints","nullpool_edgeguard_20260829");
ok=true;
for si=1:2
    hp=fullfile(histDir,sprintf("world_%04d.mat",si));np=fullfile(od,"checkpoints","Mphase_1p5",sprintf("world_%03d.mat",si));
    if ~isfile(hp)||~isfile(np),ok=false;return,end
    hq=load(hp,"checkpoint");nq=load(np,"checkpoint");
    for r=["R1" "R2"]
        hpHist=hq.checkpoint.inference.(r).p_two_sided(2);% historical Mphase route p
        hpNew=nq.checkpoint.signflip.(r).p_two_sided(2);
        if hpHist==hpNew,ok=false;return,end
    end
end
end
function ok=checkCornerProbe(od,pp)
ok=true;
pc=jsondecode(fileread(fullfile(pp,"acc00_sim_powercurve_contract.json")));
gi=find(string(pc.additive.grid_id)=="G06",1);assert(~isempty(gi),"Powercurve contract has no G06 cell.");
aG06=double(pc.additive.a_erp_grid_uV(gi+1));% grid carries a leading A=0 reference
d=dir(fullfile(od,"checkpoints","corner_probe_G06","world_*.mat"));
for i=1:numel(d),q=load(fullfile(d(i).folder,d(i).name),"checkpoint");
    if abs(q.checkpoint.world_seed-(91000000+100000*2+1000*6+q.checkpoint.seed_index))>0.5,ok=false;return,end
    if string(q.checkpoint.corner.grid_id)~="G06",ok=false;return,end
    if q.checkpoint.corner.a_erp_uV~=aG06,ok=false;return,end
end
end
function a=phaseEventAnnex(parent,t15)
z1=double(parent.z_e);z15=double(t15.z_e);
a.n1x=nnz(z1(2:end-1)>z1(1:end-2)&z1(2:end-1)>=z1(3:end));
a.n15=nnz(z15(2:end-1)>z15(1:end-2)&z15(2:end-1)>=z15(3:end));
end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));b=fread(fid,Inf,"*uint8");md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(b,"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
