function out = summarize_acc00_p0_6(root,requireComplete)
%SUMMARIZE_ACC00_P0_6 Paired FPR contrast, verdict, corner-probe table.
% MATLAB-side aggregation by design: historical checkpoints are v7.3 tables
% whose string columns h5py cannot decode directly, so the paired contrast
% is computed here; the independent Python verifier cross-checks the flat
% TSVs dumped by this script.
% Reference pool: nullpool_edgeguard_20260829 (M-REF 1.5/1.5/0.75 config),
% NOT the older nullpool; see contract crn_and_seeding.matched_seed_domain.
arguments
    root (1,1) string = pwd
    requireComplete (1,1) logical = true
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));paperPath=fullfile(root,"scripts","ieeg","paper");
c=jsondecode(fileread(fullfile(paperPath,"acc00_sim_p0_6_perturbation_contract.json")));
od=fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation");
armDir=fullfile(od,"checkpoints");
m15=loadArm(armDir,"Mphase_1p5");m05=loadArm(armDir,"Mphase_0p5");p15=loadArm(armDir,"Mpower_1p5");cp=loadArm(armDir,"corner_probe_G06");
if requireComplete
    assert(numel(m15)==120&&numel(m05)==120&&numel(p15)==120&&numel(cp)==50,"P0-6 arm checkpoint counts must be 120/120/120/50.");
end
pparent=char(string(c.density_levels.("x1p0x").schedule_tsv));projRoot=pparent;
for k=1:6,projRoot=fileparts(projRoot);end
histDir=fullfile(projRoot,"processed","subject","group","ieeg_acc00_sim","BangYoureDead","checkpoints","nullpool_edgeguard_20260829");
h=loadArm(histDir,"");assert(numel(h)>=120,"Reference pool has %d worlds; 120 matched seeds required.",numel(h));
t15=worldTable("Mphase_1p5",m15);t05=worldTable("Mphase_0p5",m05);tp15=worldTable("Mpower_1p5",p15);tcp=worldTable("corner_probe_G06",cp);th=worldTableHist("historical_1p0",h);
writetable([t15;t05;tp15;tcp;th],fullfile(od,"per_world_metrics.tsv"),"FileType","text","Delimiter",char(9));
% primary contrast: Mphase 1.5x vs 1x (verdict); descriptive: 0.5x and Mpower
rowsPrimary=pairedContrast(t15,th,"Mphase","Mphase_reject","Mphase_reject","Mphase","1p5x");
rowsDesc=[pairedContrast(t05,th,"Mphase","Mphase_reject","Mphase_reject","Mphase","0p5x");pairedContrast(tp15,th,"Mpower","Mpower_reject","Mpower_reject","Mpower","1p5x")];
writetable([rowsPrimary;rowsDesc],fullfile(od,"paired_fpr_contrast.tsv"),"FileType","text","Delimiter",char(9));
verdictItems=repmat(struct("run","","label","","text","","Delta",0,"CI_lower",0,"CI_upper",0,"n",0),0,1);
for i=1:height(rowsPrimary)
    z=rowsPrimary(i,:);
    if z.CI_upper<=-0.05
        label="attenuation";text="Attenuation: the paired 95% Newcombe confidence interval for FPR at the high-density schedule minus FPR at the original schedule has an upper bound no greater than -0.05.";
    elseif z.CI_lower>=-0.02
        label="persistence";text="Persistence: the paired 95% Newcombe confidence interval for FPR at the high-density schedule minus FPR at the original schedule has a lower bound no less than -0.02.";
    else
        label="indeterminate";text="The density perturbation did not yield a decisive attenuation or persistence verdict under the preregistered precision rule.";
    end
    verdictItems(end+1)=struct("run",char(z.run),"label",char(label),"text",char(text),"Delta",z.Delta,"CI_lower",z.CI_lower,"CI_upper",z.CI_upper,"n",z.n); %#ok<AGROW>
end
verdict=struct("analysis_id",string(c.analysis_id),"interval_rule",string(c.verdict_rules.interval),"items",{verdictItems},"status","complete");
fid=fopen(fullfile(od,"primary_verdict.json"),"w");fprintf(fid,"%s\n",jsonencode(verdict,"PrettyPrint",true));fclose(fid);
crows=cornerTable(tcp,root);
writetable(crows,fullfile(od,"corner_probe_G06.tsv"),"FileType","text","Delimiter",char(9));
out=struct("world_rows",height([t15;t05;tp15;tcp;th]),"paired_contrast",[rowsPrimary;rowsDesc],"verdict",{verdictItems},"corner",crows);
end
function t=worldTable(name,cells)
R1R2=string(["R1";"R2"]);n=numel(cells)*2;
armArr=repmat(string(name),n,1);runArr=strings(n,1);si=zeros(n,1);mp=nan(n,1);mw=nan(n,1);fr=nan(n,1);fp=nan(n,1);
for i=1:numel(cells)
    q=cells(i).checkpoint;sf=q.signflip;z=q.subject_metrics;
    for j=1:2
        r=R1R2(j);k=(i-1)*2+j;
        runArr(k)=r;si(k)=q.seed_index;
        keep=z.run==r;
        mp(k)=mean(z.Mphase(keep),"omitnan");mw(k)=mean(z.Mpower(keep),"omitnan");
        fp(k)=sf.(char(r)).p_two_sided(1)<=0.05;% col1 = Mpower route
        fr(k)=sf.(char(r)).p_two_sided(2)<=0.05;% col2 = Mphase route
    end
end
t=table(armArr,runArr,si,mw,mp,fp,fr,'VariableNames',["arm" "run" "seed_index" "Mpower_subject_mean" "Mphase_subject_mean" "Mpower_reject" "Mphase_reject"]);
end
function T=pairedContrast(tNew,tOld,metricCol,colNew,colOld,fam,density)
R1R2=string(["R1";"R2"]);nRuns=2;
runArr=strings(nRuns,1);famArr=repmat(string(fam),nRuns,1);denArr=repmat(string(density),nRuns,1);n=zeros(nRuns,1);b=zeros(nRuns,1);cc=zeros(nRuns,1);
f1=nan(nRuns,1);f0=nan(nRuns,1);d=nan(nRuns,1);lo=nan(nRuns,1);hi=nan(nRuns,1);cp1lo=nan(nRuns,1);cp1hi=nan(nRuns,1);cp0lo=nan(nRuns,1);cp0hi=nan(nRuns,1);
for j=1:nRuns
    r=R1R2(j);
    a=tNew(tNew.run==r,:);t2=tOld(tOld.run==r,:);
    [ia,ib]=intersect(a.seed_index,t2.seed_index,"stable");
    x=a.(colNew)(ia);y=t2.(colOld)(ib);
    assert(numel(x)==height(a)&&numel(y)==numel(x),"P0-6 %s %s: matched pair count %d != expected %d.",metricCol,r,numel(x),height(a));
    nb=sum(x&~y);nc=sum(y&~x);nW=numel(x);eW=sum(x&y);
    [delta,l,u]=newcombePaired(nb,nc,eW,nW);
    runArr(j)=r;n(j)=nW;b(j)=nb;cc(j)=nc;f1(j)=mean(x);f0(j)=mean(y);d(j)=delta;lo(j)=l;hi(j)=u;
    [cp1lo(j),cp1hi(j)]=clopperPearson(nW,sum(x));[cp0lo(j),cp0hi(j)]=clopperPearson(nW,sum(y));
end
T=table(runArr,famArr,denArr,n,b,cc,f1,f0,d,lo,hi,cp1lo,cp1hi,cp0lo,cp0hi,'VariableNames',["run" "family" "density" "n" "b" "c" "FPR_1p5" "FPR_1p0" "Delta" "CI_lower" "CI_upper" "FPR_1p5_CP_low" "FPR_1p5_CP_high" "FPR_1p0_CP_low" "FPR_1p0_CP_high"]);
end
function [d,lo,hi]=newcombePaired(b,c,e,n)
% Newcombe (1998b) Stat Med 17:2635-2650, method 10 (paired difference):
% d=(b-c)/n; CI = d +- [sqrt(dl2^2 - 2*phi*dl2*du3 + du3^2), sqrt(du2^2 - 2*phi*du2*dl3 + dl3^2)]
% where (dl2,du2) are the Wilson deviations of p_b=(e+b)/n, (dl3,du3) of p_c=(e+c)/n,
% and phi = (e*h - b*c)/sqrt((e+b)(c+h)(e+c)(b+h)), h=n-b-c-e, with continuity
% correction on the numerator when e*h > b*c (max(e*h-b*c-n/2,0)); phi=0 if
% the denominator vanishes. In the phi=0 limit this is Newcombe (1998a)
% method 10 for independent proportions.
z=1.959963984540054;h=n-b-c-e;
[lb,ub]=wilson(e+b,n,z);[lc,uc]=wilson(e+c,n,z);
p1=(e+b)/n;p2=(e+c)/n;dl2=p1-lb;du2=ub-p1;dl3=p2-lc;du3=uc-p2;
den=sqrt((e+b)*(c+h)*(e+c)*(b+h));
num=e*h-b*c;
if den>0
    if num>0,num=min(num,max(num-n/2,0));end% method 10 continuity correction, applied only when e*h>f*g
    phi=num/den;
else
    phi=0;
end
d=(b-c)/n;lo=d-sqrt(dl2^2-2*phi*dl2*du3+du3^2);hi=d+sqrt(du2^2-2*phi*du2*dl3+dl3^2);
end
function [lo,hi]=wilson(k,n,z),p=k/n;den=1+z^2/n;cen=p+z^2/(2*n);rad=z*sqrt(p*(1-p)/n+z^2/(4*n^2));lo=(cen-rad)/den;hi=(cen+rad)/den;end
function [lo,hi]=clopperPearson(n,k)
% Clopper-Pearson exact binomial interval (plan 8.5: report per-level CP intervals).
a=0.025;
if k==0,lo=0;else,lo=betainv(a,k,n-k+1);end
if k==n,hi=1;else,hi=betainv(1-a,k+1,n-k);end
end
function t=cornerTable(cp,root)
% Descriptive contrast vs the P0-2 historical G06 cell (265 worlds).
R1R2=string(["R1";"R2"]);runs=strings(2,1);tot=zeros(2,1);cnt=zeros(2,1);rr=nan(2,1);htot=zeros(2,1);hcnt=zeros(2,1);hrr=nan(2,1);
histFile=fullfile(root,"processed","subject","group","ieeg_methods_paper","p0_2_sim_powercurve","tables","p0_2_powercurve_worlds.tsv");
assert(isfile(histFile),"P0-2 powercurve worlds table not found: %s",histFile);
w=readtable(histFile,"FileType","text","Delimiter",char(9),"TextType","string");
g=w(w.grid_id=="G06",:);
for j=1:2
    r=R1R2(j);z=cp(cp.run==r,:);
    runs(j)=r;tot(j)=height(z);cnt(j)=sum(z.Mpower_reject);rr(j)=mean(z.Mpower_reject);
    hg=g(g.run==r,:);htot(j)=height(hg);hcnt(j)=sum(hg.reject_alpha_0_05==1);hrr(j)=mean(hg.reject_alpha_0_05==1);
end
t=table(runs,tot,cnt,rr,htot,hcnt,hrr,'VariableNames',["run" "worlds_new" "rejections_new" "rejection_rate_new" "worlds_hist_G06" "rejections_hist_G06" "rejection_rate_hist_G06"]);
end
function cells=loadArm(armDir,arm)
d=dir(fullfile(armDir,arm,"world_*.mat"));
if strlength(arm)==0,d=dir(fullfile(armDir,"world_*.mat"));end
[~,ix]=sort({d.name});d=d(ix);
cells=repmat(struct("checkpoint",[]),numel(d),1);
for i=1:numel(d),q=load(fullfile(d(i).folder,d(i).name),"checkpoint");cells(i).checkpoint=q.checkpoint;end
end
function t=worldTableHist(name,cells)
R1R2=string(["R1";"R2"]);n=numel(cells)*2;
armArr=repmat(string(name),n,1);runArr=strings(n,1);si=zeros(n,1);mp=nan(n,1);mw=nan(n,1);fr=nan(n,1);fp=nan(n,1);
for i=1:numel(cells)
    q=cells(i).checkpoint;
    for j=1:2
        r=R1R2(j);k=(i-1)*2+j;
        runArr(k)=r;si(k)=q.seed_index;
        mp(k)=NaN;mw(k)=NaN;% historical table columns are not decoded (v7.3); route p is authoritative
        fp(k)=q.inference.(char(r)).p_two_sided(1)<=0.05;% col1 = Mpower route
        fr(k)=q.inference.(char(r)).p_two_sided(2)<=0.05;% col2 = Mphase route
    end
end
t=table(armArr,runArr,si,mw,mp,fp,fr,'VariableNames',["arm" "run" "seed_index" "Mpower_subject_mean" "Mphase_subject_mean" "Mpower_reject" "Mphase_reject"]);
end
