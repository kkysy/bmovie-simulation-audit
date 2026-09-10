function world = stagec_generate_gainflat_world(root, contract, mode, seed)
%STAGEC_GENERATE_GAINFLAT_WORLD Stage-C gain-flattened P0F synthesis.
% Identical base synthesis to Stage B (same seeds => bitwise-paired
% backgrounds/nuisance); the only difference is the injection gain slope b=0
% (read from the Stage C contract) and a single arm P0F (kappa=0).
arguments
    root (1,1) string
    contract (1,1) struct
    mode (1,1) string {mustBeMember(mode,["smoke" "benchmark"])}
    seed (1,1) double {mustBeInteger,mustBePositive}
end
base=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
stream=RandStream("Threefry","Seed",seed);
covariates=readtable(inputPath(root,base.inputs.label_covariates.path), ...
    "FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
if mode=="smoke"
    covariates=covariates(1:12,:); covariates.label_time_s=linspace(2,13,12)';
    covariates.z_e=sin(2*pi*(0:11)'/3); covariates.source_complete(:)=true;
    p=table(repelem(["S01";"S02"],2),repelem(["S01_R1";"S02_R1"],2), ...
        repmat("R1",4,1),"pair"+(1:4)',repmat("A",4,1),repmat("B",4,1), ...
        repmat(20,4,1),'VariableNames',["subject" "session_id" "run" ...
        "pair_name" "contact_a" "contact_b" "T"]);
    p.bad_segments_s=repmat({zeros(0,2)},height(p),1);
else
    p0=readtable(inputPath(root,base.inputs.i00_acc_pairs.path), ...
        "FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
    p0=p0(logical(p0.qc_eligible)&contains(p0.location,"ACC"),:);
    assert(height(p0)==351,"Expected 351 ACC pairs.");
    timing=readtable(inputPath(root,base.inputs.i00_timing.path), ...
        "FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
    bad=readtable(inputPath(root,base.inputs.i00_bad_segments.path), ...
        "FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
    T=nan(height(p0),1); badCell=cell(height(p0),1);
    for i=1:height(p0)
        tr=timing(timing.session_id==p0.session_id(i),:); assert(height(tr)==1,"Timing join failed.");
        T(i)=double(tr.usable_movie_stop_s);
        br=bad(bad.session_id==p0.session_id(i)&bad.series=="LFP_macro"& ...
            ismember(bad.origchannel_name,[p0.contact_a(i),p0.contact_b(i)]),:);
        movieZero=double(tr.macro_common_start_s)+(double(tr.macro_movie_first_sample)-1)/double(tr.macro_rate_hz);
        iv=[double(br.start_common_time_s),double(br.stop_common_time_s)]-movieZero;
        iv(:,1)=max(iv(:,1),0);iv(:,2)=min(iv(:,2),T(i));badCell{i}=mergeIntervals(iv(iv(:,2)>iv(:,1),:));
    end
    p=table(p0.subject,p0.session_id,p0.run,p0.pair_name,p0.contact_a,p0.contact_b,T, ...
        'VariableNames',["subject" "session_id" "run" "pair_name" "contact_a" "contact_b" "T"]);
    p.bad_segments_s=badCell;
end
uniqueSubjects=unique(p.subject,"stable");uniqueSessions=unique(p.session_id,"stable");
etaSubject=base.background.log_scale_sigma.subject*randn(stream,numel(uniqueSubjects),1);
etaSession=base.background.log_scale_sigma.session*randn(stream,numel(uniqueSessions),1);
etaPair=base.background.log_scale_sigma.pair*randn(stream,height(p),1);
template=pairTemplate();pairsP0F=repmat(template,height(p),1);
source=covariates(logical(covariates.source_complete),:);
[blp,alp]=butter(base.sampling.anti_alias_order,base.sampling.anti_alias_lowpass_hz/(base.sampling.synthetic_hz/2));
contractSha=sha256File(fullfile(fileparts(mfilename("fullpath")),"stagec_mphase_gainflat_contract.json"));
peakMemory=currentMemoryBytes();clock=tic;
for i=1:height(p)
    si=find(uniqueSubjects==p.subject(i),1);qi=find(uniqueSessions==p.session_id(i),1);
    m=exp(etaSubject(si)+etaSession(qi)+etaPair(i));polarity=2*(rand(stream)>.5)-1;phi=2*pi*rand(stream);
    n1000=round(double(p.T(i))*base.sampling.synthetic_hz)+1;
    background=randomPhaseOneOverF(n1000,base.sampling.synthetic_hz,base.background.f_min_hz,stream);
    background250=filtfilt(blp,alp,background);background250=background250(1:4:end);
    background250=background250*(base.background.b0_uV*m/rms(background250));
    pairKey=p.subject(i)+"|"+p.session_id(i)+"|"+p.pair_name(i);% pair_name alone repeats across subjects
    [injF,epsF,nF,dF,sF]=phaseInjection(n1000,source,contract,0,phi,seed,"P0F",pairKey,contractSha);
    injF=filtfilt(blp,alp,injF);
    fp=sha256Text(jsonencode(struct("subject",p.subject(i),"session",p.session_id(i), ...
        "pair",p.pair_name(i),"T",double(p.T(i)),"eta_subject",etaSubject(si), ...
        "eta_session",etaSession(qi),"eta_pair",etaPair(i),"phi",phi)));
    q=template;q.subject=p.subject(i);q.session_id=p.session_id(i);q.run=p.run(i);q.pair_id=p.pair_name(i);
    q.contact_a=p.contact_a(i);q.contact_b=p.contact_b(i);q.T=double(p.T(i));q.bad_segments_s=p.bad_segments_s{i};
    q.eta_subject=etaSubject(si);q.eta_session=etaSession(qi);q.eta_pair=etaPair(i);q.background_scale_m=m;
    q.polarity=polarity;q.phi_u_rad=phi;q.base_fingerprint=fp;
    q.signal=background250+injF(1:4:end);q.epsilon_event_rad=epsF;q.n_injected_events=nF;q.injection_draw_count=dF;q.injection_stream_seed=sF;
    pairsP0F(i)=q;
    if mod(i,10)==0||i==height(p),peakMemory=max(peakMemory,currentMemoryBytes());end
end
world=struct("pairs",struct("P0F",pairsP0F),"covariates",covariates, ...
    "base_stream_state",stream.State,"seed",seed,"mode",mode,"contract_sha256",contractSha, ...
    "parent_contract_sha256",sha256File(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")), ...
    "truth",struct("A_G06_uV",double(contract.injection.A_G06_uV),"grid_id",string(contract.injection.grid_id), ...
    "kappa",0,"gain_slope_b",0.0,"frequency_hz",6), ...
    "timing",struct("synthesis_s",toc(clock),"peak_memory_bytes",peakMemory));
end

function q=pairTemplate()
q=struct("subject","","session_id","","run","","pair_id","","contact_a","","contact_b","", ...
    "T",NaN,"signal",[],"bad_segments_s",zeros(0,2),"eta_subject",NaN,"eta_session",NaN, ...
    "eta_pair",NaN,"background_scale_m",NaN,"polarity",NaN,"phi_u_rad",NaN, ...
    "epsilon_event_rad",[],"n_injected_events",0,"injection_draw_count",0,"injection_stream_seed",0,"base_fingerprint","");
end
function [y,epsEvent,nInjected,drawCount,streamSeed]=phaseInjection(n,ev,c,kappa,phi,worldSeed,arm,pairKey,contractSha)
fs=1000;off=round(.10*fs):round(.50*fs)-1;t=off/fs;h=sqrt(2)*sin(pi*(t-.10)/.40);
streamSeed=seedFromKey(contractSha+"|"+string(worldSeed)+"|"+arm+"|"+pairKey+"|phase_injection");
rs=RandStream("Threefry","Seed",streamSeed);[epsEvent,drawCount]=vonMises(rs,kappa,height(ev));y=zeros(n,1);nInjected=0;
for e=1:height(ev)
    ix=round(double(ev.label_time_s(e))*fs)+1+off;if ix(1)<1||ix(end)>n,continue,end
    g=exp(double(c.injection.gain_intercept_a)+double(c.injection.gain_slope_b)*double(ev.z_e(e)));
    y(ix)=y(ix)+double(c.injection.A_G06_uV)*g*h'.*sin(2*pi*double(c.injection.frequency_hz)*t'+phi+epsEvent(e));
    nInjected=nInjected+1;
end
end
function [x,drawCount]=vonMises(rs,kappa,n)
if kappa<=1e-12,x=2*pi*rand(rs,n,1)-pi;drawCount=n;return,end
a=1+sqrt(1+4*kappa^2);b=(a-sqrt(2*a))/(2*kappa);r=(1+b^2)/(2*b);x=zeros(n,1);drawCount=0;
for j=1:n
    while true,z=cos(pi*rand(rs));f=(1+r*z)/(r+z);cc=kappa*(r-f);u=rand(rs);drawCount=drawCount+2;if u<cc*(2-cc)||log(cc/u)+1>=cc,break,end,end
    x(j)=sign(rand(rs)-.5)*acos(f);drawCount=drawCount+1;
end
end
function x=randomPhaseOneOverF(n,fs,fmin,rs)
X=complex(zeros(n,1));k=(1:floor((n-1)/2))';f=k*fs/n;X(k+1)=1./max(f,fmin).*exp(1i*2*pi*rand(rs,numel(k),1));X(n-k+1)=conj(X(k+1));
if mod(n,2)==0,X(n/2+1)=(2*(rand(rs)>.5)-1)/(fs/2);end;x=ifft(X,"symmetric");x=x-mean(x);
end
function out=mergeIntervals(iv)
if isempty(iv),out=zeros(0,2);return,end;iv=sortrows(iv,[1 2]);out=iv(1,:);
for k=2:size(iv,1),if iv(k,1)<=out(end,2)+1e-12,out(end,2)=max(out(end,2),iv(k,2));else,out(end+1,:)=iv(k,:);end,end %#ok<AGROW>
end
function p=inputPath(root,rel)
% Contract input paths are absolute in the working tree and relative in the
% release snapshot; fullfile(root, absolute) would concatenate, so resolve
% absolute paths as-is first.
p=string(rel);if ~isfile(p),p=string(fullfile(root,p));end
end
function s=seedFromKey(key)
md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(uint8(char(key)),"int8"));b=typecast(md.digest(),"uint8");
s=double(b(1))+256*double(b(2))+65536*double(b(3))+16777216*double(b(4));if s<1,s=1;end
end
function bytes=currentMemoryBytes()
try
    u=memory;bytes=double(u.MemUsedMATLAB);
catch
    bytes=NaN;
end
end
function h=sha256Text(txt),md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(uint8(char(txt)),"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
function h=sha256File(path)
fid=fopen(path,"rb");assert(fid>=0,"Cannot open %s",path);cl=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance("SHA-256");
while true,b=fread(fid,1048576,"*uint8");if isempty(b),break,end;md.update(typecast(b,"int8"));end;h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end
