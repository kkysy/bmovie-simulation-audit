function world = acc00_sim_generate_world(root, contract, mode, seed, scenario, effect)
%ACC00_SIM_GENERATE_WORLD Synthesize traces; metrics stay in acc01_run_metrics.
arguments
    root (1,1) string
    contract (1,1) struct
    mode (1,1) string
    seed (1,1) double
    scenario (1,1) string = "null"
    effect (1,1) struct = struct()
end
scenario=lower(scenario); assert(ismember(scenario,["null" "additive" "aperiodic" "phase"]),"Unknown scenario.");
effect=fillEffect(effect); stream=RandStream("Threefry","Seed",seed);
covariates=readtable(contract.inputs.label_covariates.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
if mode=="smoke"
    covariates=covariates(1:12,:); covariates.label_time_s=linspace(2,13,12)';
    covariates.z_e=sin(2*pi*(0:11)'/3); covariates.source_complete(:)=true;
    p=table(repelem(["S01";"S02"],2),repelem(["S01_R1";"S02_R1"],2),repmat("R1",4,1), ...
        "pair"+(1:4)',repmat("A",4,1),repmat("B",4,1),repmat(20,4,1), ...
        'VariableNames',["subject" "session_id" "run" "pair_name" "contact_a" "contact_b" "T"]);
    p.bad_segments_s=repmat({zeros(0,2)},height(p),1);
else
    p0=readtable(contract.inputs.i00_acc_pairs.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
    p0=p0(logical(p0.qc_eligible)&contains(p0.location,"ACC"),:); assert(height(p0)==351,"Expected 351 ACC pairs.");
    timing=readtable(contract.inputs.i00_timing.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
    bad=readtable(contract.inputs.i00_bad_segments.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
    T=nan(height(p0),1); badCell=cell(height(p0),1);
    for i=1:height(p0)
        tr=timing(timing.session_id==p0.session_id(i),:); assert(height(tr)==1,"Timing join failed for %s.",p0.session_id(i));
        T(i)=double(tr.usable_movie_stop_s);
        br=bad(bad.session_id==p0.session_id(i)&bad.series=="LFP_macro"&ismember(bad.origchannel_name,[p0.contact_a(i),p0.contact_b(i)]),:);
        % I00 bad intervals use common_time. The first movie macro sample is at
        % common_movie_zero = macro_common_start_s +
        % (macro_movie_first_sample-1)/macro_rate_hz, so movie_time equals
        % common_time-common_movie_zero.
        movieZeroCommon=double(tr.macro_common_start_s)+(double(tr.macro_movie_first_sample)-1)/double(tr.macro_rate_hz);
        iv=[double(br.start_common_time_s),double(br.stop_common_time_s)]-movieZeroCommon;
        iv(:,1)=max(iv(:,1),0); iv(:,2)=min(iv(:,2),T(i)); badCell{i}=mergeIntervals(iv(iv(:,2)>iv(:,1),:));
    end
    p=table(p0.subject,p0.session_id,p0.run,p0.pair_name,p0.contact_a,p0.contact_b,T, ...
        'VariableNames',["subject" "session_id" "run" "pair_name" "contact_a" "contact_b" "T"]); p.bad_segments_s=badCell;
end
uniqueSubjects=unique(p.subject,"stable"); uniqueSessions=unique(p.session_id,"stable");
etaSubject=contract.background.log_scale_sigma.subject*randn(stream,numel(uniqueSubjects),1);
etaSession=contract.background.log_scale_sigma.session*randn(stream,numel(uniqueSessions),1);
etaPair=contract.background.log_scale_sigma.pair*randn(stream,height(p),1);
template=struct('subject',"",'session_id',"",'run',"",'pair_id',"",'contact_a',"",'contact_b',"",'T',NaN,'signal',[], ...
    'bad_segments_s',zeros(0,2),'eta_subject',NaN,'eta_session',NaN,'eta_pair',NaN,'background_scale_m',NaN, ...
    'polarity',NaN,'phi_u_rad',NaN,'epsilon_event_rad',[]); pairs=repmat(template,height(p),1);
peakMemory=currentMemoryBytes(); synthClock=tic;
for i=1:height(p)
    si=find(uniqueSubjects==p.subject(i),1); qi=find(uniqueSessions==p.session_id(i),1);
    m=exp(etaSubject(si)+etaSession(qi)+etaPair(i)); polarity=2*(rand(stream)>.5)-1; phi=2*pi*rand(stream);
    n1000=round(double(p.T(i))*contract.sampling.synthetic_hz)+1;
    background=randomPhaseOneOverF(n1000,contract.sampling.synthetic_hz,contract.background.f_min_hz,stream);
    [blp,alp]=butter(contract.sampling.anti_alias_order,contract.sampling.anti_alias_lowpass_hz/(contract.sampling.synthetic_hz/2));
    background250=filtfilt(blp,alp,background); background250=background250(1:4:end);
    background250=background250*(contract.background.b0_uV*m/rms(background250));
    injection=zeros(n1000,1); epsEvent=[]; source=covariates(logical(covariates.source_complete),:);
    switch scenario
        case "additive", injection=injectAdditive(injection,source,contract,effect.a_erp_uV,polarity);
        case "aperiodic", injection=injectAperiodic(injection,source,contract,effect.lambda_a_uV,stream);
        case "phase", [injection,epsEvent]=injectPhase(injection,source,contract,effect.a_erp_uV,effect.kappa,phi,stream);
    end
    if any(injection), injection=filtfilt(blp,alp,injection); signal=background250+injection(1:4:end); else, signal=background250; end
    pairs(i)=struct('subject',p.subject(i),'session_id',p.session_id(i),'run',p.run(i),'pair_id',p.pair_name(i), ...
        'contact_a',p.contact_a(i),'contact_b',p.contact_b(i),'T',double(p.T(i)),'signal',signal,'bad_segments_s',p.bad_segments_s{i}, ...
        'eta_subject',etaSubject(si),'eta_session',etaSession(qi),'eta_pair',etaPair(i),'background_scale_m',m, ...
        'polarity',polarity,'phi_u_rad',phi,'epsilon_event_rad',epsEvent);
    if mod(i,10)==0||i==height(p), peakMemory=max(peakMemory,currentMemoryBytes()); end
end
truth=struct('scenario',scenario,'a_erp_uV',effect.a_erp_uV,'lambda_a_uV',effect.lambda_a_uV, ...
    'kappa',effect.kappa,'provisional',logical(effect.provisional),'gain_intercept_a',contract.additive.gain_intercept_a, ...
    'gain_slope_b',contract.additive.gain_slope_b);
world=struct('pairs',pairs,'covariates',covariates,'stream',stream,'seed',seed,'mode',mode,'truth',truth, ...
    'contract_sha256',sha256File(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")), ...
    'timing',struct('synthesis_s',toc(synthClock),'peak_memory_bytes',peakMemory));

function e=fillEffect(e)
d=struct('a_erp_uV',0,'lambda_a_uV',0,'kappa',0,'provisional',false);
for f=string(fieldnames(d))', if ~isfield(e,f),e.(f)=d.(f);end,end
end
function x=randomPhaseOneOverF(n,fs,fmin,rs)
X=complex(zeros(n,1)); k=(1:floor((n-1)/2))'; f=k*fs/n; X(k+1)=1./max(f,fmin).*exp(1i*2*pi*rand(rs,numel(k),1)); X(n-k+1)=conj(X(k+1));
if mod(n,2)==0,X(n/2+1)=(2*(rand(rs)>.5)-1)/(fs/2);end; x=ifft(X,"symmetric"); x=x-mean(x);
end
function y=injectAdditive(y,ev,c,A,pol)
fs=c.sampling.synthetic_hz; off=round(.10*fs):round(.50*fs)-1; t=off/fs; h=sqrt(2)*sin(pi*(t-.10)/.40);
for e=1:height(ev),ix=round(double(ev.label_time_s(e))*fs)+1+off;if ix(1)<1||ix(end)>numel(y),continue,end;g=exp(c.additive.gain_intercept_a+c.additive.gain_slope_b*double(ev.z_e(e)));y(ix)=y(ix)+A*g*pol*h';end
end
function y=injectAperiodic(y,ev,c,L,rs)
fs=c.sampling.synthetic_hz; off=round(.10*fs):round(.50*fs)-1;t=off/fs;h=sqrt(2)*sin(pi*(t-.10)/.40)';
for e=1:height(ev),ix=round(double(ev.label_time_s(e))*fs)+1+off;if ix(1)<1||ix(end)>numel(y),continue,end;q=randomPhaseOneOverF(numel(off),fs,c.background.f_min_hz,rs).*h;q=q/rms(q);g=exp(c.additive.gain_intercept_a+c.additive.gain_slope_b*double(ev.z_e(e)));y(ix)=y(ix)+L*g*q;end
end
function [y,epsEvent]=injectPhase(y,ev,c,A,kappa,phi,rs)
fs=c.sampling.synthetic_hz;off=round(.10*fs):round(.50*fs)-1;t=off/fs;h=sqrt(2)*sin(pi*(t-.10)/.40);epsEvent=vonMises(rs,kappa,height(ev));
for e=1:height(ev),ix=round(double(ev.label_time_s(e))*fs)+1+off;if ix(1)<1||ix(end)>numel(y),continue,end;y(ix)=y(ix)+A*h'.*sin(2*pi*6*t'+phi+epsEvent(e));end
end
function x=vonMises(rs,kappa,n)
if kappa<=1e-12,x=2*pi*rand(rs,n,1)-pi;return,end;a=1+sqrt(1+4*kappa^2);b=(a-sqrt(2*a))/(2*kappa);r=(1+b^2)/(2*b);x=zeros(n,1);
for j=1:n,while true,z=cos(pi*rand(rs));f=(1+r*z)/(r+z);cc=kappa*(r-f);u=rand(rs);if u<cc*(2-cc)||log(u)<=cc-1,break,end,end;x(j)=sign(rand(rs)-.5)*acos(f);end
end
function out=mergeIntervals(iv)
if isempty(iv),out=zeros(0,2);return,end;iv=sortrows(iv,[1 2]);out=iv(1,:);for k=2:size(iv,1),if iv(k,1)<=out(end,2)+1e-12,out(end,2)=max(out(end,2),iv(k,2));else,out(end+1,:)=iv(k,:);end,end %#ok<AGROW>
end
function bytes=currentMemoryBytes()
try
    u=memory;bytes=double(u.MemUsedMATLAB);
catch
    bytes=NaN;
end
end
function hash=sha256File(path)
fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance("SHA-256");while true,q=fread(fid,1048576,"*uint8");if isempty(q),break,end;md.update(typecast(q,"int8"));end;hash=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end
end
