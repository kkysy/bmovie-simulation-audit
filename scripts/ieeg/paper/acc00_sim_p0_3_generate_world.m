function [world, provenance] = acc00_sim_p0_3_generate_world(root, paper, base, cellId, seedIndex, mode)
%ACC00_SIM_P0_3_GENERATE_WORLD Synthetic-only P0-3 injection wrapper.
% P01 delegates to the immutable P0-2 generator. Other cells reproduce its
% main Threefry stream and alter injection only; no LFP reader is reachable.
arguments
    root (1,1) string
    paper (1,1) struct
    base (1,1) struct
    cellId (1,1) string
    seedIndex (1,1) double {mustBeInteger,mustBePositive}
    mode (1,1) string {mustBeMember(mode,["smoke" "synthetic"])} = "synthetic"
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
cell=powerCell(paper,cellId);assert(seedIndex<=double(cell.worlds),"seed_index exceeds cell size.");
seed=double(paper.seed.P01_to_P10.paper_base_seed)+100000*double(paper.seed.P01_to_P10.scenario_index)+1000*double(paper.seed.P01_to_P10.effect_index)+seedIndex;
effect=struct("a_erp_uV",double(paper.additive.A_uV),"provisional",false);
if cellId=="P01"
    % Canonical bitwise route: no copied synthesis or injection code is used.
    world=acc00_sim_generate_world(root,base,legacyMode(mode),seed,"additive",effect);
    provenance=makeProvenance(world.covariates,cell,paper,seed,world.pairs,[]);
    [a,s,pa,ps]=countInjectedEvents(world.covariates,provenance.sorted_selected_source_row_indices,cell,paper,world.pairs);
    provenance.kernel_attempted_event_count=a;provenance.kernel_skipped_out_of_bounds_count=s;provenance.per_pair_attempted_event_count=pa;provenance.per_pair_skipped_out_of_bounds_count=ps;
    world.p0_3=struct("cell_id",cellId,"seed_index",seedIndex,"world_seed",seed,"injection_provenance",provenance,"wrapper_path","canonical_delegate");
    return
end
covariates=readtable(fullfile(root,string(paper.inputs.label_covariates.path)),"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
[p,covariates]=syntheticInputs(root,base,mode,covariates);
stream=RandStream("Threefry","Seed",seed); mask=makeMask(covariates,cell,paper,seed);
uniqueSubjects=unique(p.subject,"stable");uniqueSessions=unique(p.session_id,"stable");
etaSubject=base.background.log_scale_sigma.subject*randn(stream,numel(uniqueSubjects),1);
etaSession=base.background.log_scale_sigma.session*randn(stream,numel(uniqueSessions),1);
etaPair=base.background.log_scale_sigma.pair*randn(stream,height(p),1);
template=struct('subject',"",'session_id',"",'run',"",'pair_id',"",'contact_a',"",'contact_b',"",'T',NaN,'signal',[], ...
    'bad_segments_s',zeros(0,2),'eta_subject',NaN,'eta_session',NaN,'eta_pair',NaN,'background_scale_m',NaN, ...
    'polarity',NaN,'phi_u_rad',NaN,'epsilon_event_rad',[]);pairs=repmat(template,height(p),1);peakMemory=currentMemoryBytes();clock=tic;
[blp,alp]=butter(base.sampling.anti_alias_order,base.sampling.anti_alias_lowpass_hz/(base.sampling.synthetic_hz/2));
attempted=0;skipped=0;pairAttempted=zeros(height(p),1);pairSkipped=zeros(height(p),1);
for i=1:height(p)
    si=find(uniqueSubjects==p.subject(i),1);qi=find(uniqueSessions==p.session_id(i),1);m=exp(etaSubject(si)+etaSession(qi)+etaPair(i));polarity=2*(rand(stream)>.5)-1;phi=2*pi*rand(stream);
    n1000=round(double(p.T(i))*base.sampling.synthetic_hz)+1;background=randomPhaseOneOverF(n1000,base.sampling.synthetic_hz,base.background.f_min_hz,stream);
    background250=filtfilt(blp,alp,background);background250=background250(1:4:end);background250=background250*(base.background.b0_uV*m/rms(background250));
    injection=zeros(n1000,1);[injection,a,s]=injectCell(injection,covariates,mask,cell,paper,polarity,seed);attempted=attempted+a;skipped=skipped+s;pairAttempted(i)=a;pairSkipped(i)=s;
    if any(injection),injection=filtfilt(blp,alp,injection);signal=background250+injection(1:4:end);else,signal=background250;end
    pairs(i)=struct('subject',p.subject(i),'session_id',p.session_id(i),'run',p.run(i),'pair_id',p.pair_name(i),'contact_a',p.contact_a(i),'contact_b',p.contact_b(i),'T',double(p.T(i)), ...
        'signal',signal,'bad_segments_s',p.bad_segments_s{i},'eta_subject',etaSubject(si),'eta_session',etaSession(qi),'eta_pair',etaPair(i),'background_scale_m',m,'polarity',polarity,'phi_u_rad',phi,'epsilon_event_rad',[]);
    if mod(i,10)==0||i==height(p),peakMemory=max(peakMemory,currentMemoryBytes());end
end
truth=struct('scenario',"additive",'a_erp_uV',effect.a_erp_uV,'lambda_a_uV',0,'kappa',0,'provisional',false,'gain_intercept_a',base.additive.gain_intercept_a,'gain_slope_b',base.additive.gain_slope_b);
world=struct('pairs',pairs,'covariates',covariates,'stream',stream,'seed',seed,'mode',mode,'truth',truth,'contract_sha256',sha256File(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")), ...
    'timing',struct('synthesis_s',toc(clock),'peak_memory_bytes',peakMemory));
provenance=makeProvenance(covariates,cell,paper,seed,pairs,mask);provenance.kernel_attempted_event_count=attempted;provenance.kernel_skipped_out_of_bounds_count=skipped;provenance.per_pair_attempted_event_count=pairAttempted;provenance.per_pair_skipped_out_of_bounds_count=pairSkipped;
world.p0_3=struct("cell_id",cellId,"seed_index",seedIndex,"world_seed",seed,"injection_provenance",provenance,"wrapper_path","synthetic_rebuilt_injection_only");
end

function cell=powerCell(paper,id)
z=paper.cells.Mpower;ix=find(string({z.id})==id,1);assert(~isempty(ix),"Unknown P0-3 Mpower cell %s.",id);cell=z(ix);
end
function out=legacyMode(mode),if mode=="smoke",out="smoke";else,out="benchmark";end,end
function [p,c]=syntheticInputs(root,base,mode,c)
if mode=="smoke"
    c=c(1:12,:);c.label_time_s=linspace(2,13,12)';c.z_e=sin(2*pi*(0:11)'/3);c.source_complete(:)=true;
    p=table(repelem(["S01";"S02"],2),repelem(["S01_R1";"S02_R1"],2),repmat("R1",4,1),"pair"+(1:4)',repmat("A",4,1),repmat("B",4,1),repmat(20,4,1), ...
        'VariableNames',["subject" "session_id" "run" "pair_name" "contact_a" "contact_b" "T"]);p.bad_segments_s=repmat({zeros(0,2)},height(p),1);return
end
p0=readtable(base.inputs.i00_acc_pairs.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");p0=p0(logical(p0.qc_eligible)&contains(p0.location,"ACC"),:);assert(height(p0)==351,"Expected 351 ACC pairs.");
timing=readtable(base.inputs.i00_timing.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");bad=readtable(base.inputs.i00_bad_segments.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");T=nan(height(p0),1);badCell=cell(height(p0),1);
for i=1:height(p0)
    tr=timing(timing.session_id==p0.session_id(i),:);assert(height(tr)==1,"Timing join failed.");T(i)=double(tr.usable_movie_stop_s);br=bad(bad.session_id==p0.session_id(i)&bad.series=="LFP_macro"&ismember(bad.origchannel_name,[p0.contact_a(i),p0.contact_b(i)]),:);
    movieZero=double(tr.macro_common_start_s)+(double(tr.macro_movie_first_sample)-1)/double(tr.macro_rate_hz);iv=[double(br.start_common_time_s),double(br.stop_common_time_s)]-movieZero;iv(:,1)=max(iv(:,1),0);iv(:,2)=min(iv(:,2),T(i));badCell{i}=mergeIntervals(iv(iv(:,2)>iv(:,1),:));
end
p=table(p0.subject,p0.session_id,p0.run,p0.pair_name,p0.contact_a,p0.contact_b,T,'VariableNames',["subject" "session_id" "run" "pair_name" "contact_a" "contact_b" "T"]);p.bad_segments_s=badCell;
end
function mask=makeMask(covariates,cell,paper,seed)
source=find(logical(covariates.source_complete));N=numel(source);r=double(paper.event_mask.levels.(char(cell.retention)));K=round(r*N);rs=RandStream("Threefry","Seed",seed);rs.Substream=double(paper.event_mask.rng.substream);take=sort(randperm(rs,N,K));selected=source(take);vector=false(height(covariates),1);vector(selected)=true;
mask=struct('source_rows',source,'vector',vector,'selected_rows',selected,'N',N,'K',K,'retention_fraction',r,'sha256',sha256Bytes(uint8(vector)));
end
function [y,attempted,skipped]=injectCell(y,covariates,mask,pCell,paper,polarity,seed)
fs=double(paper.sampling.synthetic_hz);selected=covariates(mask.vector,:);attempted=height(selected);skipped=0;phase=[];
if string(pCell.spectral_shape)=="theta_centered_tapered_random_phase_carrier"
    rs=RandStream("Threefry","Seed",seed);rs.Substream=double(paper.spectral_shape.P10.phase_rng.substream);phase=2*pi*rand(rs,height(selected),1);
end
for e=1:height(selected)
    [components,amps]=kernelComponents(pCell,paper,phase,e);valid=true;ixs=cell(numel(components),1);
    for k=1:numel(components)
        off=components{k}.off;ix=round(double(selected.label_time_s(e))*fs)+1+off;if ix(1)<1||ix(end)>numel(y),valid=false;break,end;ixs{k}=ix;
    end
    if ~valid,skipped=skipped+1;continue,end
    gain=exp(double(paper.additive.gain_intercept_a)+double(paper.additive.gain_slope_b)*double(selected.z_e(e)));
    for k=1:numel(components),y(ixs{k})=y(ixs{k})+amps(k)*gain*polarity*components{k}.h;end
end
end
function [components,amps]=kernelComponents(cell,paper,phase,e)
fs=double(paper.sampling.synthetic_hz);post=round(.10*fs):round(.50*fs)-1;pre=round(-.50*fs):round(-.10*fs)-1;
hpost=sqrt(2)*sin(pi*(post/fs-.10)/.40)';hpre=sqrt(2)*sin(pi*(pre/fs+.50)/.40)';A=double(paper.additive.A_uV);
if string(cell.spectral_shape)=="theta_centered_tapered_random_phase_carrier"
    tau=post/fs-.10;v=sin(pi*tau/.40).*cos(2*pi*double(paper.spectral_shape.P10.carrier_hz)*tau+phase(e));hpost=(v/rms(v))';
end
if string(cell.kernel_position)=="post",components={struct('off',post,'h',hpost)};amps=A;
elseif string(cell.kernel_position)=="pre",components={struct('off',pre,'h',hpre)};amps=A;
else,components={struct('off',pre,'h',hpre),struct('off',post,'h',hpost)};amps=[A A]/sqrt(2);end
end
function p=makeProvenance(covariates,cell,paper,seed,pairs,mask)
if isempty(mask),mask=makeMask(covariates,cell,paper,seed);end
times=double(covariates.label_time_s(mask.selected_rows));g=diff(sort(times));q=nan(1,5);if ~isempty(g),q=prctile(g,[0 25 50 75 100],"Method","inclusive");end
if string(cell.kernel_position)=="pre",overlap=NaN;num=NaN;den=NaN;reason="no_post_component";else,overlap=NaN;num=NaN;den=NaN;reason="runner_pools_analysis_power_events";end
p=struct('kernel_position',string(cell.kernel_position),'spectral_shape',string(cell.spectral_shape),'A_uV',double(paper.additive.A_uV), ...
    'component_amplitude_uV',componentAmplitude(cell,paper),'mask_rng_algorithm',"Threefry",'mask_rng_seed',seed,'mask_rng_substream',double(paper.event_mask.rng.substream), ...
    'N_source_complete',mask.N,'retention_fraction',mask.retention_fraction,'K_retained',mask.K,'mask_sha256',mask.sha256,'sorted_selected_source_row_indices',double(mask.selected_rows(:)), ...
    'selected_gap_n',numel(g),'selected_gap_min_s',q(1),'selected_gap_q25_s',q(2),'selected_gap_median_s',q(3),'selected_gap_q75_s',q(4),'selected_gap_max_s',q(5), ...
    'predecessor_in_pre_overlap_n',num,'predecessor_in_pre_overlap_denominator',den,'predecessor_in_pre_overlap_fraction',overlap,'predecessor_in_pre_overlap_reason',reason, ...
    'kernel_attempted_event_count',NaN,'kernel_skipped_out_of_bounds_count',NaN,'per_pair_attempted_event_count',[],'per_pair_skipped_out_of_bounds_count',[],'pair_count',numel(pairs));
end
function [attempted,skipped,pairAttempted,pairSkipped]=countInjectedEvents(covariates,selectedRows,pCell,paper,pairs)
fs=double(paper.sampling.synthetic_hz);selected=covariates(double(selectedRows),:);attempted=0;skipped=0;pairAttempted=zeros(numel(pairs),1);pairSkipped=zeros(numel(pairs),1);
for i=1:numel(pairs)
    for e=1:height(selected)
        [components,~]=kernelComponents(pCell,paper,zeros(height(selected),1),e);pairAttempted(i)=pairAttempted(i)+1;valid=true;
        for k=1:numel(components),ix=round(double(selected.label_time_s(e))*fs)+1+components{k}.off;if ix(1)<1||ix(end)>round(double(pairs(i).T)*fs)+1,valid=false;break,end,end
        if ~valid,pairSkipped(i)=pairSkipped(i)+1;end
    end
end
attempted=sum(pairAttempted);skipped=sum(pairSkipped);
end
function a=componentAmplitude(cell,paper),if string(cell.kernel_position)=="symmetric",a=double(paper.additive.A_uV)/sqrt(2);else,a=double(paper.additive.A_uV);end,end
function out=ternary(c,a,b),if c,out=a;else,out=b;end,end
function x=randomPhaseOneOverF(n,fs,fmin,rs),X=complex(zeros(n,1));k=(1:floor((n-1)/2))';f=k*fs/n;X(k+1)=1./max(f,fmin).*exp(1i*2*pi*rand(rs,numel(k),1));X(n-k+1)=conj(X(k+1));if mod(n,2)==0,X(n/2+1)=(2*(rand(rs)>.5)-1)/(fs/2);end;x=ifft(X,"symmetric");x=x-mean(x);end
function out=mergeIntervals(iv),if isempty(iv),out=zeros(0,2);return,end;iv=sortrows(iv,[1 2]);out=iv(1,:);for k=2:size(iv,1),if iv(k,1)<=out(end,2)+1e-12,out(end,2)=max(out(end,2),iv(k,2));else,out(end+1,:)=iv(k,:);end,end,end %#ok<AGROW>
function bytes=currentMemoryBytes(),try,u=memory;bytes=double(u.MemUsedMATLAB);catch,bytes=NaN;end,end
function h=sha256File(path),fid=fopen(path,"rb");assert(fid>=0);cl=onCleanup(@()fclose(fid));h=sha256Bytes(fread(fid,Inf,"*uint8"));end
function h=sha256Bytes(bytes),md=java.security.MessageDigest.getInstance("SHA-256");md.update(typecast(uint8(bytes(:)),"int8"));h=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));end
