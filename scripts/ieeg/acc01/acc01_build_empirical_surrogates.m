function out = acc01_build_empirical_surrogates(pair, events, contract, stream)
% Analytic D_u interval construction followed by exact 4-ms quantization.
T=double(pair.T); q=double(contract.surrogate.delta_quantum_s);
% The domain's movie-end guards are deliberately separate from the unchanged
% 0.75-s filter support and bad-segment dilation.  They are expressed in the
% contract and applied directly on the same closed-interval grid convention.
edgeGuardLower=double(contract.surrogate.domain_guard_edge_lower_s);
edgeGuardUpper=double(contract.surrogate.domain_guard_edge_upper_s);
badDilation=double(contract.surrogate.domain_guard_bad_dilation_s);
bad=[0 edgeGuardLower;T-edgeGuardUpper T];
for k=1:size(pair.bad_segments_s,1),bad=[bad;pair.bad_segments_s(k,:)+[-badDilation badDilation]];end %#ok<AGROW>
bad(:,1)=max(bad(:,1),0);bad(:,2)=min(bad(:,2),T);bad=mergeIntervals(bad(bad(:,2)>bad(:,1),:));
forbidden=zeros(0,2); et=double(events.label_time_s(:));
for e=1:numel(et)
    for b=1:size(bad,1)
        a=bad(b,1)-et(e); z=bad(b,2)-et(e); len=z-a;
        if len>=T-1e-12,forbidden=[0 T];break,end
        a=mod(a,T);z=a+len;
        if z<=T,forbidden(end+1,:)=[a z];else,forbidden(end+1,:)=[a T];forbidden(end+1,:)=[0 z-T];end %#ok<AGROW>
    end
    if isequal(forbidden,[0 T]),break,end
end
forbidden=mergeIntervals(forbidden); feasibleIntervals=complementIntervals(forbidden,T);
analyticDuration=sum(feasibleIntervals(:,2)-feasibleIntervals(:,1));
% Frozen grid convention is exactly grid=(0:q:T-q), and B_u intervals are
% closed. Mark analytic forbidden intervals on that grid without rejection.
nGrid=floor((T-q)/q)+1; forbiddenGrid=false(nGrid,1);
for k=1:size(forbidden,1)
    j0=max(0,ceil(forbidden(k,1)/q-1e-12));j1=min(nGrid-1,floor(forbidden(k,2)/q+1e-12));
    if j1>=j0,forbiddenGrid(j0+1:j1+1)=true;end
end
gridIndex=find(~forbiddenGrid)-1;feasible=gridIndex*q;n=numel(feasible);domainDuration=n*q;
if n==0,d=[];mode="not_estimable";
elseif n<=contract.surrogate.enumerate_when_quantized_feasible_count_lte,d=feasible;mode="enumerated";
else,d=feasible(randperm(stream,n,contract.surrogate.max_surrogates_per_estimable_pair));mode="sampled";
end
out=struct('delta_s',d,'quantized_feasible_count',n,'domain_fraction',n/nGrid, ...
    'domain_duration_s',domainDuration,'analytic_domain_duration_s',analyticDuration, ...
    'interval_count',size(feasibleIntervals,1), ...
    'feasible_intervals_s',feasibleIntervals,'sample_count',numel(d),'mode',mode);

function out=mergeIntervals(iv)
if isempty(iv),out=zeros(0,2);return,end;iv=sortrows(iv,[1 2]);out=iv(1,:);
for j=2:size(iv,1),if iv(j,1)<=out(end,2)+1e-12,out(end,2)=max(out(end,2),iv(j,2));else,out(end+1,:)=iv(j,:);end,end %#ok<AGROW>
end
function out=complementIntervals(iv,total)
if isempty(iv),out=[0 total];return,end;out=zeros(0,2);cursor=0;
for j=1:size(iv,1),if iv(j,1)>cursor+1e-12,out(end+1,:)=[cursor iv(j,1)];end;cursor=max(cursor,iv(j,2));end %#ok<AGROW>
if cursor<total-1e-12,out(end+1,:)=[cursor total];end
end
end
