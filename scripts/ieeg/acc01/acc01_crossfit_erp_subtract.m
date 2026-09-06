function [residual,diag] = acc01_crossfit_erp_subtract(pair, events, contract, validateReference)
% Fold-level reusable operator: one event-by-offset E, one solve per fold.
if nargin<4,validateReference=false;end
residual=pair.signal;fs=contract.sampling.analysis_hz;off=round(contract.windows_s.filter_support(1)*fs):round(contract.windows_s.filter_support(2)*fs);
n=height(events);
if ismember("crossfit_fold", string(events.Properties.VariableNames))
    fold=double(events.crossfit_fold(:));
    assert(numel(fold)==n&&all(isfinite(fold))&&all(fold>=1&fold<=contract.phase.erp_crossfit_folds), ...
        "crossfit_fold must carry a valid frozen fold identity.");
else
    fold=min(contract.phase.erp_crossfit_folds,ceil((1:n)'*contract.phase.erp_crossfit_folds/max(n,1)));
end
cols=endsWith(events.Properties.VariableNames,"_delta_z")|endsWith(events.Properties.VariableNames,"movie_time_linear_z")|endsWith(events.Properties.VariableNames,"movie_time_quadratic_z");
X=[ones(n,1),double(events.z_e),double(events{:,cols})];center=round(double(events.label_time_s)*fs)+1;
indices=center+off;E=pair.signal(indices);prediction=zeros(size(E));operatorSeconds=0;
for k=1:max(fold)
    tr=fold~=k;te=find(fold==k);if nnz(tr)<size(X,2),continue,end
    q=tic;B=X(tr,:)\E(tr,:);prediction(te,:)=X(te,:)*B;operatorSeconds=operatorSeconds+toc(q);
end
for i=1:n,residual(indices(i,:))=residual(indices(i,:))-prediction(i,:)';end
relativeError=NaN;
if validateReference
    legacyEpoch=E;
    for k=1:max(fold)
        tr=fold~=k;te=find(fold==k);if nnz(tr)<size(X,2),continue,end
        for s=1:numel(off)
            bh=X(tr,:)\E(tr,s);
            legacyEpoch(te,s)=E(te,s)-X(te,:)*bh;
        end
    end
    optimizedEpoch=E-prediction;
    relativeError=norm(optimizedEpoch-legacyEpoch,"fro")/max(norm(E,"fro"),realmin);
end
diag=struct('n_events',n,'n_offsets',numel(off),'operator_seconds',operatorSeconds, ...
    'reference_relative_error',relativeError);
end
