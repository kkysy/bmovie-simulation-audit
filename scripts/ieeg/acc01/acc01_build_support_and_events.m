function [powerEvents, phaseEvents, diag] = acc01_build_support_and_events(pair, covariates, contract)
% Shared real/synthetic event and strict-support selector; never thresholds z.
t = covariates.label_time_s; z = covariates.z_e; good = covariates.source_complete == 1;
powerGood = good & t+contract.windows_s.power_post(2) <= pair.T & t+contract.windows_s.power_pre(1) >= 0;
phaseGood = good & t+contract.windows_s.filter_support(2) <= pair.T & t+contract.windows_s.filter_support(1) >= 0;
for k=1:size(pair.bad_segments_s,1), b=pair.bad_segments_s(k,:); powerGood=powerGood & ~(t+contract.windows_s.power_pre(1)<b(2) & t+contract.windows_s.power_post(2)>b(1)); phaseGood=phaseGood & ~(t+contract.windows_s.filter_support(1)<b(2) & t+contract.windows_s.filter_support(2)>b(1)); end
mx = false(size(z)); mx(2:end-1) = z(2:end-1)>z(1:end-2) & z(2:end-1)>=z(3:end); % earliest member of a plateau
powerEvents = covariates(powerGood,:); phaseEvents = covariates(phaseGood & mx,:);
diag = struct('n_power',height(powerEvents),'n_phase',height(phaseEvents),'power_coverage',mean(powerGood),'phase_coverage',mean(phaseGood));
end
