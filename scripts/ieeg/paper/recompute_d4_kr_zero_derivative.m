function result = recompute_d4_kr_zero_derivative(root)
% D4-2 faithful time-domain finite-difference probe (diagnostic only).
% Three deterministic synthetic worlds share one seed; no neural data or
% pooled simulation is run. Pair slopes are equally averaged for stability.
root = string(char(java.io.File(char(root)).getCanonicalPath()));
addpath(fullfile(root,"scripts","ieeg","acc01"));
c = jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
seed = double(c.randomization.base_seed) + 100000*2 + 1000*3 + 1;
A = [0 0.25 0.5]; slopes = nan(size(A)); meanInvM2 = NaN;
for k = 1:numel(A)
    w = acc00_sim_generate_world(root,c,"benchmark",seed,"additive",struct("a_erp_uV",A(k),"provisional",false));
    if isnan(meanInvM2), meanInvM2 = mean(arrayfun(@(p) 1/double(p.background_scale_m)^2,w.pairs)); end
    pairSlope = nan(numel(w.pairs),1);
    for p = 1:numel(w.pairs)
        [events,~,~] = acc01_build_support_and_events(w.pairs(p),w.covariates,c);
        q = acc01_compute_power_metric(w.pairs(p),events,c);
        pairSlope(p) = q.slope;
    end
    assert(all(isfinite(pairSlope)),"Nonfinite probe slope"); slopes(k)=mean(pairSlope);
end
x=A.^2; d01=(slopes(2)-slopes(1))/(x(2)-x(1)); d12=(slopes(3)-slopes(2))/(x(3)-x(2));
cz=double(c.label_gain_diagnostics.cov_z_exp2bz_over_var_z);
kr01=d01/(cz*meanInvM2); kr12=d12/(cz*meanInvM2); kr=mean([kr01 kr12]);
fs=double(c.sampling.analysis_hz); n=round(double(c.power.hann_epoch_s)*fs); h=sqrt(2)*sin(pi*((0:n-1)'/fs)/.4); h=h-mean(h); hannw=hann(n); ph=abs(fft(h.*hannw)).^2/sum(hannw.^2);
theta=double(c.power.theta_periodogram_bins_hz(:)); idx=round(theta*n/fs)+1; bg=[1469.3;653.0]; upper=mean(ph(idx)./(log(10)*bg));
d4=1.0929e-4; ratio=abs(kr-d4)/d4;
result=struct("status",ternary(kr>0 && kr<=upper && ratio<=2,"PASS","STOP"),"seed",seed,"A_uV",A,"pair_mean_slopes",slopes,"d_slope_d_A2_01",d01,"d_slope_d_A2_12",d12,"finite_difference_relative_disagreement",abs(d01-d12)/max(abs([d01 d12])) ,"mean_m_inverse2",meanInvM2,"C_z",cz,"recomputed_K_R",kr,"K_R_01",kr01,"K_R_12",kr12,"unadjusted_upper_bound",upper,"d4_historical_K_R",d4,"relative_to_D4",ratio,"background_theta_uV2",bg,"kernel_theta_power",ph(idx),"source","acc00_sim_generate_world + acc01_compute_power_metric; contract half_sine_rms1; D4-2 design draft section 3.2.4");
out=fullfile(root,"processed","subject","group","ieeg_methods_paper","p0_2_design"); if ~isfolder(out),mkdir(out);end
fid=fopen(fullfile(out,"d4_2_kr_recompute.tsv"),"w"); fprintf(fid,"seed\tA0_slope\tA1_slope\tA2_slope\td01\td12\tK_R\tupper\td4\trelative_to_D4\tstatus\n%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%s\n",seed,slopes,d01,d12,kr,upper,d4,ratio,result.status); fclose(fid);
fid=fopen(fullfile(out,"d4_2_kr_recompute.json"),"w"); fprintf(fid,"%s\n",jsonencode(result)); fclose(fid);
if result.status=="STOP", error("K_R probe failed acceptance interval"); end
end
function y=ternary(test,a,b),if test,y=a;else,y=b;end,end
