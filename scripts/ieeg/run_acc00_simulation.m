function outputs = run_acc00_simulation(varargin)
%RUN_ACC00_SIMULATION ACC00-SIM stage runner.  Full simulation is disabled.

p = inputParser;
p.addParameter("ProjectRoot", "", @(x)ischar(x) || isstring(x));
p.addParameter("Stage", "covariates", @(x)ischar(x) || isstring(x));
p.addParameter("BenchmarkWorlds", 0, @(x)isnumeric(x) && isscalar(x));
p.parse(varargin{:});
root = string(p.Results.ProjectRoot);
if strlength(root) == 0
    root = string(fileparts(fileparts(fileparts(mfilename("fullpath")))));
end
root = string(char(java.io.File(char(root)).getCanonicalPath()));
stage = lower(string(p.Results.Stage));

switch stage
    case "covariates"
        outputs = buildCovariates(root);
    case "contract"
        outputs = writeContract(root);
    case "smoke"
        outputs = runSmoke(root);
    case "smoke_optimized"
        outputs = runSmokeOptimized(root);
    case "benchmark"
        outputs = runBenchmark(root);
    case "benchmark_optimized"
        outputs = runBenchmarkOptimized(root);
    case "nullpool"
        outputs = run_acc00_nullpool(root);
    otherwise
        error("Bmovie:ACC00SIM:UnauthorizedStage", ...
            "Stage '%s' is not enabled by this entry point.", stage);
end

function outputs = runBenchmark(root)
addpath(fullfile(root,"scripts","ieeg","acc01"));c=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
scenarios=["null";"additive";"aperiodic";"phase";"phase"];
labels=["null";"additive@Q75";"aperiodic@Q75_lambda_provisional";"phase@Q75_kappa1_provisional";"phase@Q75_kappa4_provisional"];
q75=c.additive.a_erp_grid_uV(3);lambda=c.aperiodic.benchmark_lambda_provisional_uV;
effects={struct();struct('a_erp_uV',q75);struct('a_erp_uV',q75,'lambda_a_uV',lambda,'provisional',true); ...
    struct('a_erp_uV',q75,'kappa',1,'provisional',true);struct('a_erp_uV',q75,'kappa',4,'provisional',true)};
seeds=c.randomization.base_seed+[100001;203001;303001;403001;403002];
outDir=fullfile(root,"processed","subject","group","ieeg_acc00_sim","BangYoureDead","checkpoints","benchmark");if ~isfolder(outDir),mkdir(outDir),end
sec=nan(5,1);synthesis_s=sec;power_s=sec;crossfit_s=sec;surrogate_ppc_s=sec;domain_s=sec;peak_memory_bytes=sec;status=repmat("running",5,1);
for i=1:5
    worldClock=tic;w=acc00_sim_generate_world(root,c,"benchmark",seeds(i),scenarios(i),effects{i});r=acc01_run_metrics(w,c);sec(i)=toc(worldClock);
    synthesis_s(i)=w.timing.synthesis_s;power_s(i)=r.timing.power_s;crossfit_s(i)=r.timing.crossfit_s;
    surrogate_ppc_s(i)=r.timing.surrogate_ppc_s;domain_s(i)=r.timing.domain_s;peak_memory_bytes(i)=w.timing.peak_memory_bytes;status(i)="complete";
    checkpoint=struct('scenario',scenarios(i),'label',labels(i),'seed',seeds(i),'truth',w.truth,'contract_sha256',w.contract_sha256, ...
        'elapsed_s',sec(i),'synthesis_s',synthesis_s(i),'power_s',power_s(i),'crossfit_s',crossfit_s(i), ...
        'surrogate_ppc_s',surrogate_ppc_s(i),'domain_s',domain_s(i),'peak_memory_bytes',peak_memory_bytes(i),'status',status(i),'result',r);
    tmp=fullfile(outDir,sprintf("world_%02d.tmp.mat",i));final=fullfile(outDir,sprintf("world_%02d.mat",i));save(tmp,"checkpoint","-v7.3");movefile(tmp,final,"f");
end
t=table(scenarios,labels,seeds,sec,synthesis_s,power_s,crossfit_s,surrogate_ppc_s,domain_s,peak_memory_bytes,status);
writetable(t,fullfile(outDir,"benchmark.tsv"),'FileType','text','Delimiter','\t');
phaseMean=mean(sec(scenarios=="phase"));scenario_worlds=[400;600;600;4200];scenario_sec=[sec(1);sec(2);sec(3);phaseMean];
extrapolation=table(["null";"additive";"aperiodic";"phase"],scenario_worlds,scenario_sec,scenario_worlds.*scenario_sec/22, ...
    'VariableNames',["scenario" "planned_worlds" "benchmark_seconds_per_world" "estimated_wall_s_22workers"]);
bisection=3*64*sec(3)/22;writetable(extrapolation,fullfile(outDir,"benchmark_extrapolation.tsv"),'FileType','text','Delimiter','\t');
outputs=struct('benchmark_dir',string(outDir),'table',t,'extrapolation',extrapolation, ...
    'estimated_main_wall_s_22workers',sum(extrapolation.estimated_wall_s_22workers), ...
    'estimated_bisection_wall_s_22workers',bisection,'estimated_total_wall_s_22workers',sum(extrapolation.estimated_wall_s_22workers)+bisection);
end

function outputs = runBenchmarkOptimized(root)
% The optimized benchmark is separate from the historical directory and
% intentionally does not create or modify a parallel pool.
addpath(fullfile(root,"scripts","ieeg","acc01"));
c=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
scenarios=["null";"additive";"aperiodic";"phase";"phase"];
labels=["null";"additive@Q75";"aperiodic@Q75_lambda_provisional"; ...
    "phase@Q75_kappa1_provisional";"phase@Q75_kappa4_provisional"];
q75=c.additive.a_erp_grid_uV(3);lambda=c.aperiodic.benchmark_lambda_provisional_uV;
effects={struct();struct('a_erp_uV',q75); ...
    struct('a_erp_uV',q75,'lambda_a_uV',lambda,'provisional',true); ...
    struct('a_erp_uV',q75,'kappa',1,'provisional',true); ...
    struct('a_erp_uV',q75,'kappa',4,'provisional',true)};
seeds=c.randomization.base_seed+[100001;203001;303001;403001;403002];
outDir=fullfile(root,"processed","subject","group","ieeg_acc00_sim", ...
    "BangYoureDead","checkpoints","benchmark_optimized_20260829");
if ~isfolder(outDir),mkdir(outDir);end
sec=nan(5,1);synthesis_s=sec;power_s=sec;crossfit_s=sec;
surrogate_ppc_s=sec;domain_s=sec;peak_memory_bytes=sec;
status=repmat("running",5,1);
for i=1:5
    worldClock=tic;
    w=acc00_sim_generate_world(root,c,"benchmark",seeds(i),scenarios(i),effects{i});
    r=acc01_run_metrics_optimized(w,c);
    sec(i)=toc(worldClock);synthesis_s(i)=w.timing.synthesis_s;
    power_s(i)=r.timing.power_s;crossfit_s(i)=r.timing.crossfit_s;
    surrogate_ppc_s(i)=r.timing.surrogate_ppc_s;domain_s(i)=r.timing.domain_s;
    peak_memory_bytes(i)=w.timing.peak_memory_bytes;status(i)="complete";
    checkpoint=struct('scenario',scenarios(i),'label',labels(i),'seed',seeds(i), ...
        'truth',w.truth,'contract_sha256',w.contract_sha256,'elapsed_s',sec(i), ...
        'synthesis_s',synthesis_s(i),'power_s',power_s(i), ...
        'crossfit_s',crossfit_s(i),'surrogate_ppc_s',surrogate_ppc_s(i), ...
        'domain_s',domain_s(i),'peak_memory_bytes',peak_memory_bytes(i), ...
        'status',status(i),'result',r);
    tmp=fullfile(outDir,sprintf("world_%02d.tmp.mat",i));
    final=fullfile(outDir,sprintf("world_%02d.mat",i));
    save(tmp,"checkpoint","-v7.3");movefile(tmp,final,"f");
end
t=table(scenarios,labels,seeds,sec,synthesis_s,power_s,crossfit_s, ...
    surrogate_ppc_s,domain_s,peak_memory_bytes,status);
writetable(t,fullfile(outDir,"benchmark.tsv"), ...
    'FileType','text','Delimiter','\t');
outputs=struct('benchmark_dir',string(outDir),'table',t, ...
    'status',"complete",'parpool_started',false);
end
end

function outputs = runSmoke(root)
addpath(fullfile(root,"scripts","ieeg","acc01"));
c=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
cs=c;cs.surrogate.max_surrogates_per_estimable_pair=4;cs.surrogate.enumerate_when_quantized_feasible_count_lte=4;
names=["null" "additive" "aperiodic" "phase"];q75=c.additive.a_erp_grid_uV(3);
effects={struct();struct('a_erp_uV',q75);struct('a_erp_uV',q75,'lambda_a_uV',c.aperiodic.benchmark_lambda_provisional_uV,'provisional',true);struct('a_erp_uV',q75,'kappa',1,'provisional',true)};
results=cell(numel(names),1);worlds=cell(numel(names),1);
for i=1:numel(names)
    worlds{i}=acc00_sim_generate_world(root,c,"smoke",c.randomization.base_seed+1,names(i),effects{i});results{i}=acc01_run_metrics(worlds{i},cs);
    assert(all(isfinite([results{i}.pair_rows.power]))&&all(isfinite([results{i}.pair_rows.phase])),"Smoke metrics must be finite.");
    assert(height(results{i}.subject_metrics.table)==2&&numel(unique(results{i}.subject_metrics.table.subject))==2,"Smoke subject keys must be unique.");
    sf=results{i}.signflip.R1;assert(isequal(sf.observed,sf.statistics(1,:))&&all(sf.signs(1,:)==1),"Smoke observed row identity failed.");
    assert(isfield(worlds{i},"truth")&&worlds{i}.truth.scenario==names(i),"Smoke truth metadata missing.");
    assert(all([results{i}.pair_rows.refit_surrogate_count]==min(32,[results{i}.pair_rows.surrogate_count]))&& ...
        all(string({results{i}.pair_rows.mphase_surrogate_reducer})=="mean"),"Smoke Mphase refit contract failed.");
end
% A separate 60-event fixture exercises the fold operator with an estimable
% design; the 12-label end-to-end fixture is intentionally minimal.
referenceEvents=readtable(c.inputs.label_covariates.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
referenceEvents=referenceEvents(1:60,:);referenceEvents.label_time_s=linspace(1,19,60)';
designCols=endsWith(referenceEvents.Properties.VariableNames,"_delta_z")|endsWith(referenceEvents.Properties.VariableNames,"movie_time_linear_z")|endsWith(referenceEvents.Properties.VariableNames,"movie_time_quadratic_z");
fixtureDesign=cos(2*pi*(0:59)'*(1:12)/60);referenceEvents.z_e=fixtureDesign(:,1);referenceEvents{:,designCols}=fixtureDesign(:,2:12);
[~,crossDiag]=acc01_crossfit_erp_subtract(worlds{1}.pairs(1),referenceEvents,cs,true);
assert(crossDiag.reference_relative_error<1e-12,"ACC00-SIM optimized cross-fit differs from pointwise reference.");
assert(any(abs(worlds{2}.pairs(1).signal-worlds{1}.pairs(1).signal)>0));assert(worlds{3}.truth.lambda_a_uV>0);assert(worlds{4}.truth.kappa==1);
v=validate_acc00_simulation('ProjectRoot',root);outDir=fullfile(root,"processed","subject","group","ieeg_acc00_sim","BangYoureDead","checkpoints","smoke");if ~isfolder(outDir),mkdir(outDir),end
save(fullfile(outDir,"smoke_results.mat"),"results","v","names","crossDiag","-v7.3");outputs=struct('status',"PASS",'validation',v,'checkpoint',string(outDir),'crossfit_reference_relative_error',crossDiag.reference_relative_error);fprintf("ACC00-SIM smoke: PASS\n");
end

function outputs = runSmokeOptimized(root)
% Serial optimized smoke; this stage intentionally never calls gcp/parpool.
addpath(fullfile(root,"scripts","ieeg","acc01"));
c=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
cs=c;cs.surrogate.max_surrogates_per_estimable_pair=4;
cs.surrogate.enumerate_when_quantized_feasible_count_lte=4;
names=["null" "additive" "aperiodic" "phase"];q75=c.additive.a_erp_grid_uV(3);
effects={struct();struct('a_erp_uV',q75); ...
    struct('a_erp_uV',q75,'lambda_a_uV',c.aperiodic.benchmark_lambda_provisional_uV,'provisional',true); ...
    struct('a_erp_uV',q75,'kappa',1,'provisional',true)};
results=cell(numel(names),1);
for i=1:numel(names)
    w=acc00_sim_generate_world(root,c,"smoke",c.randomization.base_seed+1,names(i),effects{i});
    results{i}=acc01_run_metrics_optimized(w,cs);
    assert(all(isfinite([results{i}.pair_rows.power]))&& ...
        all(isfinite([results{i}.pair_rows.phase])),"Optimized smoke metrics must be finite.");
    assert(height(results{i}.subject_metrics.table)==2&& ...
        numel(unique(results{i}.subject_metrics.table.subject))==2, ...
        "Optimized smoke subject keys must be unique.");
    sf=results{i}.signflip.R1;
    assert(isequal(sf.observed,sf.statistics(1,:))&&all(sf.signs(1,:)==1), ...
        "Optimized smoke observed row identity failed.");
    assert(all([results{i}.pair_rows.refit_surrogate_count]== ...
        min(32,[results{i}.pair_rows.surrogate_count]))&& ...
        all(string({results{i}.pair_rows.mphase_surrogate_reducer})=="mean"), ...
        "Optimized smoke Mphase refit contract failed.");
end
outDir=fullfile(root,"processed","subject","group","ieeg_acc00_sim", ...
    "BangYoureDead","checkpoints","edgeguard_fix_20260829","smoke_optimized");
if ~isfolder(outDir),mkdir(outDir);end
v=validate_acc00_simulation('ProjectRoot',root,'OutputDir',fullfile(outDir,"validator"));
save(fullfile(outDir,"smoke_results.mat"),"results","v","names","-v7.3");
outputs=struct('status',"PASS",'validation',v,'checkpoint',string(outDir), ...
    'parpool_started',false);
fprintf("ACC00-SIM optimized smoke: PASS\n");
end

function outputs = writeContract(root)
% Write only values which are mechanically rederived from the frozen sources.
% The validator independently repeats these calculations; this function never
% accepts a hand-entered ERP amplitude level.
outPath = fullfile(root, "scripts", "ieeg", "acc00_sim_contract.json");
covPath = fullfile(root, "processed", "subject", "group", "ieeg_acc00_sim", ...
    "BangYoureDead", "task-bangyouredead_desc-acc00sim-label-covariates.tsv");
i02Path = fullfile(root, "processed", "subject", "group", "ieeg_i02_acc", ...
    "BangYoureDead", "task-bangyouredead_desc-i02acc-subject-erp.tsv");
inferencePath = fullfile(root, "processed", "subject", "group", "ieeg_i02_acc", ...
    "BangYoureDead", "task-bangyouredead_desc-i02acc-erp-inference.tsv");
labelPath = fullfile(root, "processed", "subject", "group", "ieeg_acc00_vjepa2_preflight", ...
    "BangYoureDead", "task-bangyouredead_desc-acc00vjepa2-prediction-difficulty.tsv");
pairPath = fullfile(root, "processed", "subject", "group", "ieeg_i00", ...
    "BangYoureDead", "task-bangyouredead_desc-i00-macro-bipolar-pairs.tsv");
channelQcPath = fullfile(root, "processed", "subject", "group", "ieeg_i00", ...
    "BangYoureDead", "task-bangyouredead_desc-i00-channel-qc.tsv");
badPath = fullfile(root, "processed", "subject", "group", "ieeg_i00", ...
    "BangYoureDead", "task-bangyouredead_desc-i00-bad-segments.tsv");
timingPath = fullfile(root, "processed", "subject", "group", "ieeg_i00", ...
    "BangYoureDead", "task-bangyouredead_desc-i00-timing-transforms.tsv");
requireFiles([covPath; i02Path; inferencePath; labelPath; pairPath; channelQcPath; badPath; timingPath]);

anchor = recomputeI02Anchor(i02Path, inferencePath);
backgroundAnchor = recomputeBackgroundAnchor(pairPath, channelQcPath);
labels = readtable(labelPath, "FileType", "text", "Delimiter", "\t", "VariableNamingRule", "preserve");
z = double(labels.prediction_difficulty_z);
b = 0.3055533803;
energy = exp(2 * b * z);
covZE = cov(z, energy, 1);
betaFactor = covZE(1, 2) ./ var(z, 1);
g = exp(b * z);

contract = struct();
contract.analysis_id = "ACC00-SIM-v1.0.0";
contract.status = "stage1_implementation_only_full_simulation_not_authorized";
contract.no_acc_neural_results_read = true;
contract.sampling = struct("synthetic_hz", 1000, "analysis_hz", 250, ...
    "anti_alias_lowpass_hz", 100, "anti_alias_order", 3);
contract.background = struct("recipe", "random_phase_frequency_domain_A(f)=1/max(f,0.5Hz)_DC_zero_Hermitian_IFFT", ...
    "power_exponent_chi", 2, "f_min_hz", 0.5, "b0_uV", backgroundAnchor.b0_uV, ...
    "b0_definition", "median over 351 eligible ACC pairs of sqrt(median_rms_uv(contact_a)^2+median_rms_uv(contact_b)^2)", ...
    "log_scale_sigma", struct("subject",0.25,"session",0.15,"pair",0.35), ...
    "normalization_domain", "background-only full 250-Hz trace RMS after anti-alias and decimation", ...
    "scales_injection", false);
contract.windows_s = struct("power_pre", [-0.50 -0.10], "power_post", [0.10 0.50], ...
    "phase_response", [0.10 0.50], "filter_support", [-0.75 0.75]);
contract.power = struct("band_hz", [4 8], "theta_periodogram_bins_hz", [5 7.5], ...
    "hann_epoch_s", 0.40, "periodogram_bins_hz", 2.5:2.5:30, ...
    "aperiodic_fit_bins_hz", [2.5 10:2.5:30], "log_base", 10, "robust_fit", "robustfit bisquare", ...
    "robustfit_iteration_limit_warning_id", "stats:statrobustfit:IterationLimit", ...
    "robustfit_iteration_limit_policy", "suppress stdout and count each fit per world; estimator settings unchanged");
contract.phase = struct("band_hz", [4 8], "filter", "Butterworth", "filter_order", 3, ...
    "zero_phase", true, "analytic_signal", "Hilbert", "erp_crossfit_folds", 5, ...
    "fold_assignment", "contiguous_movie_time_blocks", "ppc", "pairwise_phase_consistency", ...
    "mphase_surrogate_reducer", "mean", ...
    "refit_surrogates", struct("J", 32, ...
        "selection", "uniform_without_replacement_from_existing_pair_surrogates_same_stream", ...
        "when_B_lt_J", "use_all"), ...
    "filter_scope", "per_continuous_clean_block", ...
    "corr_series_diagnostic_surrogates", 1024);
contract.surrogate = struct("type", "whole_schedule_circular_time_shift", ...
    "signal_moved", false, "event_attributes_carried", true, "delta_grid_hz", 250, ...
    "delta_quantum_s", 0.004, "max_surrogates_per_estimable_pair", 1024, ...
    "enumerate_when_quantized_feasible_count_lte", 1024, ...
    "empty_domain", "not_estimable_no_substitution", ...
    "domain_guard_bad_dilation_s", 0.75, ...
    "domain_guard_edge_lower_s", 1.5, ...
    "domain_guard_edge_upper_s", 1.5, ...
    "domain_guard_edge_basis", "closed intervals [0,1.5] and [T-1.5,T] on the 250-Hz grid; 1.5 s * fs is an integer sample count, so the existing closed-grid quantization requires no rounding adjustment", ...
    "feasible_domain", "[0,T) minus union_e(((edge_lower_guarded_by_1.5s or edge_upper_guarded_by_1.5s or bad_dilated_by_0.75s) - t_e) mod T)");
contract.additive = struct("kernel", "half_sine_rms1", "kernel_formula", ...
    "h(t)=sqrt(2)*sin(pi*(t-0.10)/0.40), t in [0.10,0.50]s, zero outside, RMS=1", ...
    "kernel_support_s", [0.10 0.50], "kernel_rms", 1, ...
    "amplitude_relation", "amplitude_e_uV=A_ERP_uV*exp(a+b*z_e)*pair_polarity", ...
    "gain_intercept_a", 0, "gain_slope_b", b, ...
    "a_erp_grid_uV", anchor.a_erp_grid_uV, "a_erp_grid_definition", ...
    "[0,prctile(max(D_i,0),p,Method=""inclusive"") for p=50,75,90]; D_i=RMS(hard_cut,0.04:0.40)-mean(RMS(shot_matched_pseudo),RMS(circular_shift)) for I02 R1 bipolar subject waveform", ...
    "percentile_method", "MATLAB prctile Method=""inclusive"" (linear interpolation, position 1+(n-1)*p/100; never version default)", ...
    "d_i_uV", anchor.d_i_V .* 1e6, ...
    "anchor_mean_uV", anchor.mean_uV, "anchor_median_uV", anchor.median_uV, ...
    "beta_gen_uV2_per_z", (anchor.a_erp_grid_uV .^ 2) * betaFactor);
contract.label_gain_diagnostics = struct("z_quantiles", prctile(z, [10 50 90]), ...
    "z_range", [min(z) max(z)], "gain_quantiles", exp(b * prctile(z, [10 50 90])), ...
    "gain_range", [min(g) max(g)], "gain_q90_over_q10", exp(b * (prctile(z,90)-prctile(z,10))), ...
    "cov_z_exp2bz_over_var_z", betaFactor, ...
    "top_1pct_event_energy_share", sum(maxk(energy, ceil(0.01*numel(energy)))) / sum(energy));
contract.aperiodic_calibration = struct("method", "monotone_bisection_common_random_numbers", ...
    "calibration_worlds", 64, "target", "match additive median direct 4-8Hz post-window log10 power increment", ...
    "residual_tolerance_dB", 0.001, "residual_tolerance_rule", "absolute_residual_dB <= min(0.001, MCSE_dB/10)", ...
    "bracket_relative_width_tolerance", 1e-6, "max_iterations", 32, ...
    "on_max_iterations", "record_lambda_residual_mcse_iterations_and_not_converged", ...
    "initial_lambda_guess", "benchmark_lambda_provisional_uV * A_level_uV / A_Q75plus_uV", ...
    "initial_bracket_relative_to_guess", [1/8 8], ...
    "bracket_expansion", struct("factor",8,"max_steps_per_side",4, ...
        "direction","expand_only_endpoint_indicated_by_same_sign_residual", ...
        "record_each_step",true,"on_failure","record_not_converged_without_silent_fallback"), ...
    "bisection_space", "log10_lambda_geometric_midpoint", ...
    "calibration_scenario_index", 3, "calibration_seed_index", [1001 1064], ...
    "calibration_seed_formula", "base_seed + 100000*3 + 1000*effect_index + seed_index; seed_index=1001:1064", ...
    "crn", "within calibration world, synthesize and cache once; additive target and every lambda evaluation share bit-identical background and event random variates", ...
    "increment", struct("event_formula", "mean(theta_periodogram_bins_log10_power(post+injection)) - mean(theta_periodogram_bins_log10_power(post|base))", ...
        "pair_aggregation", "mean over the same eligible power-event table passed to acc01_compute_power_metric", ...
        "world_aggregation", "equal-weight mean over 351 eligible pairs", ...
        "signal_path", "paired raw_pair_signal post window [0.10,0.50] s; demeaned Hann epoch; log10(abs(fft).^2/sum(hann.^2)); theta bins only", ...
        "stops_before", "robustfit_aperiodic_residual", ...
        "excluded", ["pre_window" "crossfit" "per_block_filtering" "robustfit_aperiodic_adjustment" "covariate_regression" "surrogate_reference" "inference"]), ...
    "residual", "median_over_64_CRN_paired_world_differences_aperiodic_direct_minus_additive_direct_log10_units; multiply_by_10_for_dB", ...
    "mcse", struct("method","world_level_bootstrap_paired_median","bootstrap_replicates",10000, ...
        "rng","RandStream('mt19937ar','Seed',20260827+777000)","seed",21037827,"statistic","sd_of_bootstrap_medians"), ...
    "parallel_workers",8,"parallel_workers_rationale","calibration cache loads peak synchronously at about 4 GiB per worker; 8 workers keep the calibration peak below host RAM alongside MATLAB root, rocketstation, and OS headroom");
contract.aperiodic = struct("event_recipe", "independent random-phase 1/f segment times h(t)", ...
    "normalization", "each event 0.10-0.50s product normalized to RMS=1 before lambda_A*exp(a+b*z_e)", ...
    "polarity", "none", "lambda_unit", "window uV RMS", ...
    "benchmark_lambda_provisional_uV", anchor.mean_uV, "benchmark_value_status", "provisional_timing_only");
contract.design_matrix = struct("formula", "1+z_e+9_delta_z+movie_time_linear_z+movie_time_quadratic_z", ...
    "health_diagnostics", ["condition_number_2norm" "correlation_z_with_each_covariate" "source_complete_coverage" "power_phase_common_support_coverage"]);
contract.randomization = struct("base_seed", 20260827, "world_seed", "base_seed + 100000*scenario_index + 1000*effect_index + seed_index", ...
    "worker_rng", "independent deterministic RandStream per world", "common_random_numbers", "same background/event random variates within aperiodic calibration bracket");
contract.world_plan = struct("world_count_upper_bound",5800, "run_is_not_world_axis",true, ...
    "null_seeds",400,"additive_positive_a_erp_seeds",200,"aperiodic_positive_a_erp_seeds",200, ...
    "phase_per_a_erp_kappa_seeds",200,"phase_kappa_target_null_sd_multipliers",[0 1 2 4 8 16 32], ...
    "phase_max_grid_points",7,"phase_r2_hard_cap",0.64,"phase_kappa_hard_cap_approx",3.7, ...
    "null_pool_role","zero-level reference and sigma_null(Mphase) estimation; no duplicated zero worlds", ...
    "aperiodic_bisection_crn_worlds",64,"bisection_worlds_excluded_from_main_upper_bound",true);
contract.sign_flip = struct("runs", ["R1" "R2"], "subject_counts", [16 13], ...
    "enumeration_counts", [65536 8192], "subject_order", "lexicographic subject key", ...
    "metric_order", ["Mpower" "Mphase"], "observed_sign_row", 1, ...
    "observed_row", "all_plus_one row of same sign matrix and same matrix multiplication", ...
    "p_value", "exhaustive two-sided max-stat with enumerated rows as denominator; no floating tolerance", ...
    "nan_policy", "complete subject rows within run and metric family; report exclusions", ...
    "studentization", struct("enabled", true, ...
        "scale", "per_metric_flip_distribution_sd", ...
        "apply_before", "max_stat", ...
        "per_metric_two_sided_uses", "raw_statistics"));
contract.inputs = makeInputManifest(["label" "label_covariates" "i02_subject_erp" "i02_erp_inference" "i00_acc_pairs" "i00_channel_qc" "i00_bad_segments" "i00_timing"], ...
    [labelPath covPath i02Path inferencePath pairPath channelQcPath badPath timingPath]);
contract.inputs.i00_channel_qc.role = "unit-anchor only";
contract.implementation_boundary = struct("full_simulation", "authorized_route_B_mpower_only_additive_effect_index_2_to_4_seed_index_1_to_100; aperiodic_effect_index_2_previous_negative_target_withdrawal_superseded_for_direct_scale_calibration_only; no_aperiodic_main_worlds_authorized", ...
    "checkpoint_directory", "processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints", ...
    "parallel_workers", 12, "parallel_workers_rationale", ...
    "host RAM 61.65 GiB; benchmark observed 3.1-3.7 GiB peak per world; 22 workers imply about 80 GiB and exceed RAM; 12 workers imply about 44 GiB and retain OS/agent headroom", ...
    "long_run", "detached_process_with_completion_notification");

fid = fopen(outPath, "w");
if fid < 0, error("Bmovie:ACC00SIM:ContractWrite", "Cannot write %s", outPath); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", jsonencode(contract, "PrettyPrint", true));
clear cleanup
outputs = struct("contract_path", string(outPath), "sha256", sha256File(outPath), "a_erp_grid_uV", anchor.a_erp_grid_uV);
hashPath=fullfile(root,"scripts","ieeg","acc00_sim_contract.sha256");hf=fopen(hashPath,"w");assert(hf>=0,"Cannot write %s",hashPath);hc=onCleanup(@()fclose(hf));fprintf(hf,"%s  acc00_sim_contract.json\n",outputs.sha256);clear hc
fprintf("ACC00-SIM contract: SHA-256 %s\n", outputs.sha256);
end

function anchor = recomputeBackgroundAnchor(pairPath, channelQcPath)
p=readtable(pairPath,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
p=p(logical(p.qc_eligible)&contains(p.location,"ACC"),:);q=readtable(channelQcPath,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
v=nan(height(p),1);
for i=1:height(p)
    a=q(q.session_id==p.session_id(i)&q.series=="LFP_macro"&q.origchannel_name==p.contact_a(i),:);
    b=q(q.session_id==p.session_id(i)&q.series=="LFP_macro"&q.origchannel_name==p.contact_b(i),:);
    assert(height(a)==1&&height(b)==1,"ACC00-SIM B0 channel-QC join failed for %s %s.",p.session_id(i),p.pair_name(i));
    v(i)=hypot(double(a.median_rms_uv),double(b.median_rms_uv));
end
assert(height(p)==351&&all(isfinite(v)),"ACC00-SIM B0 requires 351 finite eligible ACC pairs.");
anchor=struct("b0_uV",median(v),"pair_anchor_uV",v);
fprintf("ACC00-SIM background B0 (uV): %.12g\n",anchor.b0_uV);
end

function anchor = recomputeI02Anchor(i02Path, inferencePath)
t = readtable(i02Path, "FileType", "text", "Delimiter", "\t", "TextType", "string", "VariableNamingRule", "preserve");
required = ["subject" "run" "reference_type" "condition" "time_s" "erp_value"];
requireColumns(t, required);
subjects = unique(t.subject(t.run == "R1" & t.reference_type == "bipolar"), "stable");
d = nan(numel(subjects), 1);
for i = 1:numel(subjects)
    rows = t.subject == subjects(i) & t.run == "R1" & t.reference_type == "bipolar" & t.time_s >= 0.04 & t.time_s <= 0.40;
    rmsFor = @(condition) sqrt(mean(double(t.erp_value(rows & t.condition == condition)).^2));
    hard = rmsFor("hard_cut");
    d(i) = hard - mean([rmsFor("shot_matched_pseudo"), rmsFor("circular_shift")]);
end
if any(~isfinite(d)) || numel(d) ~= 16
    error("Bmovie:ACC00SIM:I02Anchor", "I02 anchor requires 16 finite R1 bipolar D_i values.");
end
inf = readtable(inferencePath, "FileType", "text", "Delimiter", "\t", "TextType", "string", "VariableNamingRule", "preserve");
r1 = inf(inf.run == "R1", :);
meanU = mean(d) * 1e6; medianU = median(d) * 1e6;
if height(r1) ~= 1 || abs(mean(d) - double(r1.mean_hard_minus_control_erp_rms)) > 1e-18 || ...
        abs(median(d) - double(r1.median_hard_minus_control_erp_rms)) > 1e-18
    error("Bmovie:ACC00SIM:I02AnchorMismatch", "Recomputed I02 R1 mean/median do not match erp-inference.tsv.");
end
anchor = struct("d_i_V", d, "a_erp_grid_uV", [0 prctile(max(d, 0), [50 75 90], Method="inclusive") * 1e6], ...
    "mean_uV", meanU, "median_uV", medianU);
fprintf("ACC00-SIM I02 anchor D_i (uV, sorted): %s\\n", sprintf("%.12g ", sort(d * 1e6)));
fprintf("ACC00-SIM I02 anchor prctile grid (uV): %s\\n", sprintf("%.12g ", anchor.a_erp_grid_uV));
end

function manifest = makeInputManifest(roles, paths)
manifest = struct();
for i = 1:numel(roles)
    manifest.(char(roles(i))) = struct("path", char(paths(i)), "sha256", char(sha256File(paths(i))));
end
end

function outputs = buildCovariates(root)
stimPath = fullfile(root, "processed", "stim", "BangYoureDead", "low_level_features.csv");
audioPath = fullfile(root, "processed", "subject", "group", "ieeg_i02_acc", ...
    "BangYoureDead", "task-bangyouredead_desc-i02acc-audio-envelope.tsv");
labelPath = fullfile(root, "processed", "subject", "group", "ieeg_acc00_vjepa2_preflight", ...
    "BangYoureDead", "task-bangyouredead_desc-acc00vjepa2-prediction-difficulty.tsv");
requireFiles([stimPath; audioPath; labelPath]);

visual = readtable(stimPath, "TextType", "string", "VariableNamingRule", "preserve");
audio = readtable(audioPath, "FileType", "text", "Delimiter", "\t", ...
    "TextType", "string", "VariableNamingRule", "preserve");
labels = readtable(labelPath, "FileType", "text", "Delimiter", "\t", ...
    "TextType", "string", "VariableNamingRule", "preserve");
visualColumns = ["luminance", "rms_contrast", "global_motion", "SF_low", ...
    "SF_mid", "SF_high", "edge_density", "novelty"];
requireColumns(visual, ["timestamp_seconds", visualColumns]);
requireColumns(audio, ["movie_time_s", "audio_rms_envelope"]);
requireColumns(labels, ["label_index", "label_time_s", "prediction_difficulty_z"]);

n = height(labels);
labelIndex = double(labels.label_index);
labelTime = double(labels.label_time_s);
z = double(labels.prediction_difficulty_z);
preStart = labelTime - 0.50; preStop = labelTime - 0.10;
postStart = labelTime + 0.10; postStop = labelTime + 0.50;
sourceComplete = true(n, 1);
out = table(labelIndex, labelTime, z, 'VariableNames', ["label_index", "label_time_s", "z_e"]);
for c = 1:numel(visualColumns)
    name = visualColumns(c);
    x = double(visual.timestamp_seconds); y = double(visual.(name));
    [pre, preOk] = pwlMeans(x, y, preStart, preStop);
    [post, postOk] = pwlMeans(x, y, postStart, postStop);
    sourceComplete = sourceComplete & preOk & postOk;
    out.(name + "_pre") = pre;
    out.(name + "_post") = post;
    out.(name + "_delta") = post - pre;
end
[pre, preOk] = pwlMeans(double(audio.movie_time_s), double(audio.audio_rms_envelope), preStart, preStop);
[post, postOk] = pwlMeans(double(audio.movie_time_s), double(audio.audio_rms_envelope), postStart, postStop);
sourceComplete = sourceComplete & preOk & postOk;
out.audio_rms_envelope_pre = pre;
out.audio_rms_envelope_post = post;
out.audio_rms_envelope_delta = post - pre;
out.source_complete = sourceComplete;

deltaNames = [visualColumns + "_delta", "audio_rms_envelope_delta"];
timeCentered = labelTime - mean(labelTime(sourceComplete));
out.movie_time_linear = timeCentered;
out.movie_time_quadratic = timeCentered .^ 2;
standardized = [deltaNames, "movie_time_linear", "movie_time_quadratic"];
for name = standardized
    x = double(out.(name));
    mu = mean(x(sourceComplete)); sigma = std(x(sourceComplete), 0);
    if ~isfinite(sigma) || sigma == 0
        error("Bmovie:ACC00SIM:ZeroCovariateScale", "Covariate %s has invalid scale.", name);
    end
    y = nan(n, 1); y(sourceComplete) = (x(sourceComplete) - mu) ./ sigma;
    out.(name + "_z") = y;
end

outDir = fullfile(root, "processed", "subject", "group", "ieeg_acc00_sim", "BangYoureDead");
if ~isfolder(outDir), mkdir(outDir); end
outPath = fullfile(outDir, "task-bangyouredead_desc-acc00sim-label-covariates.tsv");
writetable(out, outPath, "FileType", "text", "Delimiter", "\t");
hash = sha256File(outPath);
manifest = table(["label"; "visual_continuous"; "audio_continuous"; "label_covariates"], ...
    string([labelPath; stimPath; audioPath; outPath]), ...
    [sha256File(labelPath); sha256File(stimPath); sha256File(audioPath); hash], ...
    'VariableNames', ["role", "path", "sha256"]);
writetable(manifest, fullfile(outDir, "task-bangyouredead_desc-acc00sim-covariate-manifest.tsv"), ...
    "FileType", "text", "Delimiter", "\t");
outputs = struct("output_path", string(outPath), "sha256", hash, ...
    "n_rows", height(out), "n_source_complete", nnz(sourceComplete), ...
    "n_source_incomplete", nnz(~sourceComplete));
fprintf("ACC00-SIM covariates: %d rows (%d complete), SHA-256 %s\n", ...
    outputs.n_rows, outputs.n_source_complete, outputs.sha256);
end

function [means, okay] = pwlMeans(time, values, starts, stops)
time = double(time(:)); values = double(values(:));
if any(diff(time) <= 0) || any(~isfinite(values))
    error("Bmovie:ACC00SIM:InvalidContinuousSource", "Continuous source must have unique ordered finite values.");
end
means = nan(numel(starts), 1); okay = starts >= time(1) & stops <= time(end);
for k = find(okay)'
    knots = [starts(k); time(time > starts(k) & time < stops(k)); stops(k)];
    samples = interp1(time, values, knots, "linear");
    means(k) = trapz(knots, samples) ./ (stops(k) - starts(k));
end
end

function requireColumns(t, names)
present = string(t.Properties.VariableNames);
missing = names(~ismember(names, present));
if ~isempty(missing)
    error("Bmovie:ACC00SIM:MissingColumn", "Missing required column(s): %s", strjoin(missing, ", "));
end
end

function requireFiles(paths)
for p = string(paths(:))'
    if ~isfile(p), error("Bmovie:ACC00SIM:MissingFile", "Missing file: %s", p); end
end
end

function hash = sha256File(path)
fid = fopen(path, "r");
if fid < 0, error("Bmovie:ACC00SIM:OpenFailed", "Cannot open %s", path); end
cleaner = onCleanup(@() fclose(fid));
digest = java.security.MessageDigest.getInstance("SHA-256");
while true
    bytes = fread(fid, 1024 * 1024, "*uint8");
    if isempty(bytes), break; end
    digest.update(typecast(bytes, "int8"));
end
hash = upper(string(reshape(dec2hex(typecast(digest.digest(), "uint8"), 2)', 1, [])));
end
