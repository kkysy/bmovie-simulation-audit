function report = validate_acc00_simulation(varargin)
p=inputParser;
p.addParameter('ProjectRoot',pwd);
p.addParameter('OutputDir',"",@(x)ischar(x)||isstring(x));
p.addParameter('RunOptimizationEquivalence',false,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});root=string(p.Results.ProjectRoot);
addpath(fullfile(root,"scripts","ieeg","acc01"));contractPath=fullfile(root,"scripts","ieeg","acc00_sim_contract.json");c=jsondecode(fileread(contractPath));
names=strings(0,1);status=strings(0,1);detail=strings(0,1);
add("contract kernel name",c.additive.kernel=="half_sine_rms1",c.additive.kernel);
add("contract kernel formula",contains(c.additive.kernel_formula,"sqrt(2)*sin"),c.additive.kernel_formula);
fs=c.sampling.synthetic_hz;t=(.10:1/fs:.50-1/fs)';h=sqrt(2)*sin(pi*(t-.10)/.40);
add("kernel discrete RMS",abs(rms(h)-1)<1e-12,sprintf("rms=%.17g",rms(h)));
outsideT=[0 .099 .501 .8]';hout=zeros(size(outsideT));inside=outsideT>=.10&outsideT<=.50;hout(inside)=sqrt(2)*sin(pi*(outsideT(inside)-.10)/.40);
add("kernel outside support zero",all(hout==0),mat2str(hout'));

contractHash=sha256File(contractPath);sidecar=split(strtrim(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.sha256"))));
add("contract SHA sidecar",upper(string(sidecar(1)))==contractHash,contractHash);
roles=string(fieldnames(c.inputs));
for role=roles'
    z=c.inputs.(role);add("input hash "+role,sha256File(z.path)==upper(string(z.sha256)),string(z.path));
end
add("channel-QC unit-anchor role",c.inputs.i00_channel_qc.role=="unit-anchor only",c.inputs.i00_channel_qc.role);
add("Mphase reducer",c.phase.mphase_surrogate_reducer=="mean",c.phase.mphase_surrogate_reducer);
add("Mphase refit surrogate J",c.phase.refit_surrogates.J==32,sprintf("%d",c.phase.refit_surrogates.J));
add("family studentization contract",c.sign_flip.studentization.enabled&& ...
    c.sign_flip.studentization.scale=="per_metric_flip_distribution_sd"&& ...
    c.sign_flip.studentization.apply_before=="max_stat",c.sign_flip.studentization.scale);
add("per-block filter contract",c.phase.filter_scope=="per_continuous_clean_block",c.phase.filter_scope);

vfs=double(c.sampling.analysis_hz);
voff=round(double(c.windows_s.filter_support(1))*vfs):round(double(c.windows_s.filter_support(2))*vfs);
add("domain guard fields",c.surrogate.domain_guard_bad_dilation_s==0.75&& ...
    c.surrogate.domain_guard_edge_lower_s==1.5&&c.surrogate.domain_guard_edge_upper_s==1.5&& ...
    c.surrogate.domain_guard_edge_lower_s*vfs==round(c.surrogate.domain_guard_edge_lower_s*vfs)&& ...
    c.surrogate.domain_guard_edge_upper_s*vfs==round(c.surrogate.domain_guard_edge_upper_s*vfs), ...
    sprintf("bad=%.3f edge=[%.3f,%.3f]",c.surrogate.domain_guard_bad_dilation_s, ...
    c.surrogate.domain_guard_edge_lower_s,c.surrogate.domain_guard_edge_upper_s));

% Sample-safety regression for the feasible domain (stage-3 retry crash):
% with a worst-case session length (nSig=ceil(T*fs)), every feasible grid
% delta must keep the cross-fit window round(t'*fs)+1+off inside the trace.
gp=struct('T',478.852206788966,'bad_segments_s',zeros(0,2));
gev=table((50:37:430)','VariableNames',"label_time_s");
gc2=c;gc2.surrogate.max_surrogates_per_estimable_pair=1e9;gc2.surrogate.enumerate_when_quantized_feasible_count_lte=1e9;
gs=acc01_build_empirical_surrogates(gp,gev,gc2,RandStream('Threefry','Seed',1));
gnSig=ceil(gp.T*vfs);gcen=round(mod(double(gev.label_time_s)+double(gs.delta_s)',gp.T)*vfs)+1;
gSafe=~isempty(gs.delta_s)&&all(gcen+voff(1)>=1&gcen+voff(end)<=gnSig,'all');
add("domain guard keeps support window in range",gSafe,sprintf("%d deltas checked",numel(gs.delta_s)));

src=readtable(c.inputs.i02_subject_erp.path,'FileType','text','Delimiter','\t','TextType','string','VariableNamingRule','preserve');
ss=unique(src.subject(src.run=="R1"&src.reference_type=="bipolar"),'stable');d=nan(numel(ss),1);
for i=1:numel(ss),r=src.subject==ss(i)&src.run=="R1"&src.reference_type=="bipolar"&src.time_s>=.04&src.time_s<=.40;q=@(z)sqrt(mean(double(src.erp_value(r&src.condition==z)).^2));d(i)=q("hard_cut")-mean([q("shot_matched_pseudo") q("circular_shift")]);end
grid=[0 prctile(max(d,0),[50 75 90],Method="inclusive")*1e6];
add("I02 anchor grid",max(abs(grid(:)-c.additive.a_erp_grid_uV(:)))<1e-12,mat2str(grid,12));
add("I02 D_i",max(abs(d*1e6-c.additive.d_i_uV(:)))<1e-12,sprintf("n=%d",numel(d)));
ct=readtable(c.inputs.label_covariates.path,'FileType','text','Delimiter','\t','TextType','string','VariableNamingRule','preserve');
add("covariate rows",height(ct)==1901&&nnz(ct.source_complete)==1900,sprintf("%d/%d",nnz(ct.source_complete),height(ct)));
z=double(ct.z_e);b=c.additive.gain_slope_b;ce=cov(z,exp(2*b*z),1);factor=ce(1,2)/var(z,1);
add("label gain",abs(factor-c.label_gain_diagnostics.cov_z_exp2bz_over_var_z)<1e-10,sprintf("%.12g",factor));

q=c.aperiodic_calibration;
add("direct calibration target definition",contains(string(q.target),"direct")&&contains(string(q.increment.event_formula),"post+injection")&&contains(string(q.increment.signal_path),"post window [0.10,0.50]")&&contains(string(q.residual),"aperiodic_direct_minus_additive_direct"),string(q.target));
calDir=fullfile(root,"processed","subject","group","ieeg_acc00_sim","BangYoureDead","checkpoints","recovery_mpower_only_20260829","calibration");
dw=readtable(fullfile(calDir,"additive_direct_theta_worlds.tsv"),"FileType","text","Delimiter","\t","TextType","string");
ds=readtable(fullfile(calDir,"additive_direct_theta_summary.tsv"),"FileType","text","Delimiter","\t","TextType","string");
add("D4 direct anchor world schema",height(dw)==192&&all(ismember(["effect_index" "seed_index" "seed" "direct_theta_increment_log10"],string(dw.Properties.VariableNames)))&&all(isfinite(dw.direct_theta_increment_log10)),sprintf("n=%d",height(dw)));
for effectIndex=2:4
    z=dw(dw.effect_index==effectIndex,:);r=ds(ds.effect_index==effectIndex,:);v=z.direct_theta_increment_log10;
    ok=height(z)==64&&height(r)==1&&abs(mean(v)-r.mean_direct_theta_increment_log10)<1e-15&&abs(median(v)-r.median_direct_theta_increment_log10)<1e-15&&abs(prctile(v,10,Method="inclusive")-r.q10_direct_theta_increment_log10)<1e-15&&abs(prctile(v,90,Method="inclusive")-r.q90_direct_theta_increment_log10)<1e-15;
    add("D4 direct anchor recompute effect "+effectIndex,ok,sprintf("n=%d mean=%.12g",height(z),mean(v)));
end

pairs=readtable(c.inputs.i00_acc_pairs.path,'FileType','text','Delimiter','\t','TextType','string','VariableNamingRule','preserve');
pairs=pairs(logical(pairs.qc_eligible)&contains(pairs.location,"ACC"),:);
channel=readtable(c.inputs.i00_channel_qc.path,'FileType','text','Delimiter','\t','TextType','string','VariableNamingRule','preserve');pairRms=nan(height(pairs),1);
for i=1:height(pairs)
    a=channel(channel.session_id==pairs.session_id(i)&channel.series=="LFP_macro"&channel.origchannel_name==pairs.contact_a(i),:);
    bb=channel(channel.session_id==pairs.session_id(i)&channel.series=="LFP_macro"&channel.origchannel_name==pairs.contact_b(i),:);
    if height(a)==1&&height(bb)==1,pairRms(i)=hypot(double(a.median_rms_uv),double(bb.median_rms_uv));end
end
b0=median(pairRms);add("background B0 independent",all(isfinite(pairRms))&&abs(b0-c.background.b0_uV)<1e-12,sprintf("%.12g",b0));
add("background frozen sigmas",isequal([c.background.log_scale_sigma.subject c.background.log_scale_sigma.session c.background.log_scale_sigma.pair],[.25 .15 .35]),"[.25 .15 .35]");
add("parallel worker count",c.implementation_boundary.parallel_workers==12,sprintf("%d",c.implementation_boundary.parallel_workers));
add("parallel RAM rationale",contains(c.implementation_boundary.parallel_workers_rationale,"61.65 GiB")&&contains(c.implementation_boundary.parallel_workers_rationale,"44 GiB"),c.implementation_boundary.parallel_workers_rationale);
add("robustfit warning diagnostic policy",c.power.robustfit_iteration_limit_warning_id=="stats:statrobustfit:IterationLimit"&&contains(c.power.robustfit_iteration_limit_policy,"count"),c.power.robustfit_iteration_limit_policy);
add("eligible pair count",height(pairs)==351,sprintf("%d",height(pairs)));
add("R1 subject family",numel(unique(pairs.subject(pairs.run=="R1")))==16,"16 subjects");
add("R2 subject family",numel(unique(pairs.subject(pairs.run=="R2")))==13,"13 subjects");
f1=acc01_exact_signflip_family(zeros(16,2));f2=acc01_exact_signflip_family(zeros(13,2));
add("R1 exact rows",size(f1.signs,1)==65536&&size(unique(f1.signs,"rows"),1)==65536,"65536");
add("R2 exact rows",size(f2.signs,1)==8192&&size(unique(f2.signs,"rows"),1)==8192,"8192");
add("observed matrix path R1",isequal(f1.observed,f1.statistics(1,:))&&all(f1.signs(1,:)==1),"row 1");
add("observed matrix path R2",isequal(f2.observed,f2.statistics(1,:))&&all(f2.signs(1,:)==1),"row 1");
add("individual signflip p-values",isequal(size(f1.p_two_sided),[1 2])&&all(f1.p_two_sided>=0&f1.p_two_sided<=1),mat2str(f1.p_two_sided));

% Independent family studentization check: per-metric two-sided values remain
% on the raw scale, while only the max-stat family uses flip-distribution SDs.
xFamily=[(1:16)'/16, ((-1).^(1:16)'.*(1:16)')/32];
family=acc01_exact_signflip_family(xFamily);
familySigns=1-2*(dec2bin(0:2^16-1)-'0');familySigns(1,:)=1;
familyStats=(familySigns*xFamily)/16;familySd=std(familyStats,0,1);
familyScale=familySd;familyScale(familyScale==0)=1;
familyStudentized=familyStats./familyScale;
familyMax=max(abs(familyStudentized),[],2);
familyPMax=mean(familyMax>=max(abs(familyStudentized(1,:))));
familyOk=max(abs(family.studentization_sd-familySd))<1e-14&& ...
    max(abs(family.studentized_statistics(:)-familyStudentized(:)))<1e-14&& ...
    family.p_maxstat==familyPMax&&isequal(family.observed,family.statistics(1,:));
add("family studentization independent recompute",familyOk,sprintf("p_maxstat=%.17g",family.p_maxstat));

% The zero-shift refit uses the exact same schedule helper as observed.
refitPair=struct('T',20,'signal',sin((0:5000)'/13),'bad_segments_s',zeros(0,2));
refitEvents=ct(1:60,:);
refitEvents.label_time_s=linspace(1,19,60)';
refitEvents.crossfit_fold=min(c.phase.erp_crossfit_folds,ceil((1:60)'*c.phase.erp_crossfit_folds/60));
[refitObserved,~]=acc01_phase_ppc_for_schedule(refitPair,refitEvents,c);
zeroEvents=refitEvents;zeroEvents.label_time_s=mod(double(refitEvents.label_time_s)+0,refitPair.T);
[refitZero,~]=acc01_phase_ppc_for_schedule(refitPair,zeroEvents,c);
add("refit(0) equals observed",abs(refitZero-refitObserved)<1e-12, ...
    sprintf("abs_diff=%.17g",abs(refitZero-refitObserved)));

% An impulse before a bad segment must not leak into the following clean
% block. A whole-trace filtfilt would fail this independent check.
blockPair=struct('T',20,'signal',zeros(5001,1),'bad_segments_s',[8 12]);
blockPair.signal(501)=1;
[blockFiltered,blockDiag]=acc01_filter_clean_blocks(blockPair.signal,blockPair,c);
postBlock=find(blockDiag.block_start_index>ceil(12*c.sampling.analysis_hz),1);
postIx=blockDiag.block_start_index(postBlock):blockDiag.block_stop_index(postBlock);
blockOk=blockDiag.n_blocks==2&&~isempty(postIx)&&max(abs(blockFiltered(postIx)))<1e-15;
add("per-block filter does not cross bad segment",blockOk, ...
    sprintf("blocks=%d post_block_max=%.17g",blockDiag.n_blocks,max(abs(blockFiltered(postIx)))));

[census,emptyKeys]=feasibleCensus(pairs,ct,c);
add("census nonempty",nnz(census.quantized_feasible_count>0)==351&&isempty(emptyKeys),sprintf("%d/351",nnz(census.quantized_feasible_count>0)));
add("census macro-no-bad pairs",nnz(census.n_bad_segments==0)==324,sprintf("%d/351",nnz(census.n_bad_segments==0)));
add("census minimum quanta",min(census.quantized_feasible_count)==62,sprintf("%d",min(census.quantized_feasible_count)));
add("census below 1pct",nnz(census.domain_fraction<.01)==351,sprintf("%d",nnz(census.domain_fraction<.01)));
add("census median fraction",abs(median(census.domain_fraction)-0.00958148859744382)<2e-5,sprintf("%.17g",median(census.domain_fraction)));
add("census median duration",abs(median(census.domain_duration_s)-4.588)<=.008,sprintf("%.17g",median(census.domain_duration_s)));
add("census all below 5pct",all(census.domain_fraction<.05),sprintf("max=%.9f",max(census.domain_fraction)));

expectedSession=["P41CS_R1" 14 0 0.00961466173264391 0.00961466173264391 0.00961466173264391;"P41CS_R2" 14 0 0.00960630842097349 0.00312413856473399 0.00960630842097349; ...
"P42CS_R1" 14 0 0.00958984211845293 0.00958984211845293 0.00958984211845293;"P42CS_R2" 14 0 0.00958984211845293 0.00958984211845293 0.00958984211845293;"P43CS_R1" 14 0 0.00958984211845293 0.00958984211845293 0.00958984211845293; ...
"P43CS_R2" 14 0 0.00962293466094191 0.00962293466094191 0.00962293466094191;"P44CS_R1" 14 0 0.00957321504648773 0.00957321504648773 0.00957321504648773;"P47CS_R1" 14 0 0.00791578806439736 0.00791578806439736 0.00791578806439736; ...
"P47CS_R2" 14 0 0.00790742042373732 0.00790742042373732 0.00790742042373732;"P48CS_R1" 14 0 0.00958148859744382 0.00958148859744382 0.00958148859744382;"P48CS_R2" 14 0 0.00953184521820491 0.00953184521820491 0.00953184521820491; ...
"P49CS_R1" 12 0 0.00958148859744382 0.00958148859744382 0.00958148859744382;"P49CS_R2" 14 0 0.00661598863921143 0.000517918302564531 0.00958148859744382;"P51CS_R1" 14 0 0.00958148859744382 0.00958148859744382 0.00958148859744382; ...
"P51CS_R2" 12 0 0.00787421341545053 0.00787421341545053 0.00787421341545053;"P53CS_R1" 12 0 0.00958984211845293 0.00261465207584997 0.00958984211845293;"P53CS_R2" 8 0 0.00959811546140288 0.00104418140354688 0.00959811546140288; ...
"P54CS_R1" 7 0 0.00103584525808419 0.00103584525808419 0.00365051917566766;"P54CS_R2" 7 0 0.00958984211845293 0.00958984211845293 0.00958984211845293;"P55CS_R1" 14 0 0.00958148859744382 0.00958148859744382 0.00958148859744382; ...
"P55CS_R2" 14 0 0.00955666752988547 0.00381765477373921 0.00955666752988547;"P56CS_R1" 12 0 0.00958984211845293 0.00958984211845293 0.00958984211845293;"P56CS_R2" 12 0 0.00956494135730277 0.00956494135730277 0.00956494135730277; ...
"P57CS_R1" 12 0 0.00958984211845293 0.00958984211845293 0.00958984211845293;"P57CS_R2" 12 0 0.00958984211845293 0.00104419012613817 0.00958984211845293;"P58CS_R1" 7 0 0.00958148859744382 0.00958148859744382 0.00958148859744382; ...
"P60CS_R1" 14 0 0.00957321504648773 0.00957321504648773 0.00957321504648773;"P62CS_R1" 7 0 0.00956494135730277 0.00956494135730277 0.00956494135730277;"P62CS_R2" 7 0 0.00956494135730277 0.00956494135730277 0.00956494135730277];
expectedSession=table(expectedSession(:,1),double(expectedSession(:,2)),double(expectedSession(:,3)),double(expectedSession(:,4)),double(expectedSession(:,5)),double(expectedSession(:,6)), ...
    'VariableNames',["session_id" "pairs" "infeasible" "median_fraction" "min_fraction" "max_fraction"]);
sessionRows=summarizeSessions(census);tol=1.1e-5;
for i=1:height(expectedSession)
    a=sessionRows(sessionRows.session_id==expectedSession.session_id(i),:);ok=height(a)==1&&a.pairs==expectedSession.pairs(i)&&a.infeasible==expectedSession.infeasible(i)&& ...
        max(abs([a.median_fraction a.min_fraction a.max_fraction]-[expectedSession.median_fraction(i) expectedSession.min_fraction(i) expectedSession.max_fraction(i)]))<=tol;
    add("census session "+expectedSession.session_id(i),ok,ternary(height(a)==1,sprintf("%.6f [%.6f,%.6f]",a.median_fraction,a.min_fraction,a.max_fraction),"missing"));
end

if p.Results.RunOptimizationEquivalence
    [optimizationOk,optimizationDetail]=optimizationEquivalence(root,c);
    add("optimized two-world equivalence",optimizationOk,optimizationDetail);
end

checks=table(names,status,detail,'VariableNames',["check" "status" "detail"]);summary=struct('pass',nnz(status=="PASS"),'warn',nnz(status=="WARN"),'fail',nnz(status=="FAIL"));
out=string(p.Results.OutputDir);
if strlength(out)==0,out=fullfile(root,'processed','subject','group','ieeg_acc00_sim','BangYoureDead');end
if ~isfolder(out),mkdir(out),end
writetable(checks,fullfile(out,'task-bangyouredead_desc-acc00sim-validation.tsv'),'FileType','text','Delimiter','\t');
writetable(census,fullfile(out,'task-bangyouredead_desc-acc00sim-feasible-domain-census.tsv'),'FileType','text','Delimiter','\t');
writetable(sessionRows,fullfile(out,'task-bangyouredead_desc-acc00sim-feasible-domain-sessions.tsv'),'FileType','text','Delimiter','\t');
report=struct('checks',checks,'summary',summary,'contract_sha256',contractHash,'census',census,'session_census',sessionRows);
fprintf("ACC00-SIM validator: %d pass / %d warn / %d fail\n",summary.pass,summary.warn,summary.fail);assert(summary.fail==0,"ACC00-SIM validator failed %d checks.",summary.fail);

function add(name,ok,why),names(end+1,1)=string(name);if ok,status(end+1,1)="PASS";else,status(end+1,1)="FAIL";end;detail(end+1,1)=string(why);end
end

function [ok,detail]=optimizationEquivalence(root,c)
% Compare the frozen runner with the new runner on the two archived null
% seeds. The gate is per-pair: all pair-row metric/diagnostic fields and the
% top-level robustfit diagnostics are compared. Timing and downstream family
% aggregates are excluded because the family implementation is unchanged and
% studentization can amplify sub-ulp pair-level differences.
addpath(fullfile(root,"scripts","ieeg","acc01"));
seedIndex=[1;2];maxDiff=zeros(numel(seedIndex),1);messages=strings(numel(seedIndex),1);
for ii=1:numel(seedIndex)
    seed=c.randomization.base_seed+100000+seedIndex(ii);
    fprintf("ACC00-SIM optimized equivalence: generating seed_index %d serially\n",seedIndex(ii));
    w=acc00_sim_generate_world(root,c,"benchmark",seed,"null",struct());
    reference=acc01_run_metrics(w,c);clear w
    w=acc00_sim_generate_world(root,c,"benchmark",seed,"null",struct());
    optimized=acc01_run_metrics_optimized(w,c);clear w
    referenceComparable=struct("pair_rows",reference.pair_rows, ...
        "diagnostics",reference.diagnostics);
    optimizedComparable=struct("pair_rows",optimized.pair_rows, ...
        "diagnostics",optimized.diagnostics);
    [maxDiff(ii),messages(ii)]=compareEquivalenceValue( ...
        referenceComparable,optimizedComparable,"seed_"+seedIndex(ii));
end
ok=all(maxDiff<1e-12);
detail=strjoin("seed_"+string(seedIndex)+" max_abs_diff="+compose("%.17g",maxDiff),"; ");
if any(messages~=""),detail=detail+"; "+strjoin(messages(messages~="")," | ");end
end

function [d,msg]=compareEquivalenceValue(a,b,label)
d=0;msg="";
if isstruct(a)&&isstruct(b)
    if ~isequal(size(a),size(b)),msg=label+" size mismatch";return,end
    fa=string(fieldnames(a));fb=string(fieldnames(b));
    if ~isequal(fa,fb),msg=label+" field mismatch";return,end
    for k=1:numel(a)
        for f=fa'
            [q,m]=compareEquivalenceValue(a(k).(f),b(k).(f),label+"."+f);
            d=max(d,q);if msg==""&&m~="",msg=m;end
        end
    end
elseif istable(a)&&istable(b)
    if height(a)~=height(b)||~isequal(a.Properties.VariableNames,b.Properties.VariableNames)
        msg=label+" table schema mismatch";return
    end
    for f=string(a.Properties.VariableNames)
        [q,m]=compareEquivalenceValue(a.(f),b.(f),label+"."+f);
        d=max(d,q);if msg==""&&m~="",msg=m;end
    end
elseif iscell(a)&&iscell(b)
    if ~isequal(size(a),size(b)),msg=label+" cell size mismatch";return,end
    for k=1:numel(a)
        [q,m]=compareEquivalenceValue(a{k},b{k},label+"{"+k+"}");
        d=max(d,q);if msg==""&&m~="",msg=m;end
    end
elseif isnumeric(a)&&isnumeric(b)
    if ~isequal(size(a),size(b)),msg=label+" numeric size mismatch";return,end
    aa=double(a);bb=double(b);bad=isnan(aa)~=isnan(bb)|isinf(aa)~=isinf(bb);
    if any(bad(:)),msg=label+" nonfinite mismatch";return,end
    finite=isfinite(aa)&isfinite(bb);
    if any(finite(:)),d=max(abs(aa(finite)-bb(finite)),[],'all');end
elseif islogical(a)&&islogical(b)
    if ~isequal(a,b),msg=label+" logical mismatch";end
elseif (isstring(a)||ischar(a))&&(isstring(b)||ischar(b))
    if ~isequal(string(a),string(b)),msg=label+" text mismatch";end
else
    if ~isequaln(a,b),msg=label+" value mismatch";end
end
end

function [out,emptyKeys]=feasibleCensus(pairs,covariates,c)
timing=readtable(c.inputs.i00_timing.path,'FileType','text','Delimiter','\t','TextType','string','VariableNamingRule','preserve');
bad=readtable(c.inputs.i00_bad_segments.path,'FileType','text','Delimiter','\t','TextType','string','VariableNamingRule','preserve');
n=height(pairs);session_id=pairs.session_id;pair_id=pairs.pair_name;subject=pairs.subject;run=pairs.run;domain_fraction=nan(n,1);domain_duration_s=nan(n,1);interval_count=zeros(n,1);quantized_feasible_count=zeros(n,1);n_bad_segments=zeros(n,1);
rs=RandStream("Threefry","Seed",1);
for i=1:n
    tr=timing(timing.session_id==pairs.session_id(i),:);T=double(tr.usable_movie_stop_s);zero=double(tr.macro_common_start_s)+(double(tr.macro_movie_first_sample)-1)/double(tr.macro_rate_hz);
    br=bad(bad.session_id==pairs.session_id(i)&bad.series=="LFP_macro"&ismember(bad.origchannel_name,[pairs.contact_a(i),pairs.contact_b(i)]),:);
    n_bad_segments(i)=height(br);iv=[double(br.start_common_time_s),double(br.stop_common_time_s)]-zero;iv(:,1)=max(iv(:,1),0);iv(:,2)=min(iv(:,2),T);iv=mergeOwn(iv(iv(:,2)>iv(:,1),:));
    pair=struct('T',T,'bad_segments_s',iv);[~,ev]=acc01_build_support_and_events(pair,covariates,c);s=acc01_build_empirical_surrogates(pair,ev,c,rs);
    domain_fraction(i)=s.domain_fraction;domain_duration_s(i)=s.domain_duration_s;interval_count(i)=s.interval_count;quantized_feasible_count(i)=s.quantized_feasible_count;
end
out=table(subject,session_id,run,pair_id,n_bad_segments,domain_fraction,domain_duration_s,interval_count,quantized_feasible_count);emptyKeys=session_id(quantized_feasible_count==0);
end
function out=summarizeSessions(c)
ids=unique(c.session_id,"sorted");pairs=zeros(numel(ids),1);infeasible=pairs;median_fraction=nan(numel(ids),1);min_fraction=median_fraction;max_fraction=median_fraction;
for i=1:numel(ids),z=c(c.session_id==ids(i),:);pairs(i)=height(z);infeasible(i)=nnz(z.quantized_feasible_count==0);median_fraction(i)=median(z.domain_fraction);min_fraction(i)=min(z.domain_fraction);max_fraction(i)=max(z.domain_fraction);end
out=table(ids,pairs,infeasible,median_fraction,min_fraction,max_fraction,'VariableNames',["session_id" "pairs" "infeasible" "median_fraction" "min_fraction" "max_fraction"]);
end
function out=mergeOwn(iv)
if isempty(iv),out=zeros(0,2);return,end;iv=sortrows(iv,[1 2]);out=iv(1,:);for i=2:size(iv,1),if iv(i,1)<=out(end,2)+1e-12,out(end,2)=max(out(end,2),iv(i,2));else,out(end+1,:)=iv(i,:);end,end %#ok<AGROW>
end
function s=ternary(cond,a,b),if cond,s=a;else,s=b;end,end
function hash=sha256File(path)
fid=fopen(path,"rb");assert(fid>=0,"Cannot open %s",path);cl=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance("SHA-256");while true,q=fread(fid,1048576,"*uint8");if isempty(q),break,end;md.update(typecast(q,"int8"));end;hash=upper(string(reshape(dec2hex(typecast(md.digest(),"uint8"),2)',1,[])));
end
