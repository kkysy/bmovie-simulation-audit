function out=run_acc00_p0_6_census_sensitivity(root)
%RUN_ACC00_P0_6_CENSUS_SENSITIVITY Feasible-domain census on 4 census-only
% schedules; per-pair enumeration follows validate_acc00_simulation.m
% feasibleCensus. Outputs never enter FPR/Newcombe/verdict aggregation.
arguments
    root (1,1) string = pwd
end
root=string(char(java.io.File(char(root)).getCanonicalPath()));
addpath(fullfile(root,"scripts","ieeg"));addpath(fullfile(root,"scripts","ieeg","acc01"));
c=jsondecode(fileread(fullfile(root,"scripts","ieeg","paper","acc00_sim_p0_6_perturbation_contract.json"))); %#ok<NASGU>
base=jsondecode(fileread(fullfile(root,"scripts","ieeg","acc00_sim_contract.json")));
sdir=fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation","schedules");
names=["schedule_census_intervalperm_1.tsv" "schedule_census_intervalperm_2.tsv" "schedule_census_eps_0.025xSD.tsv" "schedule_census_eps_0.10xSD.tsv"];
labels=["intervalperm_1" "intervalperm_2" "eps_0p025SD" "eps_0p10SD"];
% frozen 351-pair inventory from the I00 pair table (same as the main runs);
% contract paths are absolute project-root paths, do NOT prefix with root
pairs0=readtable(base.inputs.i00_acc_pairs.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
pairs0=pairs0(logical(pairs0.qc_eligible)&contains(pairs0.location,"ACC"),:);assert(height(pairs0)==351,"Expected 351 ACC pairs.");
timing=readtable(base.inputs.i00_timing.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
bad=readtable(base.inputs.i00_bad_segments.path,"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
od=fullfile(root,"processed","subject","group","ieeg_p0_6_perturbation");
detail=table(string([]),string([]),string([]),string([]),double([]),double([]),double([]),'VariableNames',["schedule" "subject" "session_id" "pair_name" "domain_fraction" "domain_duration_s" "quantized_feasible_count"]);
rs=RandStream("Threefry","Seed",1);
for s=1:numel(names)
    cov=readtable(fullfile(sdir,names(s)),"FileType","text","Delimiter","\t","TextType","string","VariableNamingRule","preserve");
    n=height(pairs0);df=nan(n,1);dd=nan(n,1);qc=zeros(n,1);
    for i=1:n
        tr=timing(timing.session_id==pairs0.session_id(i),:);T=double(tr.usable_movie_stop_s);
        movieZeroCommon=double(tr.macro_common_start_s)+(double(tr.macro_movie_first_sample)-1)/double(tr.macro_rate_hz);
        br=bad(bad.session_id==pairs0.session_id(i)&bad.series=="LFP_macro"&ismember(bad.origchannel_name,[pairs0.contact_a(i),pairs0.contact_b(i)]),:);
        iv=[double(br.start_common_time_s),double(br.stop_common_time_s)]-movieZeroCommon;iv(:,1)=max(iv(:,1),0);iv(:,2)=min(iv(:,2),T);
        pair=struct('T',T,'bad_segments_s',mergeIv(iv(iv(:,2)>iv(:,1),:)));
        [~,ev]=acc01_build_support_and_events(pair,cov,base);
        sur=acc01_build_empirical_surrogates(pair,ev,base,rs);
        df(i)=sur.domain_fraction;dd(i)=sur.domain_duration_s;qc(i)=sur.quantized_feasible_count;
    end
    detail=[detail;table(repmat(labels(s),n,1),string(pairs0.subject),string(pairs0.session_id),string(pairs0.pair_name),df,dd,qc,'VariableNames',["schedule" "subject" "session_id" "pair_name" "domain_fraction" "domain_duration_s" "quantized_feasible_count"])]; %#ok<AGROW>
    fprintf("census %s: median fraction %.6g, min quanta %d\n",labels(s),median(df,"omitnan"),min(qc));
end
writetable(detail,fullfile(od,"census_sensitivity.tsv"),"FileType","text","Delimiter","\t");
summary=varfun(@median,detail,'InputVariables',"domain_fraction",'GroupingVariable',"schedule");
out=summary;
    function iv=mergeIv(iv)
    if isempty(iv),iv=zeros(0,2);return,end;iv=sortrows(iv,[1 2]);o=iv(1,:);
    for k=2:size(iv,1),if iv(k,1)<=o(end,2)+1e-12,o(end,2)=max(o(end,2),iv(k,2));else,o(end+1,:)=iv(k,:);end,end %#ok<AGROW>
    iv=o;
    end
end
