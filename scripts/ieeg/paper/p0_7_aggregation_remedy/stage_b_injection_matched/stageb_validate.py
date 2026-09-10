# Release port of the analysis-workspace tool of the same name; see README in this folder.
"""Independent Stage-B validator; reads only MATLAB bridge TSVs."""
from __future__ import annotations
import argparse,json,math,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import beta,skew
def _repo_root():
    for _q in Path(__file__).resolve().parents:
        if (_q / "MANIFEST.md").is_file():
            return _q
    raise RuntimeError("release repository root not found")
ROOT = _repo_root();VERIFY=ROOT/"scripts"/"ieeg"/"paper";sys.path.insert(0,str(VERIFY))
from verify_p0_6_paired_contrast import newcombe_method10  # noqa: E402

def cp(k,n): return (0.0 if k==0 else float(beta.ppf(.025,k,n-k+1)),1.0 if k==n else float(beta.ppf(.975,k+1,n-k)))
def signflip(x):
    n=len(x)
    if not n:return math.nan,math.nan,math.nan
    rows=np.arange(2**n,dtype=np.uint32)[:,None];bits=(rows>>np.arange(n-1,-1,-1,dtype=np.uint32))&1;signs=1.0-2.0*bits;signs[0,:]=1.0;stats=signs@x/n
    return float(stats[0]),float(np.mean(np.abs(stats)>=abs(stats[0]))),float(np.std(stats,ddof=1))
def close(a,b,tol=1e-12): return (math.isnan(a) and math.isnan(b)) or abs(a-b)<=tol
def validate(output:Path,mode:str):
    pairs=pd.read_csv(output/"validator_pair_rows.tsv",sep="\t");subs=pd.read_csv(output/"validator_subject_rows.tsv",sep="\t");worlds=pd.read_csv(output/"validator_world_rows.tsv",sep="\t");errors=[]
    flip_stats={}  # per (arm, reducer, run): one-flip budget and straddle counts
    expected=2 if mode=="smoke" else 60;seeds=sorted(pairs.seed_index.unique().tolist())
    if seeds!=list(range(1,expected+1)):errors.append(f"seed set is not 1:{expected}")
    if pairs.world_seed.min()<22360828 or pairs.world_seed.max()>22360887:errors.append("world seed outside Stage-B namespace")
    if set(pairs.arm)!= {"P0","P4"} or not (pairs.mphase_surrogate_reducer=="mean").all() or not (pairs.refit_surrogate_count==32).all():errors.append("arm/reducer/refit identity failure")
    for seed in seeds:
        p0=pairs[(pairs.seed_index==seed)&(pairs.arm=="P0")].sort_values(["subject","run","pair_id"]);p4=pairs[(pairs.seed_index==seed)&(pairs.arm=="P4")].sort_values(["subject","run","pair_id"])
        keys=["subject","run","session_id","pair_id","base_fingerprint"]
        if len(p0)!=len(p4) or not p0[keys].reset_index(drop=True).equals(p4[keys].reset_index(drop=True)):errors.append(f"seed {seed}: base pairing mismatch")
    # Injection / stream provenance checks promised in the contract validator_checklist.
    contract=json.loads((Path(__file__).resolve().parent/"stageb_mphase_power_contract.json").read_text(encoding="utf-8"))
    eps_cfg=contract["injection"]["epsilon_resultant_check"]
    if pairs.epsilon_resultant.isna().any():errors.append("epsilon_resultant has NaN rows")
    # World-mean resultant per arm (tolerances calibrated at ~1900 events/pair,
    # world-mean SE ~0.0006) plus loose per-pair corruption caps; small smoke rows
    # (12 events) are n-noise dominated and exempt, but a formal run must consist
    # entirely of large-n pairs.
    big=pairs.epsilon_count>=1000
    if mode=="formal" and not big.all():errors.append("formal run has pairs with <1000 injected events")
    pb=pairs[big]
    if (pb.loc[pb.arm=="P0","epsilon_resultant"]>float(eps_cfg["P0_pair_max_abs"])).any():errors.append("P0 per-pair epsilon resultant above corruption cap")
    if ((pb.loc[pb.arm=="P4","epsilon_resultant"]-float(eps_cfg["P4_world_mean_target_R4"])).abs()>float(eps_cfg["P4_pair_max_abs_error"])).any():errors.append("P4 per-pair epsilon resultant beyond corruption cap")
    wm=pb.groupby(["seed_index","arm"]).epsilon_resultant.mean()
    p0m=wm.xs("P0",level="arm");p4m=wm.xs("P4",level="arm")
    if (p0m>float(eps_cfg["P0_world_mean_max"])).any():errors.append("P0 world-mean epsilon resultant above uniform tolerance")
    if ((p4m-float(eps_cfg["P4_world_mean_target_R4"])).abs()>float(eps_cfg["P4_world_mean_max_abs_error"])).any():errors.append("P4 world-mean epsilon resultant deviates from R(4) beyond tolerance")
    for seed in seeds:
        pw=pairs[pairs.seed_index==seed]
        # pair_id labels repeat across subjects, so identity is (subject, session, pair).
        dup=pw.groupby(["arm","subject","session_id","pair_id"]).injection_stream_seed.nunique()
        if (dup!=1).any() or pw.injection_stream_seed.duplicated().any():errors.append(f"seed {seed}: injection stream seed not unique per arm/pair")
        j=pw.pivot_table(index=["subject","session_id","pair_id"],columns="arm",values="injection_stream_seed",aggfunc="first")
        if (j["P0"]==j["P4"]).any():errors.append(f"seed {seed}: P0/P4 share an injection stream seed")
        ww=worlds[worlds.seed_index==seed]
        wh=ww.pivot_table(index="reducer",columns="arm",values="metric_stream_initial_state_hash",aggfunc="first")
        if (wh["P0"]!=wh["P4"]).any():errors.append(f"seed {seed}: metric stream initial states differ across arms")
        if ww.base_stream_state_hash.nunique()!=1:errors.append(f"seed {seed}: base stream state hash inconsistent within world")
    if (worlds.truth_a_g06_uV!=4.61799281479342).any() or (worlds.truth_grid_id!="G06").any():errors.append("injection truth amplitude/grid_id mismatch")
    if (pairs.n_injected_phase_intersection!=pairs.n_phase_events).any() or (pairs.n_injected_events<pairs.n_phase_events).any():errors.append("injected/phase event count identity violated")
    rows=[]
    for seed in seeds:
        for arm in ("P0","P4"):
            for reducer in ("median","mean"):
                for run in ("R1","R2"):
                    pr=pairs[(pairs.seed_index==seed)&(pairs.arm==arm)&(pairs.run==run)]
                    vals=[]
                    for subject,g in pr.groupby("subject",sort=True):
                        v=g.corrected_mphase.to_numpy(float);v=v[np.isfinite(v)]
                        if not len(v):continue
                        stat=float(np.median(v) if reducer=="median" else np.mean(v));ps=float(skew(v,bias=False)) if len(v)>=3 and not np.all(v==v[0]) else math.nan
                        sr=subs[(subs.seed_index==seed)&(subs.arm==arm)&(subs.run==run)&(subs.reducer==reducer)&(subs.subject==subject)]
                        if len(sr)!=1:errors.append(f"{seed}/{arm}/{run}/{reducer}/{subject}: subject row")
                        else:
                            r=sr.iloc[0]
                            if not close(stat,float(r.Mphase)) or int(r.n_pairs_finite)!=len(v) or not close(ps,float(r.pair_corrected_skewness)):errors.append(f"{seed}/{arm}/{run}/{reducer}/{subject}: reducer mismatch")
                        vals.append(stat)
                    obs,p,sd=signflip(np.asarray(vals,float));wr=worlds[(worlds.seed_index==seed)&(worlds.arm==arm)&(worlds.run==run)&(worlds.reducer==reducer)]
                    if len(wr)!=1:errors.append(f"{seed}/{arm}/{run}/{reducer}: world row")
                    else:
                        r=wr.iloc[0]
                        fs=flip_stats.setdefault((arm,reducer,run),{"budget":0.0,"s005":0,"s050":0,"s075":0})
                        pc=float(r.p_two_sided)
                        fs["budget"]+=abs(p-pc)
                        fs["s005"]+=int((p<=.05)!=(pc<=.05))
                        fs["s050"]+=int((p<.5)!=(pc<.5))
                        fs["s075"]+=int((p>.75)!=(pc>.75))
                        if not close(obs,float(r.observed)) or abs(p-pc)>1/(2**max(len(vals),1))+1e-12 or not close(sd,float(r.studentization_sd)):errors.append(f"{seed}/{arm}/{run}/{reducer}: signflip mismatch")
                    rows.append(dict(seed_index=seed,arm=arm,reducer=reducer,run=run,observed=obs,p_two_sided=p,n_subjects=len(vals)))
    rw=pd.DataFrame(rows);fpr=[];dist=[]
    for arm in ("P0","P4"):
        for reducer in ("median","mean"):
            for run in ("R1","R2"):
                g=rw[(rw.arm==arm)&(rw.reducer==reducer)&(rw.run==run)];n=len(g);k=int((g.p_two_sided<=.05).sum());lo,hi=cp(k,n) if n else (math.nan,math.nan);fpr.append(dict(arm=arm,reducer=reducer,run=run,n_worlds=n,rejections=k,fpr=k/n if n else math.nan,cp_low=lo,cp_high=hi))
                below=int((g.p_two_sided<.5).sum());above=int((g.p_two_sided>.75).sum());bl,bh=cp(below,n) if n else (math.nan,math.nan);al,ah=cp(above,n) if n else (math.nan,math.nan);dist.append(dict(arm=arm,reducer=reducer,run=run,n_worlds=n,mean_p=float(g.p_two_sided.mean()) if n else math.nan,p_lt_0p5_count=below,p_lt_0p5_fraction=below/n if n else math.nan,p_lt_0p5_cp_low=bl,p_lt_0p5_cp_high=bh,p_gt_0p75_count=above,p_gt_0p75_fraction=above/n if n else math.nan,p_gt_0p75_cp_low=al,p_gt_0p75_cp_high=ah))
    paired=[]
    for reducer in ("mean","median"):
        for run in ("R1","R2"):
            x=rw[(rw.arm=="P4")&(rw.reducer==reducer)&(rw.run==run)][["seed_index","p_two_sided"]].rename(columns={"p_two_sided":"p4"});y=rw[(rw.arm=="P0")&(rw.reducer==reducer)&(rw.run==run)][["seed_index","p_two_sided"]].rename(columns={"p_two_sided":"p0"});z=x.merge(y,on="seed_index",validate="one_to_one");n=len(z);p4=z.p4<=.05;p0=z.p0<=.05;both=int((p4&p0).sum());p4only=int((p4&~p0).sum());p0only=int((~p4&p0).sum());neither=int((~p4&~p0).sum());d,lo,hi=newcombe_method10(p4only,p0only,both,n) if n else (math.nan,)*3;paired.append(dict(reducer=reducer,run=run,n=n,both_reject=both,p4_only_reject=p4only,p0_only_reject=p0only,neither_reject=neither,fpr_p4=(both+p4only)/n if n else math.nan,fpr_p0=(both+p0only)/n if n else math.nan,delta_p4_minus_p0=d,ci_low=lo,ci_high=hi))
    # Summary agreement with per-group tolerances: threshold-count columns may move
    # only by actually observed cross-threshold worlds in that (arm/reducer/run)
    # group; interval columns must be internally consistent with the printed counts.
    for name,frame in (("main_fpr.tsv",pd.DataFrame(fpr)),("main_p_distribution.tsv",pd.DataFrame(dist)),("main_paired_contrast.tsv",pd.DataFrame(paired))):
        got=pd.read_csv(output/name,sep="\t")
        if list(got.columns)!=list(frame.columns) or got.shape!=frame.shape:
            errors.append(f"{name}: schema mismatch");continue
        keycols=["arm","reducer","run"] if "arm" in frame.columns else ["reducer","run"]
        for _,frow in frame.iterrows():
            sel=got
            for kc in keycols:sel=sel[sel[kc]==frow[kc]]
            if len(sel)!=1:errors.append(f"{name}: row key mismatch for {tuple(frow[kc] for kc in keycols)}");continue
            mrow=sel.iloc[0]
            if "arm" in keycols:
                fs=flip_stats[(frow["arm"],frow["reducer"],frow["run"])];n_g=max(int(frow["n_worlds"]),1)
                atol_map={"mean_p":fs["budget"]/n_g+1e-12,"rejections":float(fs["s005"]),
                          "p_lt_0p5_count":float(fs["s050"]),"p_gt_0p75_count":float(fs["s075"]),
                          "fpr":fs["s005"]/n_g+1e-12,"p_lt_0p5_fraction":fs["s050"]/n_g+1e-12,
                          "p_gt_0p75_fraction":fs["s075"]/n_g+1e-12}
                cp_src={"cp_low":("rejections",0),"cp_high":("rejections",1),
                        "p_lt_0p5_cp_low":("p_lt_0p5_count",0),"p_lt_0p5_cp_high":("p_lt_0p5_count",1),
                        "p_gt_0p75_cp_low":("p_gt_0p75_count",0),"p_gt_0p75_cp_high":("p_gt_0p75_count",1)}
            else:
                f4=flip_stats[("P4",frow["reducer"],frow["run"])];f0=flip_stats[("P0",frow["reducer"],frow["run"])]
                ssum=f4["s005"]+f0["s005"];n_g=max(int(frow["n"]),1)
                atol_map={"both_reject":float(ssum),"p4_only_reject":float(ssum),"p0_only_reject":float(ssum),
                          "neither_reject":float(ssum),"fpr_p4":f4["s005"]/n_g+1e-12,"fpr_p0":f0["s005"]/n_g+1e-12}
                cp_src={}
            for col in frame.columns:
                if col in keycols:continue
                if col in cp_src:
                    ccol,cidx=cp_src[col];okc=close(float(mrow[col]),float(cp(int(mrow[ccol]),int(mrow["n_worlds"]))[cidx]))
                elif col in ("delta_p4_minus_p0","ci_low","ci_high"):
                    d2,l2,h2=newcombe_method10(int(mrow["p4_only_reject"]),int(mrow["p0_only_reject"]),int(mrow["both_reject"]),int(mrow["n"]))
                    okc=close(float(mrow[col]),{"delta_p4_minus_p0":d2,"ci_low":l2,"ci_high":h2}[col])
                elif pd.api.types.is_numeric_dtype(frame[col]):
                    okc=close(float(mrow[col]),float(frow[col]),atol_map.get(col,1e-12))
                else:
                    okc=str(mrow[col])==str(frow[col])
                if not okc:errors.append(f"{name}: mismatch in {col} for {tuple(frow[kc] for kc in keycols)}")
    if mode=="formal":
        f=pd.DataFrame(fpr);q=pd.DataFrame(paired);rule_a=all(r.cp_low<=.05<=r.cp_high for r in f[(f.arm=="P0")&(f.reducer=="mean")].itertuples());rule_b=all(r.ci_low>0 for r in q[q.reducer=="mean"].itertuples());expected_status="INCOMPLETE" if seeds!=list(range(1,61)) else ("CONFIRMED" if rule_a and rule_b else "NOT_CONFIRMED")
        d=json.loads((output/"main_decision.json").read_text(encoding="utf-8"));
        if d.get("status")!=expected_status:errors.append("decision mismatch")
    else:
        d=json.loads((output/"main_decision.json").read_text(encoding="utf-8"));
        if d.get("status")!="INCOMPLETE_SMOKE" or d.get("scientific_rules_evaluated") is not False:errors.append("smoke evaluated scientific rules")
    report={"status":"FAIL if errors else PASS","mode":mode,"matched_checkpoints":len(seeds),"errors":errors,"checks":["matched arms/base fingerprints","epsilon distribution vs contract tolerances","injection/stream provenance","injection truth amplitude+grid_id","event count identities","pair reducers","exact sign-flip column 2","FPR/CP","Newcombe method 10","summary agreement","decision"]};report["status"]="FAIL" if errors else "PASS";(output/"validator_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8");return report
def main():
    p=argparse.ArgumentParser();p.add_argument("output",type=Path);p.add_argument("--mode",choices=("smoke","formal"),required=True);a=p.parse_args();r=validate(a.output,a.mode);print(json.dumps(r,indent=2));raise SystemExit(0 if r["status"]=="PASS" else 1)
if __name__=="__main__":main()
