# Release port of the analysis-workspace tool of the same name; see README in this folder.
"""Primary Stage-B summary; scientific rules are evaluated only in formal mode."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import pandas as pd
from scipy.stats import beta

def _repo_root():
    for _q in Path(__file__).resolve().parents:
        if (_q / "MANIFEST.md").is_file():
            return _q
    raise RuntimeError("release repository root not found")
ROOT = _repo_root()
VERIFY=ROOT/"scripts"/"ieeg"/"paper"
sys.path.insert(0,str(VERIFY))
from verify_p0_6_paired_contrast import newcombe_method10  # noqa: E402

def cp(k:int,n:int):
    return (0.0 if k==0 else float(beta.ppf(.025,k,n-k+1)),1.0 if k==n else float(beta.ppf(.975,k+1,n-k)))

def summarize(output:Path,mode:str)->dict:
    w=pd.read_csv(output/"validator_world_rows.tsv",sep="\t")
    fpr=[];dist=[]
    for arm in ("P0","P4"):
        for reducer in ("median","mean"):
            for run in ("R1","R2"):
                g=w[(w.arm==arm)&(w.reducer==reducer)&(w.run==run)].sort_values("seed_index");n=len(g);rej=g.p_two_sided<=.05;k=int(rej.sum());lo,hi=cp(k,n) if n else (float("nan"),)*2
                fpr.append(dict(arm=arm,reducer=reducer,run=run,n_worlds=n,rejections=k,fpr=k/n if n else float("nan"),cp_low=lo,cp_high=hi))
                below=int((g.p_two_sided<.5).sum());above=int((g.p_two_sided>.75).sum());bl,bh=cp(below,n) if n else (float("nan"),)*2;al,ah=cp(above,n) if n else (float("nan"),)*2
                dist.append(dict(arm=arm,reducer=reducer,run=run,n_worlds=n,mean_p=float(g.p_two_sided.mean()) if n else float("nan"),p_lt_0p5_count=below,p_lt_0p5_fraction=below/n if n else float("nan"),p_lt_0p5_cp_low=bl,p_lt_0p5_cp_high=bh,p_gt_0p75_count=above,p_gt_0p75_fraction=above/n if n else float("nan"),p_gt_0p75_cp_low=al,p_gt_0p75_cp_high=ah))
    paired=[]
    for reducer in ("mean","median"):
        for run in ("R1","R2"):
            x=w[(w.arm=="P4")&(w.reducer==reducer)&(w.run==run)][["seed_index","p_two_sided"]].rename(columns={"p_two_sided":"p4"})
            y=w[(w.arm=="P0")&(w.reducer==reducer)&(w.run==run)][["seed_index","p_two_sided"]].rename(columns={"p_two_sided":"p0"})
            z=x.merge(y,on="seed_index",validate="one_to_one");p4=z.p4<=.05;p0=z.p0<=.05;n=len(z);both=int((p4&p0).sum());p4only=int((p4&~p0).sum());p0only=int((~p4&p0).sum());neither=int((~p4&~p0).sum());d,lo,hi=newcombe_method10(p4only,p0only,both,n) if n else (float("nan"),)*3
            paired.append(dict(reducer=reducer,run=run,n=n,both_reject=both,p4_only_reject=p4only,p0_only_reject=p0only,neither_reject=neither,fpr_p4=(both+p4only)/n if n else float("nan"),fpr_p0=(both+p0only)/n if n else float("nan"),delta_p4_minus_p0=d,ci_low=lo,ci_high=hi))
    fprdf=pd.DataFrame(fpr);distdf=pd.DataFrame(dist);pairdf=pd.DataFrame(paired);fprdf.to_csv(output/"main_fpr.tsv",sep="\t",index=False);distdf.to_csv(output/"main_p_distribution.tsv",sep="\t",index=False);pairdf.to_csv(output/"main_paired_contrast.tsv",sep="\t",index=False)
    complete=sorted(w.seed_index.unique().tolist())==list(range(1,61)) and all(len(w[(w.arm==a)&(w.reducer==r)&(w.run==u)])==60 for a in ("P0","P4") for r in ("mean","median") for u in ("R1","R2"))
    if mode!="formal": decision={"analysis_id":"Bmovie-MPHASE-SUBJECT-REDUCER-POWER-STAGE-B-v1.0.0","status":"INCOMPLETE_SMOKE","scientific_rules_evaluated":False,"smoke_observed_difference_role":"implementation evidence only"}
    else:
        a=all(row.cp_low<=.05<=row.cp_high for row in fprdf[(fprdf.arm=="P0")&(fprdf.reducer=="mean")].itertuples());b=all(row.ci_low>0 for row in pairdf[pairdf.reducer=="mean"].itertuples());rules={"A_P0_mean_calibration":a,"B_mean_P4_minus_P0":b,"C_median_descriptive_only":True};decision={"analysis_id":"Bmovie-MPHASE-SUBJECT-REDUCER-POWER-STAGE-B-v1.0.0","status":"INCOMPLETE" if not complete else ("CONFIRMED" if a and b else "NOT_CONFIRMED"),"complete":complete,"rules":rules}
    (output/"main_decision.json").write_text(json.dumps(decision,indent=2),encoding="utf-8");return decision

def main():
    p=argparse.ArgumentParser();p.add_argument("output",type=Path);p.add_argument("--mode",choices=("smoke","formal"),required=True);a=p.parse_args();print(json.dumps(summarize(a.output,a.mode),indent=2))
if __name__=="__main__":main()
