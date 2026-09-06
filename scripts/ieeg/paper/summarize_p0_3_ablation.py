#!/usr/bin/env python3
"""Read-only, 768-checkpoint-gated P0-3 aggregation."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, subprocess, sys, tempfile
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from scipy.special import betaincinv

TABLES = (
    "p0_3_manifest.tsv", "p0_3_mpower_worlds.tsv", "p0_3_mpower_pair_diagnostics.tsv",
    "p0_3_mpower_cell_summary.tsv", "p0_3_mpower_factor_contrasts.tsv",
    "p0_3_mphase_worlds.tsv", "p0_3_mphase_replay_shape_summary.tsv",
    "p0_3_mphase_historical_guard_summary.tsv", "p0_3_aggregation_manifest.json",
)

class AggregationError(RuntimeError): pass

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest().upper()

def cp95(k:int,n:int) -> tuple[float,float]:
    return (0. if k==0 else float(betaincinv(k,n-k+1,.025)),1. if k==n else float(betaincinv(k+1,n-k,.975)))

def write_tsv(path:Path, rows:list[dict[str,Any]]) -> None:
    if not rows: raise AggregationError(f"empty output refused: {path.name}")
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter="\t",extrasaction="raise");w.writeheader();w.writerows(rows)

def expected(root:Path,c:dict[str,Any]) -> list[Path]:
    out=[];base=root/c["outputs"]["root"]
    for arm,key in (("mpower","Mpower"),("mphase","Mphase")):
        for cell in c["cells"][key]:
            out += [base/arm/"checkpoints"/cell["id"]/f"world_{i:03d}.mat" for i in range(1,int(cell["worlds"])+1)]
    return out

def require_complete(root:Path,c:dict[str,Any],base:Path|None=None) -> list[Path]:
    if base is None:
        paths=expected(root,c);output=root/c["outputs"]["root"]
        found=list((output/"mpower"/"checkpoints").glob("**/world_*.mat"))+list((output/"mphase"/"checkpoints").glob("**/world_*.mat"))
    else:
        paths=[]
        for arm,key in (("mpower","Mpower"),("mphase","Mphase")):
            for cell in c["cells"][key]: paths += [base/arm/"checkpoints"/cell["id"]/f"world_{i:03d}.mat" for i in range(1,int(cell["worlds"])+1)]
        found=list(base.glob("**/world_*.mat"))
    if len(found)!=768: raise AggregationError(f"checkpoint cardinality {len(found)} != 768")
    missing=[p for p in paths if not p.is_file()]
    if missing: raise AggregationError(f"missing expected checkpoint: {missing[0]}")
    tmp_roots=[base] if base is not None else [root/c["outputs"]["root"]/"mpower"/"checkpoints",root/c["outputs"]["root"]/"mphase"/"checkpoints"]
    if any(list(x.glob("**/*.tmp.mat")) for x in tmp_roots): raise AggregationError("temporary checkpoint residue")
    return paths

def decode(root:Path, paths:list[Path]) -> list[dict[str,Any]]:
    records=[]
    with tempfile.TemporaryDirectory(prefix="p03_meta_",dir=root) as td:
        td=Path(td)
        for start in range(0,len(paths),128):
            part=paths[start:start+128];lst=td/f"{start}.txt";out=td/f"{start}.json";lst.write_text("\n".join(str(x.resolve()) for x in part),encoding="utf-8")
            cmd=f"cd('{root.as_posix()}'); addpath('scripts/ieeg/paper'); export_acc00_p0_3_checkpoint_json('{lst.as_posix()}','{out.as_posix()}');"
            r=subprocess.run(["matlab","-batch",cmd],text=True,capture_output=True,encoding="utf-8",errors="replace")
            if r.returncode: raise AggregationError("MATLAB metadata export failed: "+r.stderr[-500:])
            x=json.loads(out.read_text(encoding="utf-8"));records.extend(x if isinstance(x,list) else [x])
    return records

def num(x:Any)->float: return math.nan if x is None else float(x)

def write_markdown(path:Path,title:str,intro:str,rows:list[dict[str,Any]])->None:
    if not rows: raise AggregationError(f"empty Markdown refused: {path.name}")
    fields=list(rows[0]);lines=[f"# {title}","",intro,"", "|"+"|".join(fields)+"|", "|"+"|".join("---" for _ in fields)+"|"]
    for row in rows: lines.append("|"+"|".join(str(row[k]).replace("|","\\|") for k in fields)+"|")
    path.write_text("\n".join(lines)+"\n",encoding="utf-8")

def aggregate(root:Path,c:dict[str,Any]) -> None:
    paths=require_complete(root,c);rec=decode(root,paths)
    if len(rec)!=768: raise AggregationError("decoded checkpoint cardinality mismatch")
    mpw=[];pairs=[];phw=[];shape=[]
    for r in rec:
        q=r["checkpoint"];fn=Path(r["filename"]);common={"cell_id":q["cell_id"],"seed_index":q["seed_index"],"world_seed":q["world_seed"],"checkpoint_sha256":sha256(fn)}
        if q["arm"]=="Mpower":
            for sm in q["subject_metrics"]:
                run=sm["run"];inf=q["signflip"][run];mpw.append({**common,"run":run,"observed_slope":inf["observed"][0],"p_two_sided":inf["p_two_sided"][0],"reject_alpha_0_05":int(num(inf["p_two_sided"][0])<=.05),"elapsed_s":q["timing"]["elapsed_s"],"peak_memory_bytes":q["timing"]["peak_memory_bytes"],"mask_sha256":q["injection_provenance"]["mask_sha256"]})
            for p in q["pair_diagnostics"]: pairs.append({**common,**p})
        else:
            for run,s in q["mphase_shape"].items():
                phw.append({**common,"run":run,"observed":s["observed"],"p_two_sided":s["p_two_sided"],"p_maxstat":s["p_maxstat"]})
                shape.append({**common,"run":run,**{k:v for k,v in s.items() if k not in ("subjects_lexicographic","subject_mphase_vector")}})
    grouped=defaultdict(list)
    for x in mpw: grouped[(x["cell_id"],x["run"])].append(x)
    summary=[]
    for (cell,run),x in sorted(grouped.items()):
        vals=[num(z["observed_slope"]) for z in x];k=sum(z["reject_alpha_0_05"] for z in x);lo,hi=cp95(k,len(x));summary.append({"cell_id":cell,"run":run,"n_worlds":len(x),"slope_mean":sum(vals)/len(vals),"slope_sample_sd":math.sqrt(sum((v-sum(vals)/len(vals))**2 for v in vals)/(len(vals)-1)),"rejections_alpha_0_05":k,"rejection_rate":k/len(x),"cp95_low":lo,"cp95_high":hi})
    by={(x["cell_id"],x["seed_index"],x["run"]):x for x in mpw};spec=[("pre_minus_post_canonical","P04","P01"),("pre_minus_post_half","P05","P02"),("pre_minus_post_quarter","P06","P03"),("symmetric_minus_post_canonical","P07","P01"),("symmetric_minus_post_half","P08","P02"),("symmetric_minus_post_quarter","P09","P03"),("post_half_minus_canonical","P02","P01"),("post_quarter_minus_canonical","P03","P01"),("P10_minus_P01","P10","P01")];contr=[]
    for label,a,b in spec:
        for run in ("R1","R2"):
            vals=[num(by[a,i,run]["observed_slope"])-num(by[b,i,run]["observed_slope"]) for i in range(1,65)];contr.append({"contrast":label,"run":run,"n_pairs":64,"paired_difference_mean":sum(vals)/64,"paired_difference_sample_sd":math.sqrt(sum((v-sum(vals)/64)**2 for v in vals)/63)})
    hist=[]
    for key in ("Mphase_guard_0p75","Mphase_guard_1p5"):
        ref=c["historical_references"][key];f=root/ref["path"]/"nullpool_fpr.tsv"
        for row in csv.DictReader(f.open(encoding="utf-8"),delimiter="\t"):
            if row.get("test")=="Mphase":hist.append({"reference_cell":key,"run":row["run"],"n_worlds":row["n_worlds"],"fpr":row["empirical_fpr"],"historical_contract_sha256":ref["contract_sha256"],"role":"guard_FPR_reference_n400"})
    out=root/c["outputs"]["tables"];out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="p03_tables_",dir=out.parent) as td:
        td=Path(td);write_tsv(td/"p0_3_manifest.tsv",[{"analysis_id":c["analysis_id"],"contract_sha256":sha256(root/"scripts/ieeg/paper/acc00_sim_p0_3_ablation_contract.json"),"validated_checkpoints":768,"created_at_utc":datetime.now(UTC).isoformat()}]);write_tsv(td/"p0_3_mpower_worlds.tsv",mpw);write_tsv(td/"p0_3_mpower_pair_diagnostics.tsv",pairs);write_tsv(td/"p0_3_mpower_cell_summary.tsv",summary);write_tsv(td/"p0_3_mpower_factor_contrasts.tsv",contr);write_tsv(td/"p0_3_mphase_worlds.tsv",phw);write_tsv(td/"p0_3_mphase_replay_shape_summary.tsv",shape);write_tsv(td/"p0_3_mphase_historical_guard_summary.tsv",hist);write_markdown(td/"p0_3_mpower_cell_summary.md","P0-3 Mpower ablation summary","Read-only complete-world summary; all fixed cells are retained.",summary);write_markdown(td/"p0_3_mphase_replay_shape_summary.md","P0-3 Mphase replay shape summary","The n=64 replay rows are descriptive only, not the guard-FPR reference.",shape);write_markdown(td/"p0_3_mphase_historical_guard_summary.md","P0-3 historical Mphase guard reference","Historical n=400 guard-FPR reference; it is not pooled with the n=64 replays.",hist);(td/"p0_3_aggregation_manifest.json").write_text(json.dumps({"status":"complete","script_sha256":sha256(Path(__file__)),"contract_sha256":sha256(root/"scripts/ieeg/paper/acc00_sim_p0_3_ablation_contract.json"),"checkpoints":[{"filename":str(p),"sha256":sha256(p)} for p in paths]},indent=2)+"\n",encoding="utf-8")
        for p in td.iterdir():os.replace(p,out/p.name)

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--smoke-self-test", action="store_true")
    args = ap.parse_args()
    root = args.root.resolve()
    contract = root / "scripts/ieeg/paper/acc00_sim_p0_3_ablation_contract.json"
    spec = json.loads(contract.read_text(encoding="utf-8"))
    if args.validate_only:
        print(json.dumps({"status": "read_only_implementation_ready", "expected_tables": TABLES}, indent=2)); return 0
    if args.smoke_self_test:
        try: require_complete(root,spec,root/spec["outputs"]["root"]/'smoke_v2')
        except AggregationError: print("PASS smoke self-test: incomplete smoke collection refused"); return 0
        raise AssertionError("incomplete smoke collection accepted")
    try: aggregate(root,spec);print("P0-3 aggregation complete");return 0
    except (AggregationError,OSError,KeyError,ValueError) as e: print("P0-3 aggregation refused: "+str(e),file=sys.stderr);return 2

if __name__ == "__main__":
    raise SystemExit(main())
