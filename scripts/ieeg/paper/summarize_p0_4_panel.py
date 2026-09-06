#!/usr/bin/env python3
"""Read-only, opaque-only, 480-checkpoint-gated P0-4 aggregation."""
from __future__ import annotations

import argparse, csv, hashlib, json, math, subprocess, sys, tempfile
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from scipy.special import betaincinv

TABLES=("p0_4_manifest.tsv","p0_4_member_worlds.tsv","p0_4_member_fpr_summary.tsv",
        "p0_4_paired_contrast_MD2.tsv","p0_4_exploratory_MD1.tsv","p0_4_blind_decode.tsv",
        "p0_4_aggregation_manifest.json")
class AggregationError(RuntimeError): pass
def sha256(p:Path)->str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest().upper()
def cp95(k:int,n:int)->tuple[float,float]:
    return (0. if k==0 else float(betaincinv(k,n-k+1,.025)),1. if k==n else float(betaincinv(k+1,n-k,.975)))
def write_tsv(p:Path,rows:list[dict[str,Any]])->None:
    if not rows:raise AggregationError(f"empty table refused: {p.name}")
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter="\t",extrasaction="raise");w.writeheader();w.writerows(rows)
def expected(root:Path,c:dict)->list[Path]:
    out=root/c["outputs"]["root"]/"checkpoints"
    return [out/oid/f"world_{i:03d}.mat" for oid in c["panel"]["opaque_members"] for i in range(1,121)]
def require_complete(root:Path,c:dict)->list[Path]:
    p=expected(root,c);found=list((root/c["outputs"]["root"]/"checkpoints").glob("**/world_*.mat"))
    if len(found)!=480:raise AggregationError(f"checkpoint cardinality {len(found)} != 480")
    if any(not x.is_file() for x in p):raise AggregationError("expected 480-world panel member checkpoint missing")
    if list((root/c["outputs"]["root"]/"checkpoints").glob("**/*.tmp.mat")):raise AggregationError("temporary checkpoint residue")
    return p
def decode(root:Path,paths:list[Path])->list[dict]:
    with tempfile.TemporaryDirectory(prefix="p04_meta_",dir=root) as td:
        td=Path(td);result=[]
        for start in range(0,len(paths),120):
            part=paths[start:start+120];lst=td/f"{start}.txt";out=td/f"{start}.json";lst.write_text("\n".join(str(p.resolve()) for p in part),encoding="utf-8")
            cmd=f"cd('{root.as_posix()}'); addpath('scripts/ieeg/paper'); export_acc00_p0_4_checkpoint_json('{lst.as_posix()}','{out.as_posix()}');"
            r=subprocess.run(["matlab","-batch",cmd],capture_output=True,text=True,encoding="utf-8",errors="replace")
            if r.returncode:raise AggregationError("MATLAB metadata export failed: "+r.stderr[-500:])
            q=json.loads(out.read_text(encoding="utf-8"));result.extend(q if isinstance(q,list) else [q])
        return result
def p_for(q:dict,run:str)->float:
    return float(q["signflip"][run]["p_two_sided"][0]) if q["arm"]=="Mpower" else float(q["mphase_shape"][run]["p_two_sided"])
def aggregate(root:Path,c:dict)->None:
    paths=require_complete(root,c);records=decode(root,paths)
    if len(records)!=480:raise AggregationError("decoded checkpoint cardinality mismatch")
    worlds=[]
    for r in records:
        q=r["checkpoint"];common={"opaque_id":q["opaque_id"],"seed_index":q["seed_index"],"world_seed":q["world_seed"],"checkpoint_sha256":sha256(Path(r["filename"]))}
        for run in ("R1","R2"):
            if q["arm"]=="Mpower": observed=q["signflip"][run]["observed"][0]
            else: observed=q["mphase_shape"][run]["observed"]
            worlds.append({**common,"run":run,"observed":observed,"p_two_sided":p_for(q,run),"reject_alpha_0_05":int(p_for(q,run)<=.05),"domain_fraction":q["domain"]["domain_fraction"],"elapsed_s":q["timing"]["elapsed_s"],"peak_memory_bytes":q["timing"]["peak_memory_bytes"]})
    groups=defaultdict(list)
    for x in worlds:groups[x["opaque_id"],x["run"]].append(x)
    summary=[]
    for (oid,run),rows in sorted(groups.items()):
        k=sum(x["reject_alpha_0_05"] for x in rows);n=len(rows);lo,hi=cp95(k,n);ph=k/n
        summary.append({"opaque_id":oid,"run":run,"n_worlds":n,"rejections_alpha_0_05":k,"fpr":ph,"cp95_low":lo,"cp95_high":hi,"mcse":math.sqrt(ph*(1-ph)/n),"observed_mean":sum(float(x["observed"]) for x in rows)/n})
    ids=c["panel"]["pre_registered_directional_contrast_opaque_ids"];by={(x["opaque_id"],x["seed_index"],x["run"]):x for x in worlds};contrast=[]
    for run in ("R1","R2"):
        k10=k01=0
        for i in range(1,121):
            a=by[ids[0],i,run]["reject_alpha_0_05"];b=by[ids[1],i,run]["reject_alpha_0_05"];k10+=int(a and not b);k01+=int(b and not a)
        n=k10+k01;p=1. if n==0 else min(1.,2*sum(math.comb(n,j) for j in range(0,min(k10,k01)+1))/2**n)
        contrast.append({"opaque_mutated":ids[0],"opaque_reference":ids[1],"run":run,"K10":k10,"K01":k01,"risk_difference":(k10-k01)/120,"exact_mcnemar_p_two_sided":p})
    exploratory=[x for x in worlds if x["opaque_id"]==c["panel"]["pre_registered_exploratory_opaque_id"]]
    out=root/c["outputs"]["tables"];out.mkdir(parents=True,exist_ok=True)
    manifest={"analysis_id":c["analysis_id"],"contract_sha256":sha256(root/"scripts/ieeg/paper/acc00_sim_p0_4_panel_contract.json"),"validated_checkpoints":480,"created_at_utc":datetime.now(UTC).isoformat(),"opaque_only_scoring":True}
    write_tsv(out/"p0_4_manifest.tsv",[manifest]);write_tsv(out/"p0_4_member_worlds.tsv",worlds);write_tsv(out/"p0_4_member_fpr_summary.tsv",summary);write_tsv(out/"p0_4_paired_contrast_MD2.tsv",contrast);write_tsv(out/"p0_4_exploratory_MD1.tsv",exploratory)
    # Unblinding is a separate post-score operation; all preceding tables use opaque IDs only.
    seal=json.loads((root/"scripts/ieeg/paper"/c["panel"]["member_semantics_seal"]).read_text(encoding="utf-8"));decode_rows=[]
    for oid,x in seal["opaque_to_member"].items():decode_rows.append({"opaque_id":oid,"member_id":x["member_id"],"prediction":x["prediction"],"mutation":x["mutation"]})
    write_tsv(out/"p0_4_blind_decode.tsv",decode_rows)
    (out/"p0_4_aggregation_manifest.json").write_text(json.dumps({**manifest,"script_sha256":sha256(Path(__file__)),"checkpoint_sha256":[{"filename":str(p),"sha256":sha256(p)} for p in paths],"blind_seal_sha256":sha256(root/"scripts/ieeg/paper"/c["panel"]["member_semantics_seal"])},indent=2)+"\n",encoding="utf-8")
def self_test(root:Path,c:dict)->int:
    try:require_complete(root,c)
    except AggregationError:print("PASS: incomplete collection refused");return 0
    print("FAIL: self-test expected an incomplete collection",file=sys.stderr);return 1
def main()->int:
    ap=argparse.ArgumentParser();ap.add_argument("--root",type=Path,required=True);ap.add_argument("--self-test",action="store_true");args=ap.parse_args();root=args.root.resolve();c=json.loads((root/"scripts/ieeg/paper/acc00_sim_p0_4_panel_contract.json").read_text(encoding="utf-8"))
    if args.self_test:return self_test(root,c)
    try:aggregate(root,c);print("P0-4 aggregation complete");return 0
    except (AggregationError,OSError,KeyError,ValueError,TypeError) as e:print("P0-4 aggregation refused: "+str(e),file=sys.stderr);return 2
if __name__=="__main__":raise SystemExit(main())
