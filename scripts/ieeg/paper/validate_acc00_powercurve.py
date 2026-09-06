#!/usr/bin/env python3
"""Static validator for the frozen paper-only power-curve contract/profile."""
import json, math, hashlib, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
P=ROOT/"scripts/ieeg/paper/acc00_sim_powercurve_contract.json"
A=ROOT/"scripts/ieeg/acc00_sim_contract.json"
NULL=ROOT/"processed/subject/group/ieeg_acc00_sim/BangYoureDead/checkpoints/nullpool_edgeguard_20260829/nullpool_complete.json"
HIST=ROOT/"processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/contract/historical_null_contract_E0930B71.json"
PROF=ROOT/"processed/subject/group/ieeg_methods_paper/p0_2_design/profiling"
PROF_V2=ROOT/"processed/subject/group/ieeg_methods_paper/p0_2_sim_powercurve/profiling/v2"
NULL_CONTRACT_SHA="E0930B71CE91CD1641EDD9B1F7A1012C7DBB7769458ACB96386CA5892EAF8CEF"

def sha(path):
 h=hashlib.sha256(); h.update(path.read_bytes()); return h.hexdigest().upper()
def main():
 c=json.loads(P.read_text()); a=json.loads(A.read_text()); null=json.loads(NULL.read_text()); errors=[]; hist=json.loads(HIST.read_text()) if HIST.exists() else {}
 if null.get("contract_sha256") != NULL_CONTRACT_SHA:
  errors.append("null manifest contract SHA mismatch")
 for key,entry in c.get("inputs",{}).items():
  if null.get("input_hashes",{}).get(key) != entry.get("sha256"):
   errors.append(f"null/paper input hash mismatch: {key}")
 for k in ("sampling","background","windows_s","power","inputs"):
  if c.get(k)!=a.get(k): errors.append(f"section diff: {k}")
  if HIST.exists() and c.get(k)!=hist.get(k): errors.append(f"historical null diff: {k}")
 if c["sign_flip"]!=a["sign_flip"]: errors.append("section diff: sign_flip")
 if HIST.exists() and c["sign_flip"]!=hist.get("sign_flip"): errors.append("historical null diff: sign_flip")
 if c["scenario_index"]!=2 or sum(c["additive"]["grid_worlds"])!=986: errors.append("scenario/world total mismatch")
 C=c["label_gain_diagnostics"]["cov_z_exp2bz_over_var_z"]; E=math.exp(2*(.25**2+.15**2+.35**2));
 if abs(E-c["bridge"]["E_m_inverse2"])>1e-14: errors.append("E[m^-2] mismatch")
 G=C*E*c["bridge"]["K_R_log10_per_uV2"]
 if abs(G-c["bridge"]["expected_slope_gain_G_log10_per_uV2_per_z"])>1e-14: errors.append("bridge G mismatch")
 Agrid=c["additive"]["a_erp_grid_uV"][1:]; expected=[x*x*G for x in Agrid]
 if not all(abs(x-y)<1e-12 for x,y in zip([r*c["bridge"]["mde80_log10_per_z"] for r in c["additive"]["grid_r_mde"]],expected)): errors.append("grid expected slopes mismatch")
 seeds=[]
 for gi,n in enumerate(c["additive"]["grid_worlds"],1): seeds += [91000000+200000+1000*gi+i for i in range(1,n+1)]
 historical=set([20463828])
 for base,sc,eff,maxi in [(20260827,1,0,400),(20260827,2,4,100),(20260827,3,4,100)]:
  for e in range(0,eff+1):
   indices=range(1001,1065) if sc==3 else range(1,maxi+1)
   for i in indices: historical.add(base+100000*sc+1000*e+i)
 if historical.intersection(seeds): errors.append("seed collision")
 if c["theta_inversion"]["step"]<=0 or c["theta_inversion"]["hit_tolerance"]>c["theta_inversion"]["step"]/2: errors.append("theta grid invalid")
 nprof=0
 for gid in ("G01","G06","G07"):
  files=list(PROF.glob(f"world_g{int(gid[1:]):02d}_n001.mat"))
  if files: nprof+=1
 if nprof not in (0,3): errors.append("partial profiling set")
 if not PROF_V2.exists(): errors.append("missing v2 profiling directory")
 report={"status":"PASS" if not errors else "FAIL","errors":errors,"contract_sha256":sha(P),"acc00_contract_sha256":sha(A),"null_contract_sha256":null.get("contract_sha256"),"null_manifest_sha256":sha(NULL),"profiling_checkpoints":nprof,"world_total":sum(c["additive"]["grid_worlds"]),"seed_collision_count":len(historical.intersection(seeds))}
 out=ROOT/"processed/subject/group/ieeg_methods_paper/p0_2_design/validator_dry_run.json"; out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(report,indent=2)+"\n")
 print(json.dumps(report,indent=2)); return 0 if not errors else 1
if __name__=="__main__": sys.exit(main())
