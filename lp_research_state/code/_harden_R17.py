import sys; sys.path.insert(0,'.')
import json, time
from path_b_with_polymoment import solve_with_pm
N,T,R,BN,PMK = 24000,4000,10,40,20; TARGET=0.3802838
OUT="../parallel_results/harden_loadbearing_N24K.json"
d=json.load(open(OUT)); 
t0=time.time()
res=solve_with_pm(N,T,R,0.0,0.39,0.02,0.02,BN,PMK); val=float(res["value"])
d["results"]["R17_h0.0_p0.39_q0.02"]={"h":0.0,"p":0.39,"q1":0.02,"q2":0.02,"claimed_N20K":0.38063,
  "N24K_value":val,"status":res["status"],"delta_vs_N20K":val-0.38063,
  "margin_vs_target":val-TARGET,"sec":round(time.time()-t0)}
json.dump(d,open(OUT,"w"),indent=2)
print(f"R17_h0.0_p0.39_q0.02: N24K={val:.7f} margin_vs_target={val-TARGET:+.2e} status={res['status']} {time.time()-t0:.0f}s",flush=True)
print("DONE",flush=True)
