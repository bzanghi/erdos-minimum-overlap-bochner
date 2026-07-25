import warnings, sys, json, time
warnings.filterwarnings('ignore')
sys.path.insert(0, __file__.rsplit('/',1)[0])
import cvxpy as cp
from _fside_ceiling import build_fexact
REPO='/Users/benzanghi/Documents/Claude/Projects/Erdos/.claude/worktrees/minimum-overlap-problem-5de0df'
db=json.load(open(f'{REPO}/lp_research_state/parallel_results/jansson_core12_reanchored.json'))
N=int(sys.argv[1]); M=int(sys.argv[2]); R=10
print(f"f-side CEILING vs production Jansson anchor   (N={N}, M={M}, R={R})")
print(f"{'center':16s} {'p_lo (prod, bn40+pm20)':>22s} {'V_fexact (ceiling)':>19s} {'max f-side headroom':>20s}")
rows=[]
for r in [x for x in db["centers"] if x.get("label","").startswith("cde")]:
    if not r.get('ok'): continue
    Om,cons,g = build_fexact(N,M,R,r['h_c'],r['h_c'],r['p_c'],r['p_c'],r['q1'],r['q2'])
    pr=cp.Problem(cp.Minimize(Om),cons); pr.solve(solver=cp.CLARABEL)
    hr = pr.value - r['p_lo']
    rows.append((r['label'], r['p_lo'], pr.value, hr, pr.status))
    print(f"{r['label']:16s} {r['p_lo']:22.9f} {pr.value:19.9f} {hr:+20.3e}  [{pr.status}]", flush=True)
print()
rows.sort(key=lambda t:t[1])
print("sorted by p_lo (the cover minimum is set by the lowest few):")
for lb,pl,vc,hr,st in rows[:5]:
    print(f"  {lb:16s} p_lo={pl:.9f}  ceiling={vc:.9f}  headroom={hr:+.3e}")
print(f"\ncurrent certified core floor = 0.3802946016 (binding: cde_n30_iter3)")
print(f"best possible new floor <= min_c(ceiling_c) - (transport eps) = {min(r[2] for r in rows):.9f} - eps")
