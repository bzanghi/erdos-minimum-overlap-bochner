"""Run a single row's rigorous solve, save JSON. Args: row_idx (1-7)"""
import warnings; warnings.filterwarnings('ignore')
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from path_b_rigorous import solve_rigorous_at_center
from white_full_convex import WHITE_TABLE3

if __name__ == "__main__":
    row_idx = int(sys.argv[1])  # 1..7
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    T = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    R = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    bochner_n = int(sys.argv[5]) if len(sys.argv) > 5 else 20

    h, p, qm, qp, label = WHITE_TABLE3[row_idx - 1]
    print(f"=== {label}: h={h}, p={p}, q in [{qm},{qp}] | N={N} T={T} R={R} b={bochner_n} ===", flush=True)
    c = solve_rigorous_at_center(N, T, R, h, p, qm, qp, bochner_n)
    print(f"  reported : {c['reported_value']}", flush=True)
    print(f"  rigorous : {c['rigorous_dual_LB']}", flush=True)
    print(f"  delta    : {c['reported_value'] - c['rigorous_dual_LB']:.3e}", flush=True)
    print(f"  resid    : {c['dual_residual_at_LB']:.3e}", flush=True)
    print(f"  iters    : {c['best_iter']}/{c['n_iters_total']} ({c['n_eligible_iters']} eligible)", flush=True)
    print(f"  time     : {c['time']:.1f}s", flush=True)

    out_file = f"/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/path_b/{label}_rigorous.json"
    with open(out_file, 'w') as f:
        json.dump({
            'label': label,
            'h_c': h, 'p_c': p, 'q1': qm, 'q2': qp,
            'config': {'N': N, 'T': T, 'R': R, 'bochner_n': bochner_n},
            'reported': c['reported_value'],
            'rigorous_dual_LB': c['rigorous_dual_LB'],
            'dual_residual_at_LB': c['dual_residual_at_LB'],
            'best_iter': c['best_iter'],
            'n_iters_total': c['n_iters_total'],
            'n_eligible_iters': c['n_eligible_iters'],
            'status': c['status'],
            'time_s': c['time'],
            'duals': c['duals'],
        }, f, indent=2)
    print(f"  saved to {out_file}", flush=True)
