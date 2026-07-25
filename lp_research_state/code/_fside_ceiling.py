"""CEILING TEST.  White's program with the ENTIRE f-side made exact.

f enters White's Section-5 program ONLY through the 4R numbers
(a_m, b_m)_{m=1..2R}, a_m = (1/2)int_{-1}^{1} f cos(m pi x/2), b_m likewise sin.
Everything else on the f-side -- the c,d truncation at T, |c|<=2/pi,
sum(c^2+d^2)<=1/2, T5', the odd-m transfer relation and its eps/dlt tail
slacks, and EVERY Bochner / poly-moment / Lasserre cut, integer OR half-integer --
exists only to outer-approximate the achievable set

    K_R := { (a_m,b_m)_{m<=2R} : 0<=f<=1, int f = 1 }.

Here we replace that outer approximation by an INNER one: f piecewise constant
on M cells, a_m,b_m computed by exact quadrature.  The resulting optimum V_pwc(M)
satisfies   V_pwc(M) >= V_exact >= (value of White's program + ANY valid f-side
cut family, at the same N and R),  and V_pwc(M) decreases to V_exact as M grows.
So V_pwc is a CONSERVATIVE ceiling for the whole f-side lever class.
"""
import warnings, time, sys
warnings.filterwarnings('ignore')
sys.path.insert(0,'/Users/benzanghi/Documents/Claude/Projects/Erdos/.claude/worktrees/minimum-overlap-problem-5de0df/lp_research_state/code')
import numpy as np, cvxpy as cp
from white_full_convex import cos_cell_bounds_exact, sin_cell_bounds_exact

CENTERS={'row1':(0.015,0.381,-0.02,0.02),'row4':(0.004,0.3875,-0.02,0.02),
         'row7':(0.030,0.375,-0.02,0.02),'row6':(0.000,0.381,-0.02,0.02)}

def f_moment_matrices(M, R):
    """Icos[m-1], Isin[m-1]: length-M vectors with a_m = .5*Icos@g, b_m = .5*Isin@g."""
    xe = np.linspace(-1.0, 1.0, M+1)
    Ic, Is = [], []
    for m in range(1, 2*R+1):
        k = m*np.pi/2
        Ic.append((np.sin(k*xe[1:]) - np.sin(k*xe[:-1]))/k)
        Is.append((np.cos(k*xe[:-1]) - np.cos(k*xe[1:]))/k)
    return np.array(Ic), np.array(Is)

def build_fexact(N, M, R, h1, h2, p1, p2, q1, q2):
    L = 2.0/N; j = np.arange(1, N+1)
    Om = cp.Variable(); w = cp.Variable(N); v = cp.Variable(N); g = cp.Variable(M)
    cons = [w>=0, v>=0, w<=Om, v<=Om, Om<=1, g>=0, g<=1]
    cons.append(L*cp.sum(w+v) == 1)
    cons.append((2.0/M)*cp.sum(g) == 1)
    cons.append(L**2*cp.sum(cp.multiply(j,w)-cp.multiply(j-1,v)) >= h1)
    cons.append(L**3*cp.sum(cp.multiply((j-1)**2,(w+v))) <= 2.0/3 + h2**2/2)
    Ic, Is = f_moment_matrices(M, R)
    A = [0.5*(Ic[m-1] @ g) for m in range(1, 2*R+1)]
    B = [0.5*(Is[m-1] @ g) for m in range(1, 2*R+1)]
    for m in range(1, 2*R+1):
        am, bm = A[m-1], B[m-1]
        s = np.sin(np.pi*m/2)
        a_minus, _ = cos_cell_bounds_exact(j, m, L)
        cons.append((L/2)*(a_minus@(w+v)) + 2*cp.square(am) + 2*cp.square(bm)
                    - (4*s/(m*np.pi))*am <= 0)
        b_minus, b_plus = sin_cell_bounds_exact(j, m, L)
        rhs = -(4.0/(m*np.pi))*s*bm
        cons.append((L/2)*(b_minus@w - b_plus@v) <= rhs)
        cons.append((L/2)*(b_plus@w - b_minus@v) >= rhs)
    cons += [2*A[1] >= p1, 2*A[1] <= p2, 2*B[1] >= q1, 2*B[1] <= q2]   # c1=2a_2, d1=2b_2
    _, a_plus_2 = cos_cell_bounds_exact(j, 2, L)
    cons.append((L/2)*(a_plus_2@(w+v)) >= -0.5*(max(p1**2,p2**2)+max(q1**2,q2**2)))
    return Om, cons, g

if __name__ == '__main__':
    N = int(sys.argv[1]); R = 10
    for label in sys.argv[3].split(','):
        h,p,q1,q2 = CENTERS[label]
        for M in [int(x) for x in sys.argv[2].split(',')]:
            Om, cons, g = build_fexact(N, M, R, h, h, p, p, q1, q2)
            pr = cp.Problem(cp.Minimize(Om), cons)
            t0=time.time(); pr.solve(solver=cp.CLARABEL)
            print(f"{label}  N={N} M={M:5d}  V_pwc = {pr.value:.7f}  [{pr.status}] {time.time()-t0:.0f}s",
                  flush=True)
