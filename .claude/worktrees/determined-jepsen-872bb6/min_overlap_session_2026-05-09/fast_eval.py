"""
Fast vectorized evaluation of all J^+_j and J^-_j (j = 1,...,N-1) for an
array h of length N. We compute the restricted autocorrelation via FFT.

Define
   R(j) = sum_{i=1}^{N-j} h_i h_{i+j}      (j = 0, 1, ..., N-1)
   S^+(j) = sum_{i=1}^{N-j} h_{i+j}        (j ≥ 0)
   S^-(j) = sum_{i=1}^{N-j} h_i            (j ≥ 0)

Then
   J^+_j = (2/N) (S^+(j) - R(j))   = (2/N) sum h_{i+j}(1 - h_i)
   J^-_j = (2/N) (S^-(j) - R(j))   = (2/N) sum h_i (1 - h_{i+j})

R(j) is the restricted autocorrelation; we obtain it via numpy.correlate.
S^+(j) and S^-(j) are partial sums, easy from the cumulative sum.

Gradients:
   ∂J^+_j/∂h_k = (2/N)·[ 1_{k>j} - h_{k-j}·1_{k>j} - h_{k+j}·1_{k+j ≤ N} ]·...
   (separately tracked for + and - directions; vectorized below).
"""
from __future__ import annotations
import numpy as np


def J_all(h: np.ndarray):
    """Returns (Jp, Jm) where Jp[j-1] = J^+_j, Jm[j-1] = J^-_j for j = 1..N-1."""
    N = len(h)
    Δ = 2.0 / N
    R = np.correlate(h, h, mode="full")
    R = R[N - 1 :]              # R[j] = sum_{i} h_i h_{i+j}, j = 0..N-1
    cum = np.cumsum(h)          # cum[i] = h_0 + ... + h_i  (0-indexed)
    total = cum[-1]
    # S^+(j) = sum_{i=0}^{N-1-j} h_{i+j} = sum_{k=j}^{N-1} h_k = total - cum[j-1]
    # S^-(j) = sum_{i=0}^{N-1-j} h_i = cum[N-1-j]
    j_arr = np.arange(1, N)
    Sp = total - cum[j_arr - 1]
    Sm = cum[N - 1 - j_arr]
    Rj = R[j_arr]
    Jp = Δ * (Sp - Rj)
    Jm = Δ * (Sm - Rj)
    return Jp, Jm


def smooth_max_grad_fast(h: np.ndarray, β: float):
    """Vectorized smoothed max + gradient over all 2(N-1) constraints."""
    N = len(h)
    Δ = 2.0 / N
    Jp, Jm = J_all(h)
    Js = np.concatenate([Jp, Jm])
    m = float(Js.max())
    w = np.exp(β * (Js - m))
    Z = float(w.sum())
    sm = m + np.log(Z) / β
    wp = w[: N - 1] / Z
    wm = w[N - 1 :] / Z
    # ∂J^+_j/∂h_k:
    #   J^+_j = Δ (S^+_j − R_j), where S^+_j = sum_{i=j+1}^{N} h_i,  R_j = sum_{i=1}^{N-j} h_i h_{i+j}.
    #   ∂S^+_j/∂h_k = 1_{k > j}        (1-indexed; k ∈ [j+1, N])
    #   ∂R_j/∂h_k   = h_{k-j} if k > j else 0,  plus  h_{k+j} if k+j ≤ N else 0.
    #   So ∂J^+_j/∂h_k = Δ [ 1_{k>j} − (h_{k-j} 1_{k>j} + h_{k+j} 1_{k+j≤N}) ].
    # Similarly ∂J^-_j/∂h_k = Δ [ 1_{k+j ≤ N} − (h_{k-j} 1_{k>j} + h_{k+j} 1_{k+j≤N}) ].
    grad = np.zeros(N)
    # Compute A_k := sum_j wp_j * 1_{k>j}   ; this equals sum_{j=1}^{k-1} wp_j (cumulative sum of wp through j=k-1)
    cumwp = np.concatenate([[0.0], np.cumsum(wp)])  # cumwp[k] = sum_{j=1}^{k} wp_{j-1}? careful
    # Let wp index match j-1 (j=1..N-1).  For k in 1..N (1-indexed), 1_{k>j} = 1_{j<k} = 1_{j ≤ k-1}.
    #    A_k = sum_{j=1}^{k-1} wp[j-1]  = cumsum_wp[k-1]
    # cumsum_wp[m] := sum_{i=0}^{m-1} wp[i], so cumsum_wp[0]=0, cumsum_wp[N-1] = total wp.
    cumsum_wp = np.concatenate([[0.0], np.cumsum(wp)])  # length N
    A = cumsum_wp[: N]   # A[k-1] = sum_{j=1}^{k-1} wp_j ; for k=1 this is 0
    # B_k := sum_j wm_j * 1_{k+j ≤ N} = sum_{j=1}^{N-k} wm_j  = cumsum_wm[N-k]
    cumsum_wm = np.concatenate([[0.0], np.cumsum(wm)])
    B = cumsum_wm[N - np.arange(1, N + 1)]  # B[k-1] for k=1..N
    # Cross terms:
    #   ∑_j (wp_j + wm_j) [ h_{k-j} 1_{k>j} + h_{k+j} 1_{k+j ≤ N} ]
    # Let w_total = wp + wm. We need:
    #   T1_k = ∑_{j=1}^{k-1} w_total[j-1] * h_{k-j}        (1-indexed: h_{k-j})
    #   T2_k = ∑_{j=1}^{N-k} w_total[j-1] * h_{k+j}
    w_total = wp + wm   # length N-1
    # T1_k = convolution evaluated at k (1..N): T1_k = (w_total ⋆ h)_k  computed carefully
    # T1_k = sum_{j=1}^{k-1} w_total[j-1] * h_{k-j}.  Let i = k-j (i = 1..k-1), so j = k-i.
    #   T1_k = sum_{i=1}^{k-1} h_i * w_total[k-i-1]
    # This is a "lower-triangular" convolution; np.convolve handles it.
    h_pad = h            # length N
    wt_pad = w_total     # length N-1
    full = np.convolve(h_pad, wt_pad, mode="full")  # length 2N-2
    # full[m] = sum_{i+j' = m} h_i * wt[j']  for i in [0, N-1], j' in [0, N-2]. (0-indexed)
    # We want T1_k = sum_{i=1}^{k-1} h_i * w_total[k-i-1] (1-indexed)
    # Equivalent (0-indexed): h_i 0-indexed (i'=i-1), w_total[j'] = wt[j'] (j'=k-i-1, so j' = k-1-1-i'=k-2-i'), so i'+j' = k-2.
    # T1_k = full[k-2]  for k ≥ 2; T1_1 = 0.
    T1 = np.empty(N)
    T1[0] = 0.0
    T1[1:] = full[: N - 1]
    # T2_k = sum_{j=1}^{N-k} w_total[j-1] h_{k+j}.  Let i = k+j (i = k+1..N), j = i-k.
    #   T2_k = sum_{i=k+1}^{N} h_i w_total[i-k-1]
    # 0-indexed: i' = i-1, j' = i-k-1 = i'-k.  So i' - j' = k.  T2_k = sum over (i',j') with i'-j' = k, h_{i'} wt[j'].
    # Cross-correlation. Use np.correlate(h, wt) which gives: for lag L,  (h ⋆ wt)[L] = sum_i h_i wt[i-L].
    cor_full = np.correlate(h_pad, wt_pad, mode="full")
    # cor_full[m] = sum_{i' - j' = m - (N-2)} h_{i'} wt[j'] for m = 0..2N-2.
    # numpy convention: np.correlate(a, v, mode='full') returns c[m] = sum_n a[n+m-len(v)+1] v[n]
    # That's a bit awkward; let's use a direct definition.
    #   cor[L] = sum_n h[n] wt[n - L] for L in [-(M-1), N-1] (where M = len(wt)=N-1)
    #   ranges: full has length N + M - 1 = 2N - 2, indexed 0..2N-3.  Lag L = m - (M-1) = m - N + 2.
    # We want T2_k = (lag = k) value of cor: i'-j' = k  →  L = k.
    #   m index of L=k: m = L + (M-1) = k + N - 2.
    T2 = np.empty(N)
    for k in range(1, N + 1):
        m = k + N - 2
        if 0 <= m < 2 * N - 2:
            T2[k - 1] = cor_full[m]
        else:
            T2[k - 1] = 0.0

    grad = Δ * (A + B - T1 - T2)
    return sm, grad
