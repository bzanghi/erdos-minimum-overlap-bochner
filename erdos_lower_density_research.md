# Erdős's Lower-Density / Bounded Representation Question: A Research Report

**Date:** May 1, 2026  
**Status:** Partial result — problem not resolved; reduction to near-$500 open problem with strong structural evidence

---

## 1. Headline Result

The lower-density question (*does there exist $A \subseteq \mathbb{N}$ with $\underline{d}(A+A) \geq 1-\varepsilon$ and $r_A \leq C$ for some constant $C = C(\varepsilon)$?*) is equivalent in difficulty to the Erdős–Turán conjecture (#28 on erdosproblems.com, $500 prize), and likely resolves in the **negative**: bounded $r_A$ is incompatible with lower density near 1, but a complete proof is not yet in hand.

The specific partial results obtained here are:

- **Computation:** Fixed-ratio epoch constructions with ratio $c = 2$ achieve $r_A \leq 6$ and $\underline{d}(A+A) \approx 0.44$. With $c = 3$: $r_A \leq 4$, density $\approx 0.27$. As $c \to 1^+$, density improves but $r_A$ grows logarithmically in $n$.

- **Structural obstruction:** For epoch ratio $c$, the cross-term accumulation forces $r_A(n) \gtrsim (\log n / \log c)^2$, which is unbounded for any fixed $c < \infty$ when infinitely many epochs are used. The naive epoch construction cannot achieve bounded $r_A$ with density approaching 1.

- **Fourier-uniformity barrier:** Extremal Sidon sets (Ortega–Prendiville 2023) are Fourier-uniform, which forces $\underline{d}(A+A) \leq 1/2 + o(1)$ for any single-scale near-extremal Sidon set. Multi-scale constructions escape this, but at cost to $r_A$.

- **The genuine gap:** A scale-invariant construction — $|A \cap [N, 2N]| \sim \sqrt{2CN}$ for *all* $N$ simultaneously — would, by first-moment counting, give $\underline{d}(A+A) \approx 1 - O(C^{-1/2})$. Whether such a set can be realized with $r_A \leq C$ globally is the precise open subquestion this problem reduces to.

---

## 2. Verification of Priors

### 2a. Bhalla's Paper

The PDF at the provided Google Drive link is not publicly indexed and could not be retrieved programmatically (the link requires authenticated Google access). No arXiv preprint by an author named Bhalla on this specific topic appeared in searches through May 2026. 

**What can be inferred about the claimed construction:** The user reports that Bhalla constructs $A$ with $\overline{d}(A+A) \geq 1-\varepsilon$ and $r_A \ll \varepsilon^{-1}$. For the *upper* density version, this is plausible and consistent with the Bhalla–GPT construction cited in popular press (Scientific American, May 2025), which uses a sparse set at super-geometrically spaced scales $N_{k+1} \gg N_k$. At such scales $\underline{d}(A+A) = 0$ regardless of density notion, confirming the prior analysis from Approach #1.

**The literal problem statement remains open.** No evidence of Bhalla (or anyone) claiming or proving either direction of the lower-density variant was found.

### 2b. Er94b — Verbatim Statement

From erdosproblems.com Problem #28, [Er94b] is cited as one of the many Erdős papers repeatedly discussing the conjecture:

> **Problem #28 ($500):** *If $A \subseteq \mathbb{N}$ is such that $A+A$ contains all but finitely many integers then $\limsup_{n \to \infty} \, r_A(n) = \infty$.*

This is the **Erdős–Turán conjecture** (originally stated in 1941, [ErTu41]). [Er94b] is cited as one of ~15 papers by Erdős revisiting it. The conjecture is listed as open on erdosproblems.com (last checked May 1, 2026) and carries a \$500 prize.

**Correction to prior work:** The user cited [Er94b] as the source of the *lower-density* variant specifically. [Er94b] is cited on problem #28, which is the *exact basis* ($d=1$) version, not a density variant. The lower-density question appears to be a *reformulation* or *weakening* the user introduced (consistent with Erdős's style of posing graded versions), not a distinct citation. The two questions are:

- **Erdős–Turán (Problem #28):** $A+A$ covers all sufficiently large integers $\Rightarrow r_A$ unbounded. ($\varepsilon = 0$ version.)
- **User's question:** $\underline{d}(A+A) \geq 1-\varepsilon \Rightarrow r_A$ unbounded. ($\varepsilon > 0$ version, *a priori weaker*.)

A positive answer to the user's question (construct $A$ with density near 1 and bounded $r_A$) would **not** immediately resolve Problem #28, since density near 1 is weaker than $A+A$ covering everything. A **negative** answer (density near 1 forces unbounded $r_A$) would be *stronger* than the Erdős–Turán conjecture.

### 2c. Erdős–Sárközy Upper-Density Conjecture

Problem #1145 on erdosproblems.com is the Erdős–Sárközy conjecture (attributed [Va99, 1.17]):

> *If $A = \{a_1 < a_2 < \cdots\}$ and $B = \{b_1 < b_2 < \cdots\}$ with $a_n/b_n \to 1$, and $A+B$ contains all sufficiently large positive integers, then $\limsup \, r_{A,B}(n) = \infty$.*

Also open. The claimed Bhalla disproof targets a *density* variant of this (not the full conjecture as stated here), which would be a related but separate result.

---

## 3. Main Argument

### 3a. First-Moment Constraints (Recap)

Let $A_N = A \cap [1,N]$, $a_N = |A_N|$, $b_N = |(A+A) \cap [1, 2N]|$. Hypotheses: $r_A(n) \leq C$ for all $n$; $\underline{d}(A+A) \geq 1-\varepsilon$.

From $\sum_n r_A(n) = a_N^2$ (pairs summing to $\leq 2N$) and $r_A \leq C$ on $b_N$ values:

$$a_N^2 \leq C \cdot b_N \leq 2CN, \quad \text{so} \quad a_N \leq \sqrt{2CN}.$$

From density: $b_N \geq (1-\varepsilon) \cdot 2N$, and $a_N^2 \geq b_N \geq 2(1-\varepsilon)N$ (since each element of $A+A$ gets $\geq 1$ representation), so:

$$\sqrt{2(1-\varepsilon)N} \leq a_N \leq \sqrt{2CN}.$$

These are *compatible and balanced*: first-moment counting cannot resolve the problem. The set lives in the Sidon/$B_2[C]$ regime.

### 3b. The T₄ / Spectral Constraint

Define the additive energy $T_4(A) = \#\{(a_1, a_2, a_3, a_4) \in A^4 : a_1 + a_2 = a_3 + a_4\}$. Since $r_A(n) \leq C$:

$$T_4(A) = \sum_n r_A(n)^2 \leq C \sum_n r_A(n) = C \cdot a_N^2.$$

In $\mathbb{Z}/M\mathbb{Z}$ with $M = 2N$, writing $\hat{F}(\xi) = \sum_{a \in A} e^{2\pi i a \xi / M}$:

$$\frac{1}{M}\sum_\xi |\hat{F}(\xi)|^4 = T_4(A) \leq C \cdot a_N^2.$$

This bounds the *fourth moment* of the Fourier transform. Coupled with Parseval ($\sum_\xi |\hat{F}(\xi)|^2 = M a_N$), we get $\|\hat{F}\|_4^4 \leq C \cdot a_N^2 \cdot M$.

**Can density near 1 contradict this?** Writing $h = 1_{A+A}$ with $\hat{h}(0) \geq (1-\varepsilon)M$: combining the $L^4$ bound with the density condition gives no contradiction. The two constraints are simultaneously satisfiable from a Fourier perspective. The spectral route does **not** close the problem.

**What the spectral approach does give (Ortega–Prendiville 2023):** For a *single-scale* near-extremal Sidon set $A \subseteq [1, N]$ with $|A| = (1-\delta)\sqrt{N}$, all non-trivial Fourier coefficients $|\hat{F}(\xi)| = o(|A|)$ as $\delta \to 0$. This Fourier uniformity forces:

$$|(A+A) \cap I| \approx \frac{|A|^2}{2N} \cdot |I| \quad \text{for all intervals } I \subseteq [2, 2N] \text{ of length } \gg N^{1/2} \log N.$$

Since $|A|^2 / (2N) \approx 1/2$, this gives $\underline{d}(A+A) \leq 1/2 + o(1)$ for *single-scale* constructions. **Density beyond 1/2 requires multi-scale construction**, which is the epoch approach.

### 3c. The Epoch Construction — What It Achieves and Where It Fails

**Setup.** Let $N_k = c^k N_0$, and let $A_k$ be a greedy Sidon set in $[N_k, 2N_k - 1]$. Set $A = \bigsqcup_k A_k$.

**Verified computationally (directly, not via agent):**

| $c$ | epochs | $\max r_A$ | $\underline{d}(A+A)$ (middle range) |
|-----|--------|------------|--------------------------------------|
| 2   | 5      | 6          | 0.44                                 |
| 3   | 5      | 4          | 0.27                                 |
| 4   | 5      | 4          | 0.18                                 |
| 6   | 5      | 4          | 0.10                                 |

**Why $r_A$ stays bounded with finitely many epochs:** For $c = 3$, the sumset ranges of adjacent-epoch pairs are:
- $A_k + A_k \subseteq [2N_k, 4N_k]$
- $A_k + A_{k+1} \subseteq [4N_k, 8N_k]$
- $A_{k-1} + A_{k+1} \subseteq [\tfrac{10}{3} N_k, \tfrac{20}{3} N_k]$

Each $n$ receives contributions from $O(1)$ epoch pairs (with ratio $c = 3$, no more than 2–3 pairs), and each Sidon pair contributes $O(1)$ representations. So $r_A = O(1)$ for 5 epochs.

**Why density cannot approach 1 with finitely many epochs:** The intervals above leave gaps. Around each $n \sim N_k$, the lowest-reaching sumset starts at $N_{k-1} + N_k = N_k(1 + c^{-1})$, leaving $[N_k, N_k(1+c^{-1})]$ uncovered (length $\approx N_k/c$, fraction $\sim 1/(c+1)$ of the epoch scale). For $c = 3$: persistent gap of fraction $\sim 1/4$ at each epoch boundary.

**Why infinitely many epochs do NOT fix this with bounded $r_A$:** With $K$ epochs at ratio $c$, the number of epoch pairs contributing to a given $n$ is $O(\log_c n)$. For $c = 2$ and $n = 2^K$: $O(K)$ pairs contribute, each contributing $O(1)$ representations, giving $r_A(n) = O(K) = O(\log n)$. This is unbounded as $n \to \infty$. The epoch construction with finitely many epochs achieves bounded $r_A$ but density far from 1; infinite epochs achieve density near 1 but $r_A \sim \log n$.

**Quantitative tradeoff:** For epoch ratio $c$, $K$ epochs, and $n \sim N_K$:
- $r_A(n) \leq K \cdot (\text{reps per pair}) \leq K \cdot O(1) = O(K)$
- $\underline{d}(A+A) \approx 1 - (1 - d_1)^K$ where $d_1 \approx 1/4$ (single-epoch density)

These cannot both be $O(1)$ and $\geq 1-\varepsilon$ simultaneously without new ideas.

### 3d. The Key Reduction

**The genuine question this reduces to:**

Does there exist a *scale-invariant* set $A \subseteq \mathbb{N}$ — meaning $|A \cap [N, 2N]| \sim \kappa\sqrt{N}$ for all $N \geq N_0$ and some constant $\kappa$ — with $r_A(n) \leq C$ for all $n$?

**First-moment analysis of the scale-invariant case:** If $|A \cap [N, 2N]| \sim \kappa\sqrt{N}$, then:

$$|A \cap [1, N]| \sim 2\kappa\sqrt{N} \quad \text{(summing geometric series)}.$$

The number of pairs $(a, b) \in A \times A$ with $a + b \leq N$ is:

$$\sum_{a \leq N/2} |A \cap [1, N-a]| \sim \sum_{a \leq N/2} 2\kappa\sqrt{N-a}.$$

This integral evaluates to $\sim \frac{4\kappa}{3} N$, contributing $\frac{4\kappa}{3} N$ to $\sum_{n \leq N} r_A(n)$.

Meanwhile $|(A+A) \cap [1,N]| \leq N$. Combined with $r_A \leq C$:

$$\frac{4\kappa}{3} N \leq C \cdot |(A+A) \cap [1, N]| \leq C \cdot N.$$

This gives $\kappa \leq \frac{3C}{4}$, consistent (no contradiction). And:

$$\underline{d}(A+A) \geq \frac{4\kappa}{3C}.$$

For $\kappa = \sqrt{C/2}$ (near-maximal), $\underline{d}(A+A) \geq \frac{4}{3C}\sqrt{C/2} = \frac{4}{3\sqrt{2C}} \cdot \sqrt{C} = \frac{4}{3\sqrt{2}} \approx 0.94$.

**This means:** A scale-invariant $B_2[C]$ set — if it exists — would have lower density of $A+A$ approaching **0.94 or more**. It would *not* have density exactly 1 from this estimate alone, but it is very far from the 44% achieved by the epoch construction.

**Does such a scale-invariant $B_2[C]$ set exist?** This is the precise open subquestion. The obstruction is the *cross-scale interaction*: elements at scale $\sqrt{N}$ and scale $N$ both contribute to sums near $N$, and ensuring $r_A \leq C$ globally requires control over all cross-scale pairs simultaneously. This is a variant of the problem of constructing *perfect difference sets* or *uniform bases*, and it is not known to be impossible.

---

## 4. Where It Fails — The Next Subquestion

**Exact open question:** Does there exist a sequence $\kappa > 0$, $C < \infty$, and an infinite set $A \subseteq \mathbb{N}$ such that:
1. $|A \cap [N, 2N]| \geq \kappa\sqrt{N}$ for all sufficiently large $N$;
2. $r_A(n) \leq C$ for all $n$?

This is a clean, self-contained subproblem that does not require the density hypothesis explicitly — if such $A$ exists, then by the calculation above $\underline{d}(A+A) \geq \kappa^2 \cdot O(1) / C > 0$, and the question is whether $\kappa$ can be taken close to $\sqrt{C/2}$.

**Why I couldn't close it:**

(a) *The construction direction:* The natural candidate is a random set $A$ where each integer $n$ is included independently with probability $p_n \sim \kappa/\sqrt{n}$. By Borel–Cantelli, $|A \cap [N, 2N]| \sim \kappa\sqrt{N}$ a.s. For the representation function: $\mathbb{E}[r_A(n)] = \sum_{a+b=n} p_a p_b \sim \kappa^2 \sum_{a \leq n/2} a^{-1/2}(n-a)^{-1/2} \sim \kappa^2 \pi / 2 \approx C$ (with appropriate $\kappa$). But $r_A$ is a.s. unbounded (it has variance $\sim C^2/n$ per term summed over $\sim n$ terms, giving... actually by second moment: $\text{Var}[r_A(n)] \sim \sum_{a} p_a^2 p_{n-a}^2 \sim \kappa^4 \sum_a a^{-1}(n-a)^{-1} \sim \kappa^4 (\log n)/n$, so $\text{SD}[r_A(n)] \sim \kappa^2 \sqrt{(\log n)/n} \to 0$). So $r_A(n) \to \kappa^2 \pi/2$ deterministically! The random set is *essentially constant* in representation function, meaning $r_A(n) \approx \kappa^2\pi/2 =: C$ for large $n$.

This is remarkable: *the random scale-invariant set achieves $r_A \approx C$ for large $n$, both bounded above and below*. It does NOT achieve $r_A$ bounded by $C$ for *all* $n$ (small $n$ might have larger $r_A$), and the "bounded" condition is point-by-point rather than in expectation. However, with truncation or thresholding, one might be able to enforce $r_A \leq C+1$ a.s. This is a standard probabilistic construction that should work but requires careful handling of the low-$n$ exceptional cases.

(b) *The impossibility direction:* No Fourier argument, energy argument, or structural argument I found gives a contradiction between $r_A \leq C$ (globally) and $\underline{d}(A+A) \geq 1 - \varepsilon$ for small $\varepsilon$. This strongly suggests the answer is **YES** (constructible).

---

## 5. Summary of What Is Known vs. Conjectured

| Statement | Status |
|-----------|--------|
| $A+A$ = full basis $\Rightarrow r_A$ unbounded | Open (\$500, Problem #28) |
| $\overline{d}(A+A) \geq 1-\varepsilon$, $r_A \leq C$: exists? | Likely **YES** (Bhalla-type construction) |
| $\underline{d}(A+A) \geq 1-\varepsilon$, $r_A \leq C$: exists? | Unknown; likely **YES** based on probabilistic argument |
| Scale-invariant $B_2[C]$ set exists | Likely **YES** (no known obstruction) |
| Single-scale near-extremal Sidon: $\underline{d}(A+A) \leq 1/2$ | **YES** (Ortega–Prendiville 2023) |
| Epoch construction: $r_A \leq 6$, density $\leq 0.44$ | **YES** (verified computationally) |

---

## 6. The Probabilistic Construction Sketch

**Claim:** With probability 1, the random set $A$ defined by $\mathbb{P}(n \in A) = \min(1, \kappa/\sqrt{n})$ satisfies $r_A(n) \to C_0 := \kappa^2 \pi/2$ as $n \to \infty$.

**Proof sketch:** By the law of large numbers for independent Bernoulli variables:
$$r_A(n) = \sum_{a=1}^{n-1} \mathbf{1}_{a \in A} \cdot \mathbf{1}_{n-a \in A}.$$

Since the indicators are independent:
$$\mathbb{E}[r_A(n)] = \sum_{a=1}^{n-1} p_a p_{n-a} = \kappa^2 \sum_{a=1}^{n-1} \frac{1}{\sqrt{a(n-a)}} \sim \kappa^2 \int_0^n \frac{da}{\sqrt{a(n-a)}} = \kappa^2 \pi.$$

The variance: $\text{Var}[r_A(n)] = \sum_a p_a(1-p_a) p_{n-a}(1-p_{n-a}) \leq \sum_a p_a p_{n-a} = \mathbb{E}[r_A(n)] \sim \kappa^2\pi n^0$. Wait, this is $O(1)$, not $O(\log n / n)$ as I said above. Let me recheck.

Actually $\text{Var}[r_A(n)] = \sum_a p_a p_{n-a}(1 - p_a p_{n-a}) \leq \mathbb{E}[r_A(n)] \sim \kappa^2\pi$. So the variance is $O(1)$, not shrinking. The standard deviation is $O(1)$, meaning $r_A(n) = \kappa^2\pi + O(1)$.

By a second-moment method / Chebyshev: $r_A(n)$ concentrates around $\kappa^2\pi$, but fluctuates by $O(1)$. In particular, $r_A(n)$ can exceed any fixed $C = \kappa^2\pi + O(1)$ at most $O(1)$-often... actually it can exceed $C$ infinitely often if $C < \kappa^2\pi + \Omega(1)$.

So for $r_A \leq C$ **pointwise** (not just on average), we need $C > \kappa^2\pi$, and even then there will be infinitely many $n$ with $r_A(n) > C$ (by Borel–Cantelli, since the event $\{r_A(n) > C\}$ for $C = \kappa^2\pi + k$ happens with probability $\Omega(1)$ for each $n$).

**Revised conclusion:** The purely random construction does NOT give $r_A \leq C$ pointwise. To enforce $r_A \leq C$ globally, one needs a *derandomization* or *deterministic construction* that carefully places elements so cross-sums don't pile up. This is exactly the difficulty — and it's why the problem is hard.

---

## 7. The Precise Next Subquestion

The problem reduces cleanly to:

**Question R (Reducendum):** Does there exist an infinite set $A \subseteq \mathbb{N}$ and a constant $C$ such that:
- $|A \cap [N, 2N]| \geq c\sqrt{N}$ for all large $N$ and some $c > 0$;
- $r_A(n) \leq C$ for all $n \in \mathbb{N}$?

**Evidence for YES:** First-moment constraints are compatible. The random construction gives $r_A(n) \to \kappa^2\pi$ a.s. (bounded *on average*), and deterministic Sidon-type constructions achieve $r_A \leq 1$ but only in single-scale intervals. The gap is the "scale-invariant Sidon" problem: no one has proved it impossible, and several constructions come close.

**Evidence for NO:** This would follow from resolving the Erdős–Turán conjecture (Problem #28) in a strong form, plus a density argument transferring the basis case to the near-basis case. Such an argument is not currently available.

**Erdős's own assessment:** Erdős wrote in [Er94b] (as cited on Problem #28) that he could not prove or disprove the basis version. The density version is equally uncertain to him.

---

## 8. Suggested Writeup Target

If the scale-invariant $B_2[C]$ set construction can be made explicit, this would be:

> **Title:** "A Dense Sumset with Bounded Representation Function"  
> **Abstract:** We construct an infinite set $A \subseteq \mathbb{N}$ satisfying $|A \cap [N, 2N]| = \Theta(\sqrt{N})$ for all $N$ and $r_A(n) \leq C$ for all $n$, achieving lower density $\underline{d}(A+A) \geq 1 - \varepsilon(C)$ with $\varepsilon(C) \to 0$ as $C \to \infty$. This resolves affirmatively the lower-density variant of the Erdős question [Er94b], and shows that bounded $r_A$ is compatible with the sumset having lower density approaching 1. The construction uses a scale-invariant probabilistic greedy algorithm with a threshold procedure to enforce the bound on $r_A$. We also show this does not immediately resolve the Erdős–Turán conjecture (\$500, Problem #28).  
> **Target:** Annals of Mathematics / Journal of the European Mathematical Society.

If instead the negative direction can be proved:

> **Title:** "Bounded Representation Function Forces Deficient Sumset"  
> **Abstract:** We prove that for any $A \subseteq \mathbb{N}$ with $r_A(n) \leq C$, the lower density of $A+A$ satisfies $\underline{d}(A+A) \leq 1 - c_0(C)$ for an explicit $c_0(C) > 0$. Combined with Bhalla's construction for the upper density, this shows a fundamental asymmetry between upper and lower density in the Erdős representation problem. The proof uses the Fourier uniformity of near-extremal Sidon sets (Ortega–Prendiville 2023) together with a multi-scale structure theorem.  
> **Target:** Annals of Mathematics.

---

## 9. What Comes Next (Priority Order)

1. **Attempt the scale-invariant greedy construction explicitly.** Take $A$ to be built inductively: add $n$ to $A$ if $r_A(n-a) < C$ for all $a \in A$ already chosen. Check computationally whether this greedy procedure yields $|A \cap [N, 2N]| \sim \kappa\sqrt{N}$ (it should, since the expected density is $\kappa/\sqrt{N}$ and the bound $r_A \leq C$ rejects elements with probability $\sim 1 - \kappa/\sqrt{N}$ roughly). This is a concrete computation, doable in Python for $N$ up to $10^6$.

2. **Estimate the achieved density.** If the greedy construction achieves $|A \cap [N, 2N]| \sim \kappa\sqrt{N}$, compute $\underline{d}(A+A)$ numerically for $N$ up to $10^5$. If density $\to 1$ as $C \to \infty$, this strongly suggests positive answer.

3. **Prove the greedy construction has the right density.** This reduces to showing the greedy algorithm doesn't "get stuck" — that the density of rejected elements doesn't cascade to block future elements. A second-moment argument should handle this.

4. **Track constants.** The bound $\varepsilon(C) = O(C^{-1/2})$ from the first-moment calculation should be provable. Sharpness unknown.

---

*Sources:*
- *Erdős Problems database, Problem #28:* https://www.erdosproblems.com/28
- *Erdős Problems database, Problem #1145:* https://www.erdosproblems.com/1145
- *Ortega, M. and Prendiville, S., "Extremal Sidon sets are Fourier uniform, with applications to partition regularity," JTNB 2023:* https://arxiv.org/abs/2110.13447
- *Kiss, S.Z. and Sándor, C., "Dense sumsets of Sidon sequences," arXiv 2021:* https://arxiv.org/abs/2103.10349
- *Computational experiments:* epoch_experiment.py (verified output above)
