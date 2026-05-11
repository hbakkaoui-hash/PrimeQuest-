# $p$-adic Distribution of the Minimal Offset $q_{\min}$
### Bateman–Horn Predictions, Conditional Theorems, and an Anti-Persistence Sign-Reversal Phenomenon

**Author:** Hassane Bakkaoui — Independent Researcher, 2026  
**Status:** Working draft v1.0 (May 2026)  
**Companion papers:** [`paper1-qmin-bh/`](../paper1-qmin-bh/) · [`paper2-prime-families/`](../paper2-prime-families/) · [`paper4-pocklington-certificates/`](../paper4-pocklington-certificates/)  
**MSC 2020:** 11N32 · 11N13 · 11Y11 · 11M41 · 62P99

---

## Overview

This is **Paper III** of the PrimeQuest project — the third part of a trilogy on the parametric prime family

```
p = k·m(m+1) + ε + 2k·q,   k ≥ 1,  gcd(ε, 2k) = 1
```

We study the **joint modular distribution** of the minimal offset $q_{\min}(k, m, \varepsilon)$ modulo small primes $\ell$, completing the picture established in Papers I and II:

| Paper | Contribution |
|-------|-------------|
| **I** (Apr 2026) | Marginal law: $\mathbb{E}[\|q_{\min}\|] \sim \log m / C_k$ |
| **II** (Apr 2026) | Conditional law: $\mathbb{E}[\|q\| \mid m] = m/4$; exclusion of spurious $\zeta(s)$ connection |
| **III** (May 2026) | Joint modular layer — third axis of the description |

---

## Principal Contributions

### (i) Forbidden Residue
An explicit **forbidden residue** $q^*(k, m, \varepsilon; \ell) \in \mathbb{Z}/\ell\mathbb{Z}$ (Proposition 1), giving rise to the **permission frequency** $\Pi_P(r; k, \ell)$.

### (ii) Conditional Theorem (Theorem 1)
Under hypotheses H1′ and H3′:

$$F_N(r;\, k,\ell) = \frac{\Pi_P(r;\, k,\ell)}{\ell - 1} + O(N^{-1/2})$$

Validated empirically at $N = 10^6$ for $k \in \{3, 7, 15, 21\}$ — mean Pearson correlation prediction/observation $r = 0.99$.

### (iii) Unconditional Theorem (Theorem 2)
For $\ell \mid 2k$:

$$F_N(r) = \frac{1}{\ell} + O(N^{-1/2}) \quad \text{(unconditionally)}$$

Verified across all 10 tested cases with $\chi^2 < 3.1$.

### (iv) Anti-Persistence Sign-Reversal (New Phenomenon)
The lag-1 autocorrelation of $(q_n \bmod 2)$ follows a **linear law in $k$**:

```
ACF₁ ≈ 0.024·k − 0.20   (R² = 0.97,  k ∈ [3, 21],  N = 10⁶)
```

Sign-reversal threshold at $k_0 \approx 8.4$. The full spectral signature is fitted by a **damped cosine** (period $T \approx 5.89$, damping $\gamma \approx 0.39$, $R^2 = 0.999$). At $k=21$, mod 5: $\mathrm{ACF}_1 = +0.64$.

---

## Status of Results

| Result | Status |
|--------|--------|
| Theorem 2 (ℓ ∣ 2k case) | **Unconditional** |
| Theorem 1 (general case) | Conditional on H1′, H3′ |
| ACF conjectures | Empirical / heuristic |

The unconditional barrier (parity problem; uniform PNT in arithmetic progressions of large moduli) is discussed in §6.

---

## Repository Structure

```
paper3-padic-qmin/
├── tex/
│   └── Paper_III_v1_0_FINALbi.tex    Full paper — LaTeX source
├── pdf/
│   └── Paper_III_v1_0_FINALb.pdf     Compiled PDF
└── README.md
```
