# Unconditional Primality Certificates for $p = 3m(m+1)+1$
### Arithmetic Filters, Multi-Core Search, and Records to 29 998 Decimal Digits

**Author:** Hassane Bakkaoui — Independent Researcher, 2026  
**Companion papers:** [`paper1-qmin-bh/`](../paper1-qmin-bh/) · [`paper2-prime-families/`](../paper2-prime-families/) · [`paper3-padic-qmin/`](../paper3-padic-qmin/)  
**MSC 2020:** 11A41 · 11Y11 · 11N13

---

## Overview

This is **Paper IV** of the PrimeQuest project.

We study the parametric subfamily

```
p = 3m(m+1) + 1,   m = 2^a · 3^b − 1,   a, b ≥ 1
```

from the standpoint of *unconditional* primality certification.
The 3-smooth structure of `m+1 = 2^a · 3^b` makes the factor
`F = 2^a · 3^(b+1)` of `p−1` fully explicit and satisfies `F > √p`
for every `(a, b)` — reducing the Pocklington–Lehmer certificate
to a two-witness check (`q ∈ {2, 3}` only).

Three arithmetic theorems eliminate **87.1 %** of all candidates before
any large-number computation:

| Theorem | Statement | Saving |
|---------|-----------|--------|
| **Mod 6** | `p ≡ 1 (mod 6)` for all `(a, b)` | Primes 2, 3 excluded unconditionally |
| **Cubic reciprocity sieve** | Every prime `q ≡ 2 (mod 3)` satisfies `q ∤ p` | **−50 %** of trial divisions |
| **Forbidden classes mod 7** | `7 ∣ p ⟺ (a mod 3, b mod 6) ∈ F₇` | **−33 %** at zero arithmetic cost |

The **PrimeQuest** multi-core algorithm (v1–v7) applies this pipeline
and has produced four unconditional certificates, culminating in a prime
of **29 998 decimal digits** found in **3 h 06 min** on a single
consumer laptop (ARM64, 9 workers).

---

## Record — 29 998-Digit Prime

```
p = 3m(m+1) + 1,   m = 2^19435 · 3^19173 − 1   (14 999 digits)
p                                                (29 998 digits)

p = 1602087359509060348500… 16930852697403293697
```

### Pocklington–Lehmer Certificate

```
p − 1  =  F · m
F      =  2^19435 · 3^19174   (15 000 digits),   F > √p  ✓

q = 2,  w = 5 :   5^(p−1)     ≡ 1  (mod p)
                  gcd(5^((p−1)/2) − 1,  p)  =  1  ✓

q = 3,  w = 7 :   7^(p−1)     ≡ 1  (mod p)
                  gcd(7^((p−1)/3) − 1,  p)  =  1  ✓

Conclusion: p is prime  (deterministic proof).
```

---

## All Certified Primes

| Digits | `a` | `b` | `b/a` | Wall time | Workers |
|-------:|----:|----:|------:|----------:|--------:|
| 9 998 | 6 212 | 6 738 | 1.085 | — | 1 |
| 10 000 | 6 213 | 6 740 | 1.085 | — | 1 |
| 19 999 | 12 228 | 13 242 | 1.083 | 2 h 40 min | 9 |
| **29 998** | **19 435** | **19 173** | **0.987** | **3 h 06 min** | **9** |

---

## Algorithm Versions

| Version | Key optimisation | MR rounds | Sieve limit |
|---------|-----------------|-----------|-------------|
| v1 | Sequential baseline | 25 | 10⁶ |
| v2 | Sieve restricted to `q ≡ 1 (mod 3)` | 25 | 10⁶ |
| v3 | Multi-core · checkpoint · zigzag centre | 25 | 10⁶ |
| v4 | Mod-7 forbidden-class filter | 25 | 10⁶ |
| v5 | MR 25 → 3+20 rounds · sieve 10⁶ → 10⁷ | 3 + 20 | 10⁷ |
| v6 | Multi-ratio b/a simultaneous exploration | 3 + 20 | 10⁷ |
| v7 | Adaptive spread · tolerance ±50 | 3 + 20 | 10⁷ |

---

## Repository Structure

```
paper4-pocklington-certificates/
├── code/
│   ├── primequest_v1.py             Sequential baseline
│   ├── primequest_v2.py             Cubic sieve
│   ├── primequest_v3.py             Multi-core + checkpoint
│   ├── primequest_v4.py             + Mod-7 filter
│   ├── primequest_v5.py             + MR 3+20 · sieve 10M
│   ├── primequest_v6.py             + Multi-ratio
│   ├── primequest_v7.py             + Adaptive spread · tolerance ±50
│   └── migrer_checkpoint_v6_vers_v7.py
├── data/
│   ├── 01_premier_29998chiffres.txt  29 998 digits — current record
│   ├── 02_premier_19999chiffres.txt
│   ├── 03_premier_10000chiffres.txt
│   ├── 04_premier_9998chiffres_v2.txt
│   └── 05_premier_9998chiffres.txt
├── tex/
│   ├── paper_pocklington.tex         Full paper (English)
│   ├── paper_pocklington_fr.tex      Full paper (French)
│   └── Paper_IV_Pocklington_v2_text.tex
├── pdf/
│   └── Paper_IV_Pocklington_v2_0.pdf
├── figures/
│   └── fig_scalingpaper4.png
└── README.md
```

---

## Requirements

```bash
pip install gmpy2
```

Python ≥ 3.9, gmpy2 ≥ 2.1 (GMP backend).  
For optimal performance on ARM64 (Windows): install the native ARM64 build of Python and gmpy2.
