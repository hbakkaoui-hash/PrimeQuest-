# Unconditional Primality Certificates for $p = 3m(m+1)+1$
### Arithmetic Filters, Multi-Core Search, and Records to 29 998 Decimal Digits

**Author:** Hassane Bakkaoui — Independent Researcher, 2026  
**Companion paper:** [`paper2-prime-families/`](../paper2-prime-families/) · [`paper1-qmin-bh/`](../paper1-qmin-bh/)  
**MSC 2020:** 11A41 · 11Y11 · 11N13

---

## Overview

This repository contains the full computational and theoretical material for
**Paper III** of the PrimeQuest project.

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

The **PrimeQuest** multi-core algorithm (v1–v5) applies this pipeline
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

| Digits | `a` | `b` | `b/a` | MR tests | Pre-filtered | Wall time | Workers |
|-------:|----:|----:|------:|---------:|-------------:|----------:|--------:|
| 9 998 | 6 212 | 6 738 | 1.085 | — | — | — | 1 |
| 10 000 | 6 213 | 6 740 | 1.085 | — | — | — | 1 |
| 19 999 | 12 228 | 13 242 | 1.083 | 1 050 | 87.9 % | 2 h 40 min | 9 |
| **29 998** | **19 435** | **19 173** | **0.987** | **286** | **87.1 %** | **3 h 06 min** | **9** |

---

## Filter Pipeline — 29 998-Digit Record

```
2 221 candidate pairs explored
│
├─  740  (33.3 %)  ── Mod-7 filter  [O(1) lookup — no bignum arithmetic]
├─ 1195  (53.8 %)  ── Sieve mod q≡1(mod 3)  [small modular exponentiations]
└─  285  (12.8 %)  ── Miller–Rabin  [composites surviving the sieve]
     │
     └─    1  (0.045 %)  ── PRIME — Pocklington certificate issued
```

**Total pre-MR elimination: 87.1 %**

---

## Algorithm Versions

| Version | Key optimisation | MR rounds | Sieve limit |
|---------|-----------------|-----------|-------------|
| v1 | Sequential baseline | 25 | 10⁶ |
| v2 | Sieve restricted to `q ≡ 1 (mod 3)` | 25 | 10⁶ |
| v3 | Multi-core · checkpoint · zigzag centre | 25 | 10⁶ |
| v4 | Mod-7 forbidden-class filter | 25 | 10⁶ |
| v5 | MR 25 → 3+20 rounds · sieve 10⁶ → 10⁷ | 3 + 20 | 10⁷ |

> **v5 note:** Miller–Rabin is deterministic for primes — every prime
> passes every round with certainty. Reducing from 25 to 3 initial rounds
> yields an **~8× speedup** with zero risk of missing a true prime.
> Composites surviving 3 rounds (≤ 4⁻³ ≈ 1.6 %) are caught by
> 20 confirmation rounds and the Pocklington test.

---

## Repository Structure

```
paper3-pocklington-certificates/
├── tex/
│   └── paper_pocklington.tex        Full paper — LaTeX, ready to compile
├── code/
│   ├── primequest_v1.py             Sequential search (baseline)
│   ├── primequest_v2.py             Cubic sieve (Theorem 3)
│   ├── primequest_v3.py             Multi-core + checkpoint + zigzag
│   ├── primequest_v4.py             + Mod-7 filter (Theorem 2), ratio = 1.0
│   └── primequest_v5.py             + MR 3+20 rounds, sieve 10M  [50k target]
└── data/
    ├── 01_premier_29998chiffres.txt  29 998 digits — current record
    ├── 02_premier_19999chiffres.txt  19 999 digits
    ├── 03_premier_10000chiffres.txt  10 000 digits
    ├── 04_premier_9998chiffres_v2.txt
    └── 05_premier_9998chiffres.txt   9 998 digits
```

---

## Requirements

```bash
pip install gmpy2
```

Python ≥ 3.9, gmpy2 ≥ 2.1 (GMP backend).  
For optimal performance on ARM64 (Windows): install the native ARM64 build
of Python and gmpy2 — running under x86 emulation (Prism) is ~4.6× slower.

---

## Running the search

```bash
# Edit DIGITS_CIBLE, TIMEOUT_S at the top of the script, then:
python primequest_v4.py     # 30 000-digit search
python primequest_v5.py     # 50 000-digit search  (in progress)
```

Checkpoint files are saved every 60 s; restart the same command to resume.

---

## Hardware (record run)

Dell Inspiron 14 Plus 7441 — Snapdragon X Plus (ARM64, 10 Oryon cores, LPDDR5X)  
9 parallel workers · 11 146 s (3 h 06 min) total wall time
