# Paper 3 — Unconditional Primality Certificates

**Title:** Unconditional Primality Certificates for Numbers of the Form p = 3m(m+1)+1:
Arithmetic Filters, Multi-Core Search, and Records to 29 998 Decimal Digits

**Author:** Hassane Bakkaoui — 2026

## Contents

```
paper3-pocklington-certificates/
├── tex/
│   └── paper_pocklington.tex     ← Full paper (LaTeX)
├── code/
│   ├── primequest_v1.py          ← Sequential search (baseline)
│   ├── primequest_v2.py          ← Sieve restricted to q≡1(mod3) [Thm cubic]
│   ├── primequest_v3.py          ← Multi-core + checkpoint + zigzag
│   ├── primequest_v4.py          ← Mod-7 filter [Thm mod7], ratio=1.0
│   └── primequest_v5.py          ← MR 3+20 rounds, sieve limit 10M (50k target)
└── data/
    ├── 01_premier_29998chiffres.txt   ← 29 998 digits (a=19435, b=19173) — record
    ├── 02_premier_9998chiffres.txt    ← 9 998 digits (first certificate)
    ├── 03_premier_9998chiffres_v2.txt ← 9 998 digits (v2 variant)
    ├── 04_premier_10000chiffres.txt   ← 10 000 digits
    └── 05_premier_19999chiffres.txt   ← 19 999 digits
```

## Key results

| Digits | a | b | b/a | MR tests | Wall time | Workers |
|--------|---|---|-----|----------|-----------|---------|
| 9 998 | 6 212 | 6 738 | 1.085 | — | — | 1 |
| 10 000 | 6 213 | 6 740 | 1.085 | — | — | 1 |
| 19 999 | 12 228 | 13 242 | 1.083 | 1 050 | 2h40 | 9 |
| **29 998** | **19 435** | **19 173** | **0.987** | **286** | **3h06** | **9** |

## Certificate (29 998 digits)

```
p − 1 = 2^19435 · 3^19174 · m,   F = 2^19435 · 3^19174 > √p
q=2  →  w=5 :  5^(p−1) ≡ 1 (mod p),  gcd(5^((p−1)/2)−1, p) = 1
q=3  →  w=7 :  7^(p−1) ≡ 1 (mod p),  gcd(7^((p−1)/3)−1, p) = 1
```

## Pre-filtering: 87.1% eliminated before Miller-Rabin

| Stage | Eliminated | % |
|-------|-----------|---|
| Mod-7 filter (Theorem mod7) | 740 | 33.3% |
| Sieve q≡1(mod3) | 1 195 | 53.8% |
| Miller-Rabin (composites) | 285 | 12.8% |
