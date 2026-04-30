# PrimeQuest — Deux papiers sur les nombres premiers

Dépôt de recherche de Hassane BAKKAOUI.

## Structure

```
PrimeQuest-/
├── paper1-qmin-bh/          — Papier I  : famille paramétrique & constante BH
└── paper2-pocklington-primes/ — Papier II : primes de Pocklington p=3m(m+1)+1
```

---

## Papier I — Famille paramétrique de premiers & analyse BH

**Dossier :** `paper1-qmin-bh/`

Étude de la distribution de q_min dans une famille paramétrique de nombres
premiers, et vérification numérique de la conjecture de Bateman-Horn.

- `tex/main.tex`                  — source LaTeX
- `Parametric_Family_of_Primes.pdf` — article complet (PDF)
- `code/`                         — scripts Python d'analyse
- `data/TABLE_full_subsample_10e8.csv` — données numériques

---

## Papier II — Primes de Pocklington : p = 3m(m+1)+1, m = 2^a·3^b−1

**Dossier :** `paper2-pocklington-primes/`

Construction, preuve de primalité (Théorème de Pocklington N-1) et analyse
arithmétique d'une famille de grands nombres premiers.

- `tex/paper_II_revised.tex`      — source LaTeX
- `code/primequest.py`            — PrimeQuest v3 (multi-cœurs, ARM64)
- `data/`                         — certificats de primalité produits
- Théorèmes prouvés : p≡1(mod 6), classes interdites mod 7, q≡1(mod 3)
- Résultat : premier de **19 999 chiffres** prouvé sur PC ordinaire

---

## Installation

```bash
pip install gmpy2 matplotlib pandas
```

ARM64 (Snapdragon / Apple Silicon) : utiliser Python natif ARM64.
