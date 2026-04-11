# Parametric Primes 2026

**Code source associé au manuscrit :**

> Hassane BAKKAOUI (2026).  
> *Une famille paramétrique de nombres premiers via les formes quadratiques —
> Distribution heuristique du décalage minimal, validation numérique à grande
> échelle et contrôle statistique des corrélations spectrales spurieuses.*  
> arXiv : math.NT [à compléter] · MSC 2020 : 11N32, 11Y11, 11M41, 94A60, 62P99

---

## Famille étudiée

$$p_{k,m,\varepsilon,q} = k \cdot m(m+1) + \varepsilon + 2kq, \qquad
k \in \{3, 15, 21\},\; \varepsilon \in \{\pm1\},\; q \in \mathbb{Z}$$

Le **décalage minimal** $q_{\min}(k,m)$ est le plus petit $|q|$ tel que
$p_{k,m,\varepsilon,q}$ soit premier.

**Loi heuristique principale (conditionnelle à Bateman–Horn + H2') :**

$$\mathbb{E}[|q_{\min}(k,m)|] \;\sim\; \frac{\ln m}{C_k}, \qquad m \to \infty$$

---

## Structure du dépôt

```
parametric-primes-2026/
├── README.md                  ← ce fichier
├── LICENSE                    ← MIT License
├── requirements.txt           ← dépendances Python
│
├── bh_constants.py            ← (1) Constantes de Bateman–Horn
├── qmin_compute.py            ← (2) Calcul de q_min sur 10^6 valeurs
├── distrib_validation.py      ← (3) Tests KS, bootstrap, chi²
├── permutation_tests.py       ← (4) Tests de permutation T1/T2/T3
└── Q_decomposition.py         ← (5) Décomposition Q ≈ α√p + β·D_N
```

---

## Description des scripts

### 1. `bh_constants.py` — Constantes de Bateman–Horn

Calcule $C(f_{k,\varepsilon})$ pour les six paires $(k,\varepsilon)$
par produit eulérien tronqué à $p \leq 200\,000$ (17 984 premiers).

```bash
python bh_constants.py
```

**Résultats reproduits (Table 2 du papier) :**

| $k$ | $C(f_{k,+})$ | $C(f_{k,-})$ | $C_k$    | $1/C_k$  |
|-----|-------------|-------------|---------|---------|
| 3   | 3.361048    | 2.783543    | 6.14459 | 0.16274 |
| 15  | 3.034936    | 3.778642    | 6.81358 | 0.14677 |
| 21  | 3.876345    | 4.329956    | 8.20630 | 0.12186 |

---

### 2. `qmin_compute.py` — Calcul de $q_{\min}$

Calcule $q_{\min}(21, m)$ pour $m \in [0, M]$, ajuste la loi logarithmique
$\mathbb{E}[|q|] \approx a \ln m + b$ et sauvegarde les résultats en `.npy`.

```bash
python qmin_compute.py --quick      # M=10^4, rapide
python qmin_compute.py              # M=10^6, complet (~5h en Python pur)
```

**Résultats reproduits (Observation 5.1) :**
- Pente = 0.320 · [IC 95 % : 0.31, 0.33], $R^2 = 0.984$
- $\rho_{21} = 0.320 \times C_{21} = 2.626 \pm 0.08$

---

### 3. `distrib_validation.py` — Validation distributionnelle

Teste la Conjecture 3.2 : $|q_{\min}| \cdot C_{21} / \ln m \xrightarrow{d} \mathrm{Exp}(1)$.

```bash
python distrib_validation.py --quick    # données synthétiques
python distrib_validation.py            # nécessite qmin_k21_M*.npy
```

**Trois tests (Section 5.4) :**

| Test | Statistique | Résultat papier |
|------|-------------|----------------|
| KS naïf | $D_{KS}$ vs $\mathrm{Exp}(1)$ | $D=0.031$, $p=0.12$ |
| KS corrigé autocorr. | Correction Bartlett $n_\mathrm{eff}$ | $p_\mathrm{adj}=0.04$ |
| $\chi^2$ conditionnel | 15 classes, tranche $\log_{10}m \in [4, 4.5]$ | $\chi^2=13.1$, $p=0.52$ |

---

### 4. `permutation_tests.py` — Tests de permutation T1/T2/T3

Trois méthodes pour établir que les corrélations spectrales de $Q_r$ avec les
zéros de Riemann sont des **artefacts de régression spurieuse** [Granger & Newbold 1974].

```bash
python permutation_tests.py --quick     # N=500, B=200
python permutation_tests.py             # N=10000, B=5000 (~9 min)
```

**Résultats reproduits (Table 5 du papier) :**

| Méthode | $z(\gamma_1)$ | $p(\gamma_1)$ | Signif. obs. | Signif. perm. |
|---------|-------------|-------------|-------------|--------------|
| T1 Mélange | $-0.23\sigma$ | 0.550 | 8/15 | 13.5/15 |
| T2 Blocs (50) | $-0.21\sigma$ | 0.547 | 8/15 | 13.5/15 |
| T3 Theiler | $-0.02\sigma$ | 0.514 | 8/15 | 13.5/15 |

Scores $z$ négatifs → signal observé **plus faible** qu'une permutation aléatoire.

---

### 5. `Q_decomposition.py` — Décomposition $Q \approx \alpha\sqrt{p} + \beta D_N$

Calcule la décomposition structurelle sous trois conventions paramétriques
indépendantes, la décomposition de variance et l'annulation arithmétique.

```bash
python Q_decomposition.py --quick     # N=2000
python Q_decomposition.py             # N=10000 (papier)
```

**Résultats reproduits (Tables 6–8 du papier) :**

| Convention | $R^2$ | $r(Q_r, D_N)$ | Annulation | Hurst $H$ |
|-----------|------|--------------|-----------|----------|
| Canonique | 0.967 | +0.869 | 99.9964 % | 0.392 |
| Min \|q\| | 0.979 | +0.846 | 99.9899 % | 0.398 |
| Min \|c_i\| | 0.980 | +0.843 | 99.9889 % | 0.395 |

Décomposition de variance (convention canonique) :
- $\alpha\sqrt{p}$ → **51.9 %** de Var$(Q)$
- $\beta D_N$ → **44.8 %** de Var$(Q)$
- Résidu pur → **3.3 %** de Var$(Q)$

---

## Installation

```bash
git clone https://github.com/BAKKAOUI/parametric-primes-2026.git
cd parametric-primes-2026
pip install -r requirements.txt
```

## Exécution rapide (test du pipeline)

```bash
python bh_constants.py                          # ~1 s
python qmin_compute.py --quick --no-plot        # ~0.2 s
python distrib_validation.py --quick --no-plot  # ~2 s
python permutation_tests.py --quick --no-plot   # ~0.3 s
python Q_decomposition.py --quick --no-plot     # ~0.1 s
```

## Reproductibilité complète (papier)

```bash
# Étape 1 : constantes BH (P_MAX=200000, ~30 min)
python bh_constants.py

# Étape 2 : calcul q_min (M=10^6, ~5h Python pur)
python qmin_compute.py --M 1000000

# Étape 3 : validation distributionnelle
python distrib_validation.py --input qmin_k21_M1000000.npy

# Étape 4 : tests de permutation (N=10000, B=5000, ~9 min)
python permutation_tests.py --N 10000 --B 5000

# Étape 5 : décomposition Q (N=10000, ~0.5 s)
python Q_decomposition.py --N 10000
```

---

## Dépendances

| Package | Version testée | Usage |
|---------|---------------|-------|
| Python | ≥ 3.10 | — |
| NumPy | ≥ 1.24 | Calcul vectoriel |
| SciPy | ≥ 1.10 | Tests statistiques (KS, χ²) |
| SymPy | ≥ 1.12 | Test de primalité BPSW, crible |
| Matplotlib | ≥ 3.7 | Figures (optionnel) |

---

## Licence

MIT License — voir [LICENSE](LICENSE).

---

## Citation

```bibtex
@misc{bakkaoui2026parametric,
  author       = {Bakkaoui, Hassane},
  title        = {Une famille paramétrique de nombres premiers via les formes
                  quadratiques},
  year         = {2026},
  note         = {Prépublication, arXiv:XXXX.XXXXX [math.NT]},
  url          = {https://github.com/BAKKAOUI/parametric-primes-2026}
}
```

---

## Contact

Hassane BAKKAOUI — Chercheur indépendant  
✉ bakkahassa@hotmail.com
