# PrimeQuest — Recherche en théorie des nombres premiers

**Auteur :** Hassane Bakkaoui — Chercheur indépendant  
**Période :** Avril – Mai 2026

Dépôt de recherche contenant quatre papiers sur les nombres premiers paramétriques,
la distribution de $q_{\min}$, et les certificats de primalité.

---

## Structure du dépôt

```
PrimeQuest-/
├── paper1-qmin-bh/                  Papier I  — Famille paramétrique & constante BH
├── paper2-prime-families/           Papier II — Loi conditionnelle & connexion ζ(s)
├── paper3-padic-qmin/               Papier III — Distribution p-adique de q_min
└── paper4-pocklington-certificates/ Papier IV — Certificats de primalité Pocklington
```

---

## Papier I — Famille paramétrique de premiers & analyse Bateman–Horn

**Dossier :** [`paper1-qmin-bh/`](paper1-qmin-bh/)  
**MSC :** 11N13 · 11Y11 · 11N32

Étude de la loi marginale du décalage minimal $q_{\min}$ dans la famille paramétrée
$p = km(m+1) + \varepsilon + 2kq$, et vérification numérique de la conjecture
de Bateman–Horn. Résultat principal : $\mathbb{E}[|q_{\min}|] \sim \log m / C_k$.

---

## Papier II — Loi conditionnelle & exclusion de la connexion spectrale

**Dossier :** [`paper2-prime-families/`](paper2-prime-families/)  
**MSC :** 11N13 · 11M41 · 11Y11

Loi conditionnelle $\mathbb{E}[|q| \mid m] = m/4$ et réfutation de la connexion
spurieuse à $\zeta(s)$. Analyse de la couche spectrale de la famille.

---

## Papier III — Distribution $p$-adique de $q_{\min}$ : phénomène d'anti-persistance

**Dossier :** [`paper3-padic-qmin/`](paper3-padic-qmin/)  
**MSC :** 11N32 · 11N13 · 11M41 · 62P99

Distribution jointe modulaire de $q_{\min} \bmod \ell$ pour les petits premiers $\ell$.
Théorèmes conditionnel et inconditionnel sur la fréquence de permission $\Pi_P$.
Découverte d'un **phénomène d'anti-persistance** : inversion de signe de l'autocorrélation
lag-1 à un seuil $k_0 \approx 8.4$, avec loi linéaire $\mathrm{ACF}_1 \approx 0.024k - 0.20$
($R^2 = 0.97$).

---

## Papier IV — Certificats de primalité inconditionnels : $p = 3m(m+1)+1$

**Dossier :** [`paper4-pocklington-certificates/`](paper4-pocklington-certificates/)  
**MSC :** 11A41 · 11Y11 · 11N13

Certification inconditionnelle (Théorème de Pocklington N−1) de la sous-famille
$m = 2^a \cdot 3^b - 1$. Trois filtres arithmétiques éliminent **87,1 %** des candidats.
Algorithme multi-cœurs PrimeQuest v1–v7. Record actuel : premier de **29 998 chiffres**
(3 h 06 min, 9 workers, ARM64).

| Chiffres | a | b | Temps | Workers |
|---------:|--:|--:|------:|--------:|
| 9 998 | 6 212 | 6 738 | — | 1 |
| 10 000 | 6 213 | 6 740 | — | 1 |
| 19 999 | 12 228 | 13 242 | 2 h 40 min | 9 |
| **29 998** | **19 435** | **19 173** | **3 h 06 min** | **9** |

---

## Installation

```bash
pip install gmpy2 matplotlib pandas
```

Python ≥ 3.9, gmpy2 ≥ 2.1 (backend GMP).  
ARM64 (Snapdragon / Apple Silicon) : utiliser Python natif ARM64 pour performances optimales.
