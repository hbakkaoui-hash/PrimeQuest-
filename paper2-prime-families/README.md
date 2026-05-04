# Papier II — Primes de Pocklington : p = 3m(m+1)+1

## Famille étudiée

```
p = 3·m·(m+1) + 1
m = 2^a · 3^b − 1
```

La primalité de p est prouvée par le **Théorème de Pocklington N-1** :
F = 2^a · 3^(b+1) divise p-1 et F > √p.

## Théorèmes prouvés

1. **p ≡ 1 (mod 6)** pour tout (a,b) — propriété structurelle
2. **Classes interdites mod 7** — 7|p ssi (a mod 3, b mod 6) ∈ F₇
3. **q|p impossible si q ≡ 2 (mod 3)** — réduit le crible de ~50%

## Résultats expérimentaux

| Chiffres | a | b | Témoins | Temps |
|---|---|---|---|---|
| 50 | — | — | — | < 1s |
| 200 | — | — | — | < 1s |
| ~10 000 | — | — | — | ~3 min |
| **19 999** | **12228** | **13242** | q=2→w=5, q=3→w=7 | **6h16** |

## Contenu

| Dossier / Fichier | Description |
|---|---|
| `code/primequest.py` | **PrimeQuest v3** — multi-cœurs, ARM64, checkpoint |
| `code/primequest_v1.py` | PrimeQuest v1 — version originale |
| `tex/paper_II_revised.tex` | Source LaTeX de l'article |
| `data/` | Certificats de primalité (valeurs de p et témoins) |
| `figures/` | Figures de l'article |

## Usage

```bash
# Lancer la recherche (20 000 chiffres, multi-cœurs)
python code/primequest.py

# Modifier les paramètres en haut du fichier :
#   DIGITS_CIBLE = 20_000
#   TIMEOUT_S    = 14_400
#   DELTA_DEPART = 0        # 0 = automatique
```

## Dépendance

```bash
pip install gmpy2
```

**ARM64** (Dell Snapdragon / Apple Silicon) : installer Python ARM64 natif
pour les performances maximales (~9× plus rapide avec 10 cœurs).
