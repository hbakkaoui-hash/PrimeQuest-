# Synthèse de session — PrimeQuest
**Date :** 21 avril 2026  
**Repo :** `hbakkaoui-hash/PrimeQuest-`  
**Branche de travail :** `claude/optimize-prime-finder-JrSNV`  
**À coller en début de nouvelle conversation Claude pour reprendre le contexte.**

---

## 1. Contexte du projet

**Auteur :** Hassane BAKKAOUI

Le dépôt `PrimeQuest-` est un mono-repo mathématique en cours de structuration. Il contient (ou contiendra) deux papiers :

### Paper 1 — Famille paramétrique / constantes BH / distribution de q_min
- Titre complet : *Une famille paramétrique de nombres premiers via les formes quadratiques — Distribution heuristique du décalage minimal, validation numérique à grande échelle et contrôle statistique des corrélations spectrales spurieuses.*
- Famille étudiée : `p_{k,m,ε,q} = k·m(m+1) + ε + 2kq`
- Loi heuristique : `E[|q_min(k,m)|] ~ ln(m) / C_k`
- Fichiers existants : `src/main.tex`, `src/*.py` (5 scripts), `figures/` (21 fichiers), `data/TABLE_full_subsample_10e8.csv`, `requirements.txt`, `Parametric Family of Primes. pdf.pdf`

### Paper 2 — Grands nombres premiers prouvés via Pocklington
- Famille : `p = 3m(m+1)+1` avec `m = 2^a · 3^b − 1`
- Preuve : théorème de Pocklington N-1
- **Statut : work in progress** — script opérationnel, papier à écrire
- Fichiers : `primequest.py`, `premier_10000chiffres.txt`

---

## 2. Structure cible du repo (validée, PAS encore exécutée)

La réorganisation en mono-repo multi-papiers a été **proposée et validée** par l'auteur mais **les `git mv` n'ont pas encore été faits**. Plan exact :

```
PrimeQuest-/
├── README.md                              ← nouveau (index des deux papiers)
├── .gitignore                             ← mis à jour
├── paper1-qmin-bh/
│   ├── README.md                          ← à créer
│   ├── Parametric Family of Primes. pdf.pdf  ← reste à la racine de paper1
│   ├── tex/
│   │   └── main.tex                       ← git mv src/main.tex
│   ├── figures/                           ← git mv figures/
│   ├── code/
│   │   ├── Q_decomposition.py             ← git mv src/Q_decomposition.py
│   │   ├── bh_constants.py                ← git mv src/bh_constants.py
│   │   ├── distrib_validation.py          ← git mv src/distrib_validation.py
│   │   ├── permutation_tests.py           ← git mv src/permutation_tests.py
│   │   └── qmin_compute.py                ← git mv src/qmin_compute.py
│   ├── data/
│   │   └── TABLE_full_subsample_10e8.csv  ← git mv data/
│   └── requirements.txt                   ← git mv
└── paper2-pocklington-primes/
    ├── README.md                          ← à créer
    ├── tex/                               ← squelette vide
    ├── figures/                           ← squelette vide
    ├── code/
    │   └── primequest.py                  ← git mv primequest.py
    └── data/
        ├── premier_10000chiffres.txt      ← git mv
        ├── premier_200chiffres.txt        ← git mv
        └── premier_50chiffres.txt         ← git mv
```

**Prochaine action : valider et exécuter les `git mv`.**

---

## 3. Le script `primequest.py`

### Description
Script Python de recherche et **preuve de primalité** pour n'importe quelle taille,  
fondé sur le **théorème de Pocklington N-1**.

### Trois méthodes (ordre d'essai)

**Méthode 1 — Famille paramétrique** (la plus rapide quand elle s'applique)
- `p = 3m(m+1)+1`, `m = 2^a · 3^b − 1`
- `p−1 = 2^a · 3^(b+1) · m` → factorisation complète connue
- `F = 2^a · 3^(b+1) > √p` → Pocklington direct
- Témoins cherchés pour q=2 et q=3 uniquement

**Méthode 2 — Nombres de Proth** (fallback)
- `p = k · 2^n + 1`, k impair, k < 2^n
- Preuve : `∃a : a^((p-1)/2) ≡ −1 (mod p)` (théorème de Proth, 1878)

**Méthode 3 — Pocklington récursif (Maurer)** (garantie universelle)
- Certifie un sous-premier `q` de ≈ `bits/2` bits (récursivement)
- Cherche `p = 2·k·q + 1`, F = 2q > √p
- Témoins pour q=2 et q (grand premier certifié)
- Profondeur ≈ 7 niveaux pour 10 000 chiffres

### Paramètre principal
```python
DIGITS_CIBLE = 1000  # modifier ici ou : python3 primequest.py 10000
```

### Dépendance
```
gmpy2  (binding Python pour GNU MP — arithmétique grande précision)
```

---

## 4. Résultat obtenu — Premier de 10 000 chiffres PROUVÉ

```
a = 10094
b = 4110
m = 2^10094 · 3^4110 − 1       (5 000 chiffres)
p = 3m(m+1) + 1                (10 000 chiffres)

Preuve Pocklington :
  F = 2^10094 · 3^4111          (5 001 chiffres) > √p  ✅
  q = 2  →  témoin w = 11       ✅
  q = 3  →  témoin w =  5       ✅

Temps de calcul : 2h 05min
  dont énumération des paires : 20min 35s  ← goulot identifié
  dont tests MR + Pocklington  : ~35min
```

Le nombre complet est dans `premier_10000chiffres.txt` (commité sur GitHub).

---

## 5. Performances et limites connues

| Taille | Temps MR | Temps total estimé | Statut |
|---|---|---|---|
| 1 000 chiffres | < 1ms | ~1s | ✅ immédiat |
| 10 000 chiffres | 0.44ms | ~2h | ✅ fait |
| 30 000 chiffres | 51s | ~100h | ⚠️ trop lent |

### Goulot identifié : l'énumération
Pour 10 000 chiffres, `_enumerer_paires()` génère **188 595 paires (a,b)** toutes en mémoire avant de tester. Cela prend 20 minutes car chaque paire implique le calcul de `2^a · 3^b` sur des milliers de chiffres.

**Optimisation à faire (non encore implémentée) :**  
Générer et tester les paires **à la volée** (streaming), sans tout stocker en mémoire.  
Cela réduirait l'empreinte mémoire et le temps d'initialisation de ~20min à quelques secondes.

---

## 6. Position dans la littérature mathématique

| Élément | Originalité |
|---|---|
| Théorème de Pocklington | Classique (1914) |
| Pocklington récursif (Maurer) | Publié (1995), dans PARI/GP |
| Proth | Classique (1878) |
| Famille `p = 3m(m+1)+1, m = 2^a·3^b−1` | **Original** — construction propre à ce projet |
| Combinaison famille + recherche systématique | **Original** — pas de script équivalent connu |

La valeur scientifique réside dans la **famille paramétrique** : elle produit des premiers dont `p−1` est complètement factorisable par construction, ce qui est rare et utile pour la cryptographie et la théorie des nombres.

---

## 7. Ce qui reste à faire (prochaines sessions)

### Priorité 1 — Réorganisation du repo
- Exécuter les `git mv` selon le plan de la section 2
- Créer les README de paper1 et paper2
- Mettre à jour le README racine

### Priorité 2 — Optimisation du script
- Réécrire `_enumerer_paires()` en mode streaming (tester à la volée)
- Objectif : réduire l'énumération de 20min à <1min pour 10 000 chiffres
- Cela rendra aussi les recherches à 30 000+ chiffres plus réalistes

### Priorité 3 — Paper 2
- L'auteur a un papier en cours (PDF + LaTeX) à intégrer dans `paper2-pocklington-primes/tex/`
- Il souhaite le partager via drag & drop dans Claude Code
- Le script `primequest.py` et le premier de 10 000 chiffres seront les résultats numériques du papier

### Priorité 4 — Recherche de premiers plus grands
- Objectif déclaré : 30 000 chiffres (nécessite l'optimisation de l'énumération)
- Contrainte : travail par tranches de ~5% du forfait hebdomadaire (≈ 5 échanges actifs)

---

## 8. Règles de travail convenues

- **Travail par tranches** : chaque tranche ≤ 5% du forfait hebdomadaire, puis arrêt et attente de confirmation
- **Push automatique** : Claude pousse sur `claude/optimize-prime-finder-JrSNV` — ne jamais push sur `main` sans accord explicite
- **Pas de PR automatique** : créer une PR uniquement si demandé explicitement
- **Commit immédiat** : tout fichier généré (premiers trouvés, scripts) est commité et poussé

---

*Fin de synthèse — coller ce document au début d'une nouvelle conversation pour reprendre sans perte de contexte.*
