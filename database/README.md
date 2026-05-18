# Base de données unifiée — Famille paramétrique de premiers

**Formule centrale :**  `p = k·m(m+1) + ε + 2k·q`

| Champ       | Valeur                          |
|-------------|----------------------------------|
| k           | 3, 15, 21, 23                   |
| ε           | +1 ou −1                        |
| m           | entier ≥ 0 (indice de couche)   |
| q           | entier signé (décalage minimal) |

Cette base est le socle commun des Papers I–IV et du programme de recherche
sur la connexion GUE / zéros de Riemann (H. Bakkaoui, 2026).

---

## Structure des fichiers

```
prime_database/
├── parquet/
│   ├── k=3/      part_000000.parquet, part_000001.parquet, …
│   ├── k=15/
│   ├── k=21/
│   └── k=23/
├── csv/          k3_000000.csv.gz, k3_000001.csv.gz, …
├── raw/          raw_k3.parquet, raw_k15.parquet, …  (intermédiaire)
├── primes_parametric.sqlite
├── checkpoint.json
└── generation.log
```

---

## Schéma exact — table `primes_parametric`

Une ligne par `(premier p, paramètre k)`.  
Tous les entiers en **int64**, tous les flottants en **float64**.

### Groupe A — Paramètres fondamentaux

| Colonne  | Type   | Définition |
|----------|--------|------------|
| `p`      | int64  | Le premier |
| `k`      | int64  | Paramètre de famille ∈ {3, 15, 21, 23} |
| `m`      | int64  | Couche canonique (voir convention ci-dessous) |
| `eps`    | int64  | ε ∈ {+1, −1}, déterminé par p mod 2k |
| `q`      | int64  | Décalage signé = q\_min (voir convention) |
| `q_min`  | int64  | Identique à q (= q par construction canonique) |
| `H_m`    | int64  | Base hexagonale généralisée = k·m·(m+1) |

### Groupe B — Structure de couche

| Colonne  | Type   | Définition |
|----------|--------|------------|
| `r`      | int64  | Rang dans la k-famille (ordre croissant de p), 1-indexé |
| `n_m`    | int64  | Nombre de premiers dans la couche m pour ce k |
| `H_m`    | int64  | (répété dans le Groupe A) |

### Groupe C — Fonctionnelles cumulatives

Calculées dans l'ordre croissant de p au sein de chaque k-famille.

| Colonne  | Type    | Définition |
|----------|---------|------------|
| `Q_r`    | float64 | Σᵢ₌₁ʳ qᵢ  (somme cumulative des décalages signés) |
| `X_r`    | float64 | Σᵢ₌₁ʳ (qᵢ + εᵢ)  (= Σcᵢ ; terme Q du papier) |
| `D_r`    | float64 | Σᵢ₌₁ʳ ({√(pᵢ'/k)} − ½) avec pᵢ' = pᵢ − εᵢ  (terme de Weyl / équidistribution) |

*Notation : {x} = partie fractionnaire de x.*

### Groupe D — Distribution modulaire

| Colonne      | Type  | Définition |
|--------------|-------|------------|
| `q_mod_2`   | int64 | q mod 2    |
| `q_mod_3`   | int64 | q mod 3    |
| `q_mod_5`   | int64 | q mod 5    |
| `q_mod_7`   | int64 | q mod 7    |
| `q_mod_11`  | int64 | q mod 11   |
| `q_mod_13`  | int64 | q mod 13   |

### Groupe E — Grandeurs analytiques par couche

| Colonne   | Type    | Définition |
|-----------|---------|------------|
| `ln_p`    | float64 | log(p) (logarithme naturel) |
| `sqrt_p`  | float64 | √p |
| `ln_m`    | float64 | log(m) si m > 0, sinon 0 |
| `nm_hat`  | float64 | Espérance théorique n̂_m (voir formule ci-dessous) |
| `T_m`     | float64 | (n_m − n̂_m) / √n̂_m  (statistique de Poisson normalisée) |

---

## Convention de signe pour q et m (q\_min)

Pour un premier p et un paramètre k :

1. **ε** est déterminé uniquement par `p mod 2k` :
   - `p ≡ 1 (mod 2k)` → ε = +1
   - `p ≡ 2k−1 (mod 2k)` → ε = −1
   - Sinon : p n'appartient pas à la k-famille.

2. Poser `N = (p − ε) / k` (entier exact, toujours pair).

3. Calculer `m_low = ⌊(−1 + √(1+4N)) / 2⌋`  (le plus grand m avec m(m+1) ≤ N).

4. Calculer :
   - `q_low = (N − m_low·(m_low+1)) / 2`  ≥ 0
   - `q_high = q_low − (m_low+1)`  < 0

5. **Choix canonique** :
   - Si `|q_high| < |q_low|` → m = m_low+1, q = q_high  (q négatif)
   - Sinon (|q_low| ≤ |q_high|) → m = m_low, q = q_low  (q positif, priorité en cas d'égalité)

**q_min = q** (identiques par construction).  
La borne inconditionnelle garantit `|q| ≤ m` (Proposition 2.3 du papier).

---

## Formule de l'espérance théorique n̂_m

```
n̂_m = 4k(m+1) / (φ(2k) · log(k·m·(m+1)))
```

Valeurs de φ(2k) :

| k  | 2k | φ(2k) |
|----|-----|-------|
| 3  | 6   | 2     |
| 15 | 30  | 8     |
| 21 | 42  | 12    |
| 23 | 46  | 22    |

---

## Tables SQLite

### `primes_parametric_sample`
Échantillon de 50 000 lignes par k (25 000 premières + 25 000 dernières)
avec toutes les colonnes ci-dessus. Pour les requêtes ad hoc.
Index : `(p, k)`.

### `layers`
Une ligne par `(k, m)` avec :
- `H_m`, `n_m`, `nm_hat`, `T_m`
- `sum_q_r` : Q_r à la dernière ligne de la couche
- `avg_ln_p` : moyenne de log(p) dans la couche
Index : `(k, m)`.

### `modular_stats`
Distribution de q mod l par blocs de 100 000 premiers, avec comptages
`cnt_modL_v` pour chaque (L ∈ {2,3,5,7,11,13}, v ∈ [0, L−1]).

---

## Utilisation rapide

```python
from query_tools import query_layer, compute_c

# Tous les premiers de la couche m=100 pour k=3
df = query_layer(m=100, k=3)
print(df[["p", "eps", "q", "n_m", "T_m"]])

# Statistique chi-deux sur les couches 10 à 1000 pour k=15
c_val, detail = compute_c(m_min=10, m_max=1000, k=15)
print(f"c = {c_val:.4f}  (≈1 si distribution de Poisson)")
```

---

## Lancement

```bash
# Test rapide (N = 10^6)
python generate_database.py --quick

# Production standard (N = 10^9)
python generate_database.py --n 1e9 --workers 4

# Production haute échelle (N = 10^12, défaut)
# Estimation : ~30 jours mono-nœud ; prévoir cluster HPC
# Reprise automatique par checkpoint toutes les 10^8 premiers
python generate_database.py --workers $(nproc)

# Vérification indépendante (5 premiers + 5 derniers par k)
python verify_database.py

# Vérification globale complète (tous les fichiers Parquet)
python verify_database.py --global

# Reprise après interruption
python generate_database.py --resume
```

### Ressources nécessaires pour N = 10^12

| Ressource    | Estimation |
|--------------|-----------|
| π(10^12) ≈   | 37.6 milliards de premiers |
| Lignes k=3   | ~37.6 G (toutes les familles) |
| Stockage Parquet | ~2–4 To (Snappy) |
| Durée (mono-nœud, 8 cœurs) | ~5–7 jours |
| RAM peak (Phase 2) | ~8 Go par k |

---

## Effet de bord (couche terminale)

La couche `m_max` (la plus grande couche dans la base) est **tronquée** : elle
ne contient que les premiers jusqu'à N_MAX, pas jusqu'à H_{m_max+1}. Par conséquent :
- `n_m < nm_hat` pour la dernière couche → `T_m` fortement négatif
- Symptôme normal : `|T_m|` est maximal pour `m = m_max`

Pour l'analyse statistique, **exclure la dernière couche** par k :
```python
df_valid = df[df["m"] < df["m"].max()]
```

## Validation intégrée

À chaque chunk, le script vérifie :

1. **Reconstruction** : `p = H_m + ε + 2k·q` pour chaque ligne.
2. **Borne inconditionnelle** : `|q| ≤ m`.
3. **Annulation arithmétique** : `2k·|Σq| / Σp < 5 %`
   (la somme cumulative des décalages doit rester petite devant la somme des premiers,
   ce qui reflète l'approximation `p ≈ H_m` pour grands m).

Toute violation est loggée dans `prime_database/generation.log`.

---

## Connexion GUE / Zéros de Riemann

Le terme de Weyl `D_r = Σ({√(p'/k)} − ½)` teste l'équidistribution de
`√(p/k) mod 1` au sens de Weyl. La décomposition structurelle
`Q_r ≈ α·√p_r + β·D_r + γ` (Théorème 5.1) est le lien analytique entre
la famille paramétrique et le comportement des corrélations spectrales
à deux points (GUE). Voir le Paper III pour les détails.
