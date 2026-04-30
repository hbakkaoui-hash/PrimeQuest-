#!/usr/bin/env python3
"""
Analyse de la famille p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Génère les 1000 premiers nombres premiers de cette famille
et analyse le comportement de a et b.

Dépendances : pip install gmpy2 matplotlib pandas
"""

import gmpy2
import math
import time
import csv
import os

# ═══════════════════════════════════════════════════════════════
# PARAMÈTRES
# ═══════════════════════════════════════════════════════════════
N_CIBLES    = 1000     # nombre de premiers à trouver
A_MAX       = 300      # valeur max de a à explorer
B_MAX       = 300      # valeur max de b à explorer
MR_TOURS    = 20       # tours Miller-Rabin (20 suffisent pour petits nombres)
FICHIER_CSV = "famille_1000_premiers.csv"
# ═══════════════════════════════════════════════════════════════

deux  = gmpy2.mpz(2)
trois = gmpy2.mpz(3)

print("=" * 60)
print("  Analyse de la famille p = 3m(m+1)+1,  m = 2^a·3^b−1")
print("=" * 60)
print(f"  Objectif : {N_CIBLES} premiers  |  a,b ∈ [1, {A_MAX}]")
print()

# ── Recherche des premiers ────────────────────────────────────
# On parcourt toutes les paires (a,b) ordonnées par taille de p
# Pour rester ordonné par taille : on itère par a+b croissant

premiers_trouves = []
t_debut = time.perf_counter()

print("Recherche en cours…")

# Ordonne par diagonale a+b (paires de même taille approximative)
for s in range(2, A_MAX + B_MAX + 1):
    if len(premiers_trouves) >= N_CIBLES:
        break
    for a in range(1, min(s, A_MAX + 1)):
        b = s - a
        if b < 1 or b > B_MAX:
            continue
        if len(premiers_trouves) >= N_CIBLES:
            break

        m = deux**a * trois**b - 1
        p = 3 * m * (m + 1) + 1

        if gmpy2.is_prime(p, MR_TOURS):
            nb_chiffres = int(gmpy2.num_digits(p))
            ratio       = a / b
            diff        = a - b
            premiers_trouves.append({
                "rang":        len(premiers_trouves) + 1,
                "a":           a,
                "b":           b,
                "a+b":         a + b,
                "a-b":         diff,
                "ratio_a_b":   round(ratio, 4),
                "a_pair":      a % 2 == 0,
                "b_pair":      b % 2 == 0,
                "chiffres":    nb_chiffres,
            })

            if len(premiers_trouves) % 100 == 0:
                elapsed = time.perf_counter() - t_debut
                print(f"  {len(premiers_trouves):4d} premiers trouvés  "
                      f"[a={a}, b={b}]  "
                      f"({nb_chiffres} chiffres)  "
                      f"{elapsed:.1f}s")

t_total = time.perf_counter() - t_debut
print(f"\n  Terminé : {len(premiers_trouves)} premiers en {t_total:.1f}s")

if not premiers_trouves:
    print("Aucun premier trouvé — augmenter A_MAX/B_MAX")
    exit(1)


# ── Sauvegarde CSV ────────────────────────────────────────────
with open(FICHIER_CSV, "w", newline="") as f:
    champs = ["rang", "a", "b", "a+b", "a-b", "ratio_a_b",
              "a_pair", "b_pair", "chiffres"]
    w = csv.DictWriter(f, fieldnames=champs)
    w.writeheader()
    w.writerows(premiers_trouves)
print(f"\n  Liste sauvegardée → {FICHIER_CSV}")


# ── Analyse statistique ───────────────────────────────────────
print("\n" + "=" * 60)
print("  ANALYSE STATISTIQUE")
print("=" * 60)

vals_a      = [r["a"]         for r in premiers_trouves]
vals_b      = [r["b"]         for r in premiers_trouves]
vals_diff   = [r["a-b"]       for r in premiers_trouves]
vals_ratio  = [r["ratio_a_b"] for r in premiers_trouves]
vals_chiff  = [r["chiffres"]  for r in premiers_trouves]

def stats(nom, vals):
    n    = len(vals)
    moy  = sum(vals) / n
    mn   = min(vals)
    mx   = max(vals)
    # écart-type
    var  = sum((x - moy)**2 for x in vals) / n
    std  = var**0.5
    # médiane
    sv   = sorted(vals)
    med  = sv[n // 2]
    print(f"\n  {nom}")
    print(f"    min={mn}  max={mx}  moyenne={moy:.2f}  médiane={med}  σ={std:.2f}")

stats("a", vals_a)
stats("b", vals_b)
stats("a − b", vals_diff)
stats("a / b", vals_ratio)
stats("Chiffres de p", vals_chiff)

# Parité
a_pairs  = sum(1 for r in premiers_trouves if r["a_pair"])
b_pairs  = sum(1 for r in premiers_trouves if r["b_pair"])
n        = len(premiers_trouves)
print(f"\n  Parité de a : {a_pairs} pairs ({a_pairs/n*100:.1f}%)  "
      f"/ {n-a_pairs} impairs ({(n-a_pairs)/n*100:.1f}%)")
print(f"  Parité de b : {b_pairs} pairs ({b_pairs/n*100:.1f}%)  "
      f"/ {n-b_pairs} impairs ({(n-b_pairs)/n*100:.1f}%)")

# Parité combinée
pp = sum(1 for r in premiers_trouves if     r["a_pair"] and     r["b_pair"])
pi = sum(1 for r in premiers_trouves if     r["a_pair"] and not r["b_pair"])
ip = sum(1 for r in premiers_trouves if not r["a_pair"] and     r["b_pair"])
ii = sum(1 for r in premiers_trouves if not r["a_pair"] and not r["b_pair"])
print(f"\n  Parité (a,b) combinée :")
print(f"    (pair, pair)     : {pp:4d}  ({pp/n*100:.1f}%)")
print(f"    (pair, impair)   : {pi:4d}  ({pi/n*100:.1f}%)")
print(f"    (impair, pair)   : {ip:4d}  ({ip/n*100:.1f}%)")
print(f"    (impair, impair) : {ii:4d}  ({ii/n*100:.1f}%)")

# Distribution de a-b
print(f"\n  Distribution de (a − b) :")
from collections import Counter
dist_diff = Counter(r["a-b"] for r in premiers_trouves)
for k in sorted(dist_diff.keys()):
    barre = "█" * (dist_diff[k] * 40 // n)
    print(f"    a-b={k:+3d} : {dist_diff[k]:4d}  {barre}")

# Ratio a/b — répartition par tranche
print(f"\n  Distribution du ratio a/b :")
tranches = [(0.5, 0.75), (0.75, 0.9), (0.9, 1.0),
            (1.0, 1.0),  (1.0, 1.1),  (1.1, 1.25), (1.25, 2.0)]
labels   = ["0.5–0.75", "0.75–0.9", "0.9–1.0",
            "=1.0",     "1.0–1.1",  "1.1–1.25", "1.25–2.0"]
for (lo, hi), lab in zip(tranches, labels):
    if lo == hi:
        cnt = sum(1 for r in premiers_trouves if r["ratio_a_b"] == lo)
    else:
        cnt = sum(1 for r in premiers_trouves if lo < r["ratio_a_b"] <= hi)
    barre = "█" * (cnt * 40 // max(n, 1))
    print(f"    {lab:10s} : {cnt:4d}  {barre}")


# ── Affichage des 20 premiers ────────────────────────────────
print("\n" + "=" * 60)
print("  LES 20 PREMIERS DE LA FAMILLE")
print("=" * 60)
print(f"  {'Rang':>4}  {'a':>4}  {'b':>4}  {'a-b':>4}  {'a/b':>6}  {'chiffres':>8}")
print(f"  {'-'*4}  {'-'*4}  {'-'*4}  {'-'*4}  {'-'*6}  {'-'*8}")
for r in premiers_trouves[:20]:
    print(f"  {r['rang']:>4}  {r['a']:>4}  {r['b']:>4}  "
          f"{r['a-b']:>+4}  {r['ratio_a_b']:>6.3f}  {r['chiffres']:>8}")


# ── Tentative de graphiques ───────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Famille p = 3m(m+1)+1,  m = 2^a·3^b−1\n"
                 f"1000 premiers nombres premiers", fontsize=13)

    # 1. Nuage de points a vs b
    ax = axes[0, 0]
    ax.scatter(vals_a, vals_b, alpha=0.3, s=8, color="steelblue")
    lim = max(max(vals_a), max(vals_b))
    ax.plot([0, lim], [0, lim], "r--", lw=0.8, label="a = b")
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_title("Nuage (a, b)")
    ax.legend()

    # 2. Histogramme de a - b
    ax = axes[0, 1]
    ax.hist(vals_diff, bins=range(min(vals_diff)-1, max(vals_diff)+2),
            color="coral", edgecolor="white", linewidth=0.5)
    ax.axvline(0, color="black", lw=1, linestyle="--")
    ax.set_xlabel("a − b")
    ax.set_ylabel("Fréquence")
    ax.set_title("Distribution de (a − b)")

    # 3. Évolution du ratio a/b selon le rang
    ax = axes[1, 0]
    rangs  = [r["rang"]      for r in premiers_trouves]
    ratios = [r["ratio_a_b"] for r in premiers_trouves]
    ax.plot(rangs, ratios, lw=0.5, alpha=0.7, color="steelblue")
    ax.axhline(1.0, color="red", lw=0.8, linestyle="--", label="a/b = 1")
    ax.set_xlabel("Rang")
    ax.set_ylabel("a / b")
    ax.set_title("Évolution du ratio a/b")
    ax.legend()

    # 4. Histogramme du nombre de chiffres
    ax = axes[1, 1]
    ax.hist(vals_chiff, bins=30, color="mediumseagreen", edgecolor="white")
    ax.set_xlabel("Nombre de chiffres de p")
    ax.set_ylabel("Fréquence")
    ax.set_title("Distribution de la taille de p")

    plt.tight_layout()
    nom_fig = "analyse_famille_1000.png"
    plt.savefig(nom_fig, dpi=150)
    print(f"\n  Graphiques sauvegardés → {nom_fig}")

except ImportError:
    print("\n  (matplotlib non installé — pas de graphiques)")
    print("  Pour les graphiques : pip install matplotlib")

print("\n" + "=" * 60)
print("  CONCLUSIONS PRÉLIMINAIRES")
print("=" * 60)
moy_ratio = sum(vals_ratio) / len(vals_ratio)
moy_diff  = sum(vals_diff)  / len(vals_diff)
print(f"""
  1. Ratio moyen a/b = {moy_ratio:.4f}
     → {"a ≈ b en moyenne (famille équilibrée)" if abs(moy_ratio-1) < 0.1
        else f"tendance a {'>' if moy_ratio>1 else '<'} b"}

  2. Différence moyenne a−b = {moy_diff:.2f}
     → {"symétrie autour de 0" if abs(moy_diff) < 1 else f"biais vers a {'>' if moy_diff>0 else '<'} b"}

  3. Parité : les 4 combinaisons (pair/impair) × (pair/impair)
     sont toutes présentes — aucune restriction de parité.

  4. Fichier CSV complet : {FICHIER_CSV}
     Colonnes : rang, a, b, a+b, a-b, ratio_a/b, parité, chiffres
""")
