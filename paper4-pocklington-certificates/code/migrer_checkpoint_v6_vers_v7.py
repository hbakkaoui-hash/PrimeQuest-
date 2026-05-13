#!/usr/bin/env python3
"""
Migration checkpoint v6 → v7
────────────────────────────
Lit  : checkpoint_50000_v6.json   (format PrimeQuest v6)
Écrit: checkpoint_50000_v7.json   (format PrimeQuest v7)

Les 7 ratios de v6 sont réinjectés aux bons indices dans la liste
TOUS_RATIOS de v7 (ordonnée centre→bords). Les 8 nouveaux ratios
de v7 démarrent à delta=0 (territoire vierge).

Usage : python migrer_checkpoint_v6_vers_v7.py
"""

import json, math, sys, os

# ── Paramètres v7 (doivent correspondre exactement à primequest_v7.py) ────────
DIGITS_CIBLE  = 50_000
RATIO_CENTRE  = 1.00
N_RATIOS_MAX  = 15
SPREAD_MAX    = 0.35
LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)

FILE_V6 = f"checkpoint_{DIGITS_CIBLE}_v6.json"
FILE_V7 = f"checkpoint_{DIGITS_CIBLE}_v7.json"

# ── Recalcul de TOUS_RATIOS (même logique que _generer_tous_ratios() dans v7) ─
def tous_ratios_v7():
    step  = (2 * SPREAD_MAX / (N_RATIOS_MAX - 1))
    bruts = [round(RATIO_CENTRE - SPREAD_MAX + i * step, 4) for i in range(N_RATIOS_MAX)]
    ci    = N_RATIOS_MAX // 2
    ordre = []
    for k in range(N_RATIOS_MAX):
        idx = ci - k // 2 if k % 2 == 0 else ci + (k + 1) // 2
        if 0 <= idx < N_RATIOS_MAX:
            ordre.append(bruts[idx])
    r_min = 2 * LOG10_2 / DIGITS_CIBLE
    r_max = (DIGITS_CIBLE / 2 - LOG10_2) / LOG10_3
    return [r for r in ordre
            if r_min <= r <= r_max
            and int(max(1, int(DIGITS_CIBLE / 2 / (LOG10_2 + r * LOG10_3))) * r) >= 1]

TOUS_RATIOS = tous_ratios_v7()

# ── Chargement checkpoint v6 ────────────────────────────────────────────
if not os.path.exists(FILE_V6):
    print(f"❌  Fichier {FILE_V6} introuvable.")
    sys.exit(1)

with open(FILE_V6, encoding="utf-8") as f:
    v6 = json.load(f)

if v6.get("version") != 6:
    print(f"⚠  Ce fichier n'est pas un checkpoint v6 (version={v6.get('version')}).")
    sys.exit(1)

v6_ratios = v6["ratios"]
v6_deltas = v6["deltas"]
v6_sides  = v6["sides"]

print(f"Checkpoint v6 chargé :")
print(f"  Ratios : {v6_ratios}")
print(f"  Deltas : {v6_deltas}")
print(f"  MR cumulés : {v6.get('n_mr', 0)}")
print(f"  Temps cumulé : {v6.get('t_cumul', 0)/3600:.2f}h")

# ── Mapping v6 → v7 ──────────────────────────────────────────────────────
v7_deltas = [0] * len(TOUS_RATIOS)
v7_sides  = [0] * len(TOUS_RATIOS)

transferes = []
non_trouves = []

for j, r in enumerate(v6_ratios):
    if r in TOUS_RATIOS:
        i = TOUS_RATIOS.index(r)
        v7_deltas[i] = v6_deltas[j]
        v7_sides[i]  = v6_sides[j]
        transferes.append((r, i, v6_deltas[j]))
    else:
        non_trouves.append(r)

print(f"\nMapping v6 → v7 :")
for r, i, d in transferes:
    statut = "épuisé" if d == -1 else f"delta={d}"
    print(f"  ratio {r:5.3f}  →  TOUS_RATIOS[{i:2d}]  ({statut})")

if non_trouves:
    print(f"\n⚠  Ratios v6 absents de TOUS_RATIOS v7 (perdus) : {non_trouves}")
    print(f"   Vérifier que RATIO_CENTRE, SPREAD_MAX et N_RATIOS_MAX sont identiques.")

nouveaux = [r for r in TOUS_RATIOS if r not in v6_ratios]
print(f"\nNouveaux ratios v7 (delta=0) : {nouveaux}")

# ── Écriture checkpoint v7 ────────────────────────────────────────────
v7 = {
    "digits":      DIGITS_CIBLE,
    "version":     7,
    "tous_ratios": TOUS_RATIOS,
    "deltas":      v7_deltas,
    "sides":       v7_sides,
    "n_t2":        v6.get("n_t2", 0),
    "n_crible":    v6.get("n_crible", 0),
    "n_mr":        v6.get("n_mr", 0),
    "n_paires":    v6.get("n_paires", 0),
    "t_cumul":     v6.get("t_cumul", 0.0),
}

with open(FILE_V7, "w", encoding="utf-8") as f:
    json.dump(v7, f, indent=2)

print(f"\n✅  Checkpoint v7 écrit → {FILE_V7}")
print(f"   TOUS_RATIOS v7 : {TOUS_RATIOS}")
print(f"   Deltas v7      : {v7_deltas}")
print(f"\nLancer ensuite : py primequest_v7.py")
