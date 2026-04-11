# =============================================================================
# Auteur  : Hassane BAKKAOUI — Chercheur indépendant
# Contact : bakkahassa@hotmail.com
# Dépôt   : https://github.com/BAKKAOUI/parametric-primes-2026
# Papier  : "Une famille paramétrique de nombres premiers via les formes
#            quadratiques — Distribution heuristique du décalage minimal,
#            validation numérique à grande échelle et contrôle statistique
#            des corrélations spectrales spurieuses"
# Auteur  : H. Bakkaoui (2026)
# arXiv   : math.NT [à compléter après soumission]
# Licence : MIT License — voir LICENSE dans le dépôt
# Date    : Avril 2026
# Version : 1.0.0
# =============================================================================
"""
bh_constants.py
===============
Calcul des constantes de Bateman–Horn par produit eulérien tronqué.

Référence :
    Bateman, P.T. & Horn, R.A. (1962). A heuristic asymptotic formula
    concerning the distribution of prime numbers. Math. Comp., 16(79), 363–367.

Famille étudiée :
    f_{k,ε}(m) = k·m·(m+1) + ε,   k ∈ {3, 15, 21},  ε ∈ {+1, -1}

Formule :
    C(f) = ∏_{p premier} (p - ω_f(p)) / (p - 1)

    ω_f(p) = #{m ∈ {0,...,p-1} : f(m) ≡ 0 (mod p)}
           = 1 + (k² - 4kε | p)   [symbole de Legendre]
      pour p impair, p ∤ k.
    ω_f(2) = 0  (car km(m+1) est pair, donc f(m) ≡ ε ≢ 0 mod 2).
    ω_f(p) = 0  pour p | k.

Résultats reproduits dans le papier (Table 2) :
    k=3  : C(f_{3,+}) = 3.361048,  C(f_{3,-}) = 2.783543,  C_3 = 6.144591
    k=15 : C(f_{15,+})= 3.034936,  C(f_{15,-})= 3.778642,  C_15= 6.813578
    k=21 : C(f_{21,+})= 3.876345,  C(f_{21,-})= 4.329956,  C_21= 8.206301

Usage :
    python bh_constants.py

Dépendances : sympy (génération des premiers), numpy
"""

import time
import math
import numpy as np
from sympy import primerange, jacobi_symbol

# ---------------------------------------------------------------------------
# Paramètres
# ---------------------------------------------------------------------------
K_VALUES   = [3, 15, 21]          # valeurs de k étudiées
EPS_VALUES = [+1, -1]             # ε ∈ {+1, -1}
P_MAX      = 200_000              # troncature du produit eulérien
# Erreur de troncature estimée < 0.1 % pour P_MAX = 200 000


# ---------------------------------------------------------------------------
# Noyau : calcul de ω_f(p) pour f_{k,ε}
# ---------------------------------------------------------------------------

def omega_f(p: int, k: int, eps: int) -> int:
    """
    Nombre de racines de f_{k,ε}(m) = k·m·(m+1) + ε dans Z/pZ.

    Règles :
      - p = 2            → ω = 0  (f(m) ≡ ε ≡ ±1 mod 2, jamais 0)
      - p | k            → ω = 0  (f(m) ≡ ε ≢ 0 mod p car ε = ±1)
      - p impair, p ∤ k  → ω = 1 + Jacobi(k² - 4kε, p)
    """
    if p == 2:
        return 0
    if k % p == 0:
        return 0
    disc = k * k - 4 * k * eps          # discriminant réduit
    leg  = int(jacobi_symbol(disc % p, p))   # symbole de Legendre/Jacobi
    return 1 + leg


# ---------------------------------------------------------------------------
# Produit eulérien tronqué
# ---------------------------------------------------------------------------

def bh_constant(k: int, eps: int, p_max: int = P_MAX) -> float:
    """
    Calcule C(f_{k,ε}) = ∏_{p ≤ p_max} (p - ω_f(p)) / (p - 1).

    Les termes pour p > p_max ont rapport (p - 1)/(p - 1) = 1 à l'ordre 1/p²,
    l'erreur de troncature est O(1/p_max) ≈ 5×10⁻⁶ pour p_max = 200 000.
    """
    log_C = 0.0
    for p in primerange(2, p_max + 1):
        omega = omega_f(p, k, eps)
        if omega >= p:          # sécurité (ne devrait pas arriver)
            continue
        # terme = (p - ω) / (p - 1)
        log_C += math.log(int(p) - omega) - math.log(int(p) - 1)
    return math.exp(log_C)


# ---------------------------------------------------------------------------
# Programme principal
# ---------------------------------------------------------------------------

def main():
    print("=" * 62)
    print("Constantes de Bateman–Horn — produit eulérien, p ≤ {:,}".format(P_MAX))
    print("=" * 62)
    print(f"{'k':>4}  {'ε':>3}  {'C(f_{k,ε})':>12}  {'1/C(f)':>10}  {'t (s)':>7}")
    print("-" * 62)

    results = {}
    for k in K_VALUES:
        C_total = 0.0
        for eps in EPS_VALUES:
            t0 = time.time()
            C  = bh_constant(k, eps, P_MAX)
            dt = time.time() - t0
            C_total += C
            sign = "+" if eps == +1 else "-"
            print(f"  k={k}  ε={sign}1  C = {C:.6f}   1/C = {1/C:.5f}   {dt:.1f}s")
            results[(k, eps)] = C
        print(f"       ↳  C_{{k={k}}} = C(f+) + C(f-) = {C_total:.6f}   "
              f"1/C_k = {1/C_total:.5f}")
        print()

    # Vérification par rapport aux valeurs du papier (Table 2)
    paper_values = {
        (3,  +1): 3.361048,  (3,  -1): 2.783543,
        (15, +1): 3.034936,  (15, -1): 3.778642,
        (21, +1): 3.876345,  (21, -1): 4.329956,
    }
    print("Vérification par rapport aux valeurs du papier (Table 2) :")
    print(f"{'(k,ε)':>8}  {'calculé':>10}  {'papier':>10}  {'écart relatif':>14}")
    print("-" * 50)
    for key, C_calc in results.items():
        C_ref  = paper_values[key]
        rel_err = abs(C_calc - C_ref) / C_ref
        k, eps = key
        sign   = "+" if eps == +1 else "-"
        status = "✓" if rel_err < 1e-4 else "✗"
        print(f"  ({k},{sign}1)  {C_calc:10.6f}  {C_ref:10.6f}  {rel_err:.2e}  {status}")

    print()
    print("Erreur de troncature estimée pour p_max = {:,} : < 0.1 %".format(P_MAX))


if __name__ == "__main__":
    main()
