#!/usr/bin/env python3
"""
Recherche optimisée — premier de 20 000 chiffres
Zigzag depuis a_centre, timeout 4 heures.
"""

import gmpy2
import math
import time
import sys

DIGITS_CIBLE = 20_000
TOLERANCE    = 10
MR_TOURS     = 25
TEMOINS_MAX  = 300
SIEVE_LIMIT  = 1_000_000
TIMEOUT_S    = 14_400     # 4 heures
RAPPORT_S    = 300        # rapport toutes les 5 min

deux  = gmpy2.mpz(2)
trois = gmpy2.mpz(3)
_LOG10_2 = math.log10(2)
_LOG10_3 = math.log10(3)

# ── Crible ───────────────────────────────────────────────────────────────────
def _eratosthene(lim):
    t = bytearray([1]) * (lim + 1)
    t[0] = t[1] = 0
    for i in range(2, int(lim**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    return [i for i in range(2, lim + 1) if t[i]]

print("Construction du crible…", end=" ", flush=True)
PREMIERS = _eratosthene(SIEVE_LIMIT)
print(f"{len(PREMIERS):,} premiers")

def crible_modulaire(a, b):
    for q in PREMIERS:
        if q <= 3:
            continue
        x     = (pow(2, a % (q-1), q) * pow(3, b % (q-1), q)) % q
        p_mod = (3 * (x-1) * x + 1) % q
        if p_mod == 0:
            return False
    return True

def b_valides(a, digits, tol):
    b_c = (digits/2 - a*_LOG10_2 - _LOG10_3/2) / _LOG10_3
    return [b for b in range(max(1, int(b_c)-3), int(b_c)+4)
            if abs(2*(a*_LOG10_2 + b*_LOG10_3) + _LOG10_3 - digits) <= tol+2]

def pocklington(p, a, b):
    p1 = p - 1
    v  = {2: None, 3: None}
    for w in range(2, 2 + TEMOINS_MAX):
        gw = gmpy2.mpz(w)
        if gmpy2.powmod(gw, p1, p) != 1:
            continue
        for q in [2, 3]:
            if v[q] is not None: continue
            if gmpy2.gcd(gmpy2.powmod(gw, p1//gmpy2.mpz(q), p) - 1, p) == 1:
                v[q] = w
        if all(x is not None for x in v.values()):
            return True, v
    return False, v

# ── Zigzag depuis le centre ───────────────────────────────────────────────────
a_centre = int((DIGITS_CIBLE / 2) / math.log10(6))
a_max    = int(DIGITS_CIBLE / (2 * _LOG10_2)) + 10

def zigzag(centre, a_min, a_max):
    yield centre
    delta = 1
    while True:
        up, down = centre + delta, centre - delta
        if up   <= a_max: yield up
        if down >= a_min: yield down
        if up > a_max and down < a_min:
            break
        delta += 1

print(f"\nRecherche de {DIGITS_CIBLE} chiffres")
print(f"a_centre={a_centre}  a_max={a_max}  timeout={TIMEOUT_S}s ({TIMEOUT_S/3600:.1f}h)")
print(f"Chaque MR ≈ 20s  —  référence 10 000 chiffres : 19 MR en 2min21s")
print(f"{'─'*60}\n")

n_crible = 0
n_mr     = 0
n_paires = 0
trouve   = None
t_debut  = time.perf_counter()
t_rapport= t_debut

for a in zigzag(a_centre, 1, a_max):

    elapsed = time.perf_counter() - t_debut

    # ── Timeout ──────────────────────────────────────────────────────────
    if elapsed >= TIMEOUT_S:
        print(f"\n  ⏱ TIMEOUT {TIMEOUT_S/3600:.1f}h atteint.")
        break

    for b in b_valides(a, DIGITS_CIBLE, TOLERANCE):
        n_paires += 1

        # ── Rapport toutes les 5 min ──────────────────────────────────────
        now = time.perf_counter()
        if now - t_rapport >= RAPPORT_S:
            elapsed = now - t_debut
            restant = TIMEOUT_S - elapsed
            print(f"  [{elapsed/60:5.1f}min]  a={a:7d}  "
                  f"paires={n_paires:8,}  crible={n_crible:7,}  "
                  f"MR={n_mr:4,}  restant={restant/60:.0f}min",
                  flush=True)
            t_rapport = now

        # ── Crible modulaire ──────────────────────────────────────────────
        if not crible_modulaire(a, b):
            n_crible += 1
            continue

        # ── Miller-Rabin ──────────────────────────────────────────────────
        n_mr += 1
        elapsed = time.perf_counter() - t_debut
        restant = TIMEOUT_S - elapsed
        print(f"  [{elapsed/60:5.1f}min]  MR #{n_mr:4d}  "
              f"a={a}, b={b}  (restant {restant/60:.0f}min) …",
              flush=True)

        m = deux**a * trois**b - 1
        p = 3 * m * (m+1) + 1
        if not gmpy2.is_prime(p, MR_TOURS):
            print(f"  → composé", flush=True)
            continue

        print(f"  *** PROBABLE — Pocklington… ***", flush=True)
        ok, temoins = pocklington(p, a, b)
        if ok:
            trouve = (a, b, m, p, int(gmpy2.num_digits(p)), temoins)
            break

    if trouve:
        break

t_total = time.perf_counter() - t_debut

# ── Compte rendu ─────────────────────────────────────────────────────────────
print(f"\n{'═'*65}")
if trouve:
    a, b, m, p, nb, temoins = trouve
    F  = deux**a * trois**(b+1)
    sp = str(p)
    print(f"  PREMIER DE {nb} CHIFFRES — PROUVÉ ✅")
    print(f"{'═'*65}")
    print(f"  a = {a},  b = {b}")
    print(f"  F = 2^{a}·3^{{{b+1}}}  ({gmpy2.num_digits(F)} chiffres) > √p  ✅")
    print(f"  q=2 → w={temoins[2]}  ✅     q=3 → w={temoins[3]}  ✅")
    print(f"\n  Temps total  : {t_total:.1f}s  ({t_total/60:.1f} min)")
    print(f"  Tests MR     : {n_mr}")
    print(f"  Crible élim. : {n_crible}")
    print(f"\n  p = {sp[:40]}…{sp[-20:]}")

    nom = f"premier_{nb}chiffres_zigzag.txt"
    with open(nom, "w") as f:
        f.write(f"Premier de {nb} chiffres — zigzag\n")
        f.write(f"a={a}, b={b}\n")
        f.write(f"Témoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\n")
        f.write(f"Temps   : {t_total:.1f}s\n\nm =\n{m}\n\np =\n{p}\n")
    print(f"\n  Sauvegardé → {nom}")

else:
    print(f"  COMPTE RENDU — Timeout {TIMEOUT_S/3600:.1f}h atteint")
    print(f"{'═'*65}")
    print(f"  Durée           : {t_total/3600:.2f}h")
    print(f"  a exploré       : centre={a_centre} ± {abs(a - a_centre)}")
    print(f"  Paires testées  : {n_paires:,}")
    print(f"  Élim. par crible: {n_crible:,}  ({n_crible/max(n_paires,1)*100:.1f}%)")
    print(f"  Tests MR        : {n_mr}")
    print(f"  Résultat        : aucun premier trouvé")
    print(f"\n  Diagnostic :")
    print(f"  → Chaque MR prend ~{t_total/max(n_mr,1):.0f}s pour {DIGITS_CIBLE} chiffres")
    print(f"  → En 4h on ne peut faire que ~{int(TIMEOUT_S/(t_total/max(n_mr,1)))} tests MR")
    print(f"  → Il faut ~{int(46052*0.125)} tests en moyenne pour trouver un premier")
    print(f"  → Temps nécessaire estimé : {int(46052*0.125)*(t_total/max(n_mr,1))/3600:.0f}h sur 1 cœur")
