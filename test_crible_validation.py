#!/usr/bin/env python3
"""
Validation du crible modulaire + recherche d'un premier de 20 000 chiffres
Deux phases :
  1. Validation : compare crible modulaire vs crible explicite sur N paires
  2. Recherche  : streaming (sans pré-énumération) + crible modulaire
"""

import gmpy2
import math
import time
import sys

DIGITS_CIBLE  = 10_000
TOLERANCE     = 10
MR_TOURS      = 25
TEMOINS_MAX   = 300
SIEVE_LIMIT   = 1_000_000   # réduit pour les tests (78 498 premiers)
N_VALIDATION  = 2_000       # paires à comparer pour la validation

deux  = gmpy2.mpz(2)
trois = gmpy2.mpz(3)
_LOG10_2 = math.log10(2)
_LOG10_3 = math.log10(3)

# ── Crible ───────────────────────────────────────────────────────────────────

def _eratosthene(limite):
    t = bytearray([1]) * (limite + 1)
    t[0] = t[1] = 0
    for i in range(2, int(limite**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    return [i for i in range(2, limite + 1) if t[i]]

print(f"Construction du crible jusqu'à {SIEVE_LIMIT:,}…", end=" ", flush=True)
t0 = time.perf_counter()
PREMIERS = _eratosthene(SIEVE_LIMIT)
print(f"{len(PREMIERS):,} premiers  ({time.perf_counter()-t0:.2f}s)")

# ── Méthode A : crible EXPLICITE (calcule p en grande précision) ──────────────

def crible_explicite(a, b):
    """Référence : calcule p entier puis teste p % q pour chaque petit premier."""
    m = deux**a * trois**b - 1
    p = 3 * m * (m + 1) + 1
    for q in PREMIERS:
        if q * q > p:
            return True
        if p % q == 0:
            return p == q
    return True

# ── Méthode B : crible MODULAIRE (sans calculer p) ───────────────────────────

def crible_modulaire(a, b):
    """
    Teste p mod q SANS calculer p.
    p = 3·(2^a·3^b−1)·(2^a·3^b) + 1
    p mod q = (3·(x−1)·x + 1) mod q  où  x = (2^a · 3^b) mod q
    Petit théorème de Fermat : 2^a mod q = 2^(a mod (q−1)) mod q  (q premier impair)
    """
    for q in PREMIERS:
        if q <= 3:
            continue   # p ≡ 1 (mod 2) et p ≡ 1 (mod 3) toujours
        a_red  = a % (q - 1)
        b_red  = b % (q - 1)
        pow2a  = pow(2, a_red, q)
        pow3b  = pow(3, b_red, q)
        x      = (pow2a * pow3b) % q
        p_mod  = (3 * (x - 1) * x + 1) % q
        if p_mod == 0:
            return False
    return True

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 1 : VALIDATION — accord crible explicite vs modulaire
# ═══════════════════════════════════════════════════════════════════════════════

print(f"\n{'─'*60}")
print(f"PHASE 1 — Validation sur {N_VALIDATION} paires (a,b)")
print(f"{'─'*60}")

# Paires autour de 1 000 chiffres (rapides à calculer explicitement)
DIGITS_TEST = 1_000

def b_approx(a, digits):
    b_c = (digits / 2 - a * _LOG10_2 - _LOG10_3 / 2) / _LOG10_3
    return [b for b in range(max(1, int(b_c) - 2), int(b_c) + 3)]

accord      = 0
desaccord   = 0
total_test  = 0
erreurs     = []

t_val = time.perf_counter()

a = 1
while total_test < N_VALIDATION:
    for b in b_approx(a, DIGITS_TEST):
        if total_test >= N_VALIDATION:
            break

        # Vérifier que p est bien dans la bonne fenêtre
        nb_approx = 2 * (a * _LOG10_2 + b * _LOG10_3) + _LOG10_3
        if abs(nb_approx - DIGITS_TEST) > 10:
            continue

        res_explicite = crible_explicite(a, b)
        res_modulaire = crible_modulaire(a, b)
        total_test   += 1

        if res_explicite == res_modulaire:
            accord += 1
        else:
            desaccord += 1
            erreurs.append((a, b, res_explicite, res_modulaire))
    a += 1

t_val = time.perf_counter() - t_val

print(f"Paires testées  : {total_test}")
print(f"Accords         : {accord}  ({accord/total_test*100:.1f}%)")
print(f"Désaccords      : {desaccord}")
print(f"Temps           : {t_val:.2f}s")

if desaccord == 0:
    print(f"\n✅ CRIBLE MODULAIRE VALIDÉ — accord parfait sur {total_test} paires")
else:
    print(f"\n❌ ERREURS DÉTECTÉES :")
    for a, b, exp, mod in erreurs[:10]:
        print(f"   a={a}, b={b} : explicite={exp}, modulaire={mod}")
    sys.exit(1)

# Validation supplémentaire : le premier de 10 000 chiffres connu
print(f"\nTest du premier connu (a=10094, b=4110) :")
r_exp = crible_explicite(10094, 4110)
r_mod = crible_modulaire(10094, 4110)
print(f"  Crible explicite  : {'PASSE ✅' if r_exp else 'ÉLIMINÉ ❌'}")
print(f"  Crible modulaire  : {'PASSE ✅' if r_mod else 'ÉLIMINÉ ❌'}")
if r_exp == r_mod:
    print(f"  Accord ✅")
else:
    print(f"  DÉSACCORD ❌ — bug !")
    sys.exit(1)

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 2 : RECHERCHE STREAMING — 20 000 chiffres
# Pas de pré-énumération : on teste chaque (a,b) à la volée
# ═══════════════════════════════════════════════════════════════════════════════

print(f"\n{'─'*60}")
print(f"PHASE 2 — Recherche d'un premier de {DIGITS_CIBLE} chiffres")
print(f"Approche streaming (sans pré-énumération)")
print(f"{'─'*60}\n")

def pocklington(p, a, b):
    p1     = p - 1
    valide = {2: None, 3: None}
    for w in range(2, 2 + TEMOINS_MAX):
        gw = gmpy2.mpz(w)
        if gmpy2.powmod(gw, p1, p) != 1:
            continue
        for q in [2, 3]:
            if valide[q] is not None:
                continue
            val = gmpy2.powmod(gw, p1 // gmpy2.mpz(q), p)
            if gmpy2.gcd(val - 1, p) == 1:
                valide[q] = w
        if all(v is not None for v in valide.values()):
            return True, valide
    return False, valide

def b_streaming(a, digits, tol):
    b_c = (digits / 2 - a * _LOG10_2 - _LOG10_3 / 2) / _LOG10_3
    return [b for b in range(max(1, int(b_c) - 3), int(b_c) + 4)
            if abs(2*(a*_LOG10_2 + b*_LOG10_3) + _LOG10_3 - digits) <= tol + 2]

a_max    = int(DIGITS_CIBLE / (2 * _LOG10_2)) + 50
n_crible   = 0
n_mr       = 0
n_paires   = 0
trouve     = None
t_debut    = time.perf_counter()
t_rapport  = t_debut
RAPPORT_S  = 60   # affiche un rapport toutes les 60 secondes

for a in range(1, a_max):
    for b in b_streaming(a, DIGITS_CIBLE, TOLERANCE):
        n_paires += 1

        # ── Rapport de progression toutes les 60 secondes ────────────────
        maintenant = time.perf_counter()
        if maintenant - t_rapport >= RAPPORT_S:
            elapsed = maintenant - t_debut
            cadence = n_paires / elapsed if elapsed > 0 else 0
            print(f"  [+{elapsed:6.0f}s]  a={a:6d}  paires={n_paires:8,}"
                  f"  crible={n_crible:7,}  MR={n_mr:5,}"
                  f"  cadence={cadence:.0f}/s",
                  flush=True)
            t_rapport = maintenant

        # ── Crible modulaire (µs) ────────────────────────────────────────
        if not crible_modulaire(a, b):
            n_crible += 1
            continue

        # ── Miller-Rabin (calcul de p requis ici) ────────────────────────
        n_mr += 1
        elapsed = time.perf_counter() - t_debut
        print(f"  [+{elapsed:.1f}s]  MR #{n_mr}  a={a}, b={b} …", flush=True)
        m = deux**a * trois**b - 1
        p = 3 * m * (m + 1) + 1
        nb = int(gmpy2.num_digits(p))

        if not gmpy2.is_prime(p, MR_TOURS):
            print(f"  → composé", flush=True)
            continue

        # ── Candidat probable ────────────────────────────────────────────
        elapsed = time.perf_counter() - t_debut
        print(f"  Candidat (a={a}, b={b}) — {nb} chiffres "
              f"[{elapsed:.1f}s] — Pocklington…")

        ok, temoins = pocklington(p, a, b)
        if ok:
            trouve = (a, b, m, p, nb, temoins)
            break

    if trouve:
        break

t_total = time.perf_counter() - t_debut

if not trouve:
    print("❌ Aucun premier trouvé.")
    sys.exit(1)

a, b, m, p, nb, temoins = trouve
F = deux**a * trois**(b+1)
sp = str(p)

print(f"""
{'═'*65}
  PREMIER DE {nb} CHIFFRES — PRIMALITÉ PROUVÉE
{'═'*65}

  a = {a},  b = {b}
  m = 2^{a}·3^{b}−1   ({gmpy2.num_digits(m)} chiffres)
  p = 3m(m+1)+1    ({nb} chiffres)
  F = 2^{a}·3^{{{b+1}}} ({gmpy2.num_digits(F)} chiffres) > √p  ✅

  Témoins Pocklington :
    q=2  → w={temoins[2]}  ✅
    q=3  → w={temoins[3]}  ✅

  Performance
  ───────────────────────────────────────────
  Éliminés par crible modulaire : {n_crible}
  Tests Miller-Rabin             : {n_mr}
  Temps total                    : {t_total:.2f}s

  p = {sp[:40]}…{sp[-20:]}
""")

nom = f"premier_{nb}chiffres.txt"
with open(nom, "w") as f:
    f.write(f"Premier de {nb} chiffres\n")
    f.write(f"Méthode : Pocklington, p=3m(m+1)+1, streaming\n")
    f.write(f"a={a}, b={b}\n")
    f.write(f"Témoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\n")
    f.write(f"Temps   : {t_total:.2f}s\n\n")
    f.write(f"m =\n{m}\n\np =\n{p}\n")
print(f"  Sauvegardé → {nom}")
