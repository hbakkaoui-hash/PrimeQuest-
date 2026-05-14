#!/usr/bin/env python3
"""
PrimeQuest v2 — Recherche et preuve de primalité AVEC REPRISE
Famille : p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Preuve  : Théorème de Pocklington N-1

Optimisations intégrées :
  1. Crible modulaire restreint aux premiers q ≡ 1 (mod 3)
  2. Zigzag depuis a_centre = (D/2)/log10(6)
  3. Checkpoint automatique — reprise exacte après interruption
"""

import gmpy2
import math
import time
import sys
import json
import os
import signal

# ═══════════════════════════════════════════════════════════════
DIGITS_CIBLE    = 20_000
TIMEOUT_S       = 14_400
TOLERANCE       = 10
MR_TOURS        = 25
TEMOINS_MAX     = 300
SIEVE_LIMIT     = 1_000_000
RAPPORT_S       = 300
CHECKPOINT_S    = 60
CHECKPOINT_FILE = f"checkpoint_{DIGITS_CIBLE}.json"
DELTA_DEPART    = 0
TEMPS_CUMUL_H   = 0.0
# ═══════════════════════════════════════════════════════════════

deux    = gmpy2.mpz(2)
trois   = gmpy2.mpz(3)
LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)

def _eratosthene_filtre(lim):
    t = bytearray([1]) * (lim + 1)
    t[0] = t[1] = 0
    for i in range(2, int(lim**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    return [i for i in range(5, lim + 1) if t[i] and i % 3 == 1]

print(f"Construction du crible (q ≡ 1 mod 3) jusqu'à {SIEVE_LIMIT:,}…",
      end=" ", flush=True)
t0 = time.perf_counter()
PREMIERS = _eratosthene_filtre(SIEVE_LIMIT)
print(f"{len(PREMIERS):,} premiers utiles  ({time.perf_counter()-t0:.2f}s)")

def crible_modulaire(a, b):
    for q in PREMIERS:
        x     = (pow(2, a % (q - 1), q) * pow(3, b % (q - 1), q)) % q
        p_mod = (3 * (x - 1) * x + 1) % q
        if p_mod == 0:
            return False
    return True

def b_valides(a, digits, tol):
    b_c = (digits / 2 - a * LOG10_2 - LOG10_3 / 2) / LOG10_3
    return [
        b for b in range(max(1, int(b_c) - 3), int(b_c) + 4)
        if abs(2 * (a * LOG10_2 + b * LOG10_3) + LOG10_3 - digits) <= tol + 2
    ]

def pocklington(p, a, b):
    p1 = p - 1
    v  = {2: None, 3: None}
    for w in range(2, 2 + TEMOINS_MAX):
        gw = gmpy2.mpz(w)
        if gmpy2.powmod(gw, p1, p) != 1:
            continue
        for q in [2, 3]:
            if v[q] is not None:
                continue
            val = gmpy2.powmod(gw, p1 // gmpy2.mpz(q), p)
            if gmpy2.gcd(val - 1, p) == 1:
                v[q] = w
        if all(x is not None for x in v.values()):
            return True, v
    return False, v

def a_depuis_position(centre, delta, side):
    if delta == 0:
        return centre
    return centre + delta if side == 0 else centre - delta

def position_suivante(delta, side, centre, a_min, a_max):
    if delta == 0:
        if centre + 1 <= a_max:
            return 1, 0
        elif centre - 1 >= a_min:
            return 1, 1
        else:
            return None, None
    if side == 0:
        if centre - delta >= a_min:
            return delta, 1
        else:
            return position_suivante(delta + 1, -1, centre, a_min, a_max)
    else:
        nd = delta + 1
        if centre + nd <= a_max:
            return nd, 0
        elif centre - nd >= a_min:
            return nd, 1
        else:
            return None, None

def sauver_checkpoint(delta, side, n_crible, n_mr, n_paires, t_cumul):
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump({"digits": DIGITS_CIBLE, "delta": delta, "side": side,
                   "n_crible": n_crible, "n_mr": n_mr, "n_paires": n_paires,
                   "t_cumul": t_cumul}, f, indent=2)

def charger_checkpoint():
    if not os.path.exists(CHECKPOINT_FILE):
        return None
    try:
        with open(CHECKPOINT_FILE) as f:
            data = json.load(f)
        if data.get("digits") != DIGITS_CIBLE:
            return None
        return data
    except Exception:
        return None

def supprimer_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        os.remove(CHECKPOINT_FILE)

_etat_global = {}

def handler_signal(sig, frame):
    print(f"\n  ⚡ Interruption — sauvegarde du checkpoint…", flush=True)
    if _etat_global:
        sauver_checkpoint(
            _etat_global["delta"], _etat_global["side"],
            _etat_global["n_crible"], _etat_global["n_mr"],
            _etat_global["n_paires"],
            _etat_global["t_cumul"] + (time.perf_counter() - _etat_global["t_debut"])
        )
    sys.exit(0)

signal.signal(signal.SIGINT,  handler_signal)
signal.signal(signal.SIGTERM, handler_signal)

a_centre = int((DIGITS_CIBLE / 2) / math.log10(6))
a_max    = int(DIGITS_CIBLE / (2 * LOG10_2)) + 10

if DELTA_DEPART > 0:
    delta = DELTA_DEPART; side = 0
    n_crible = n_mr = n_paires = 0
    t_cumul  = TEMPS_CUMUL_H * 3600
else:
    cp = charger_checkpoint()
    if cp:
        delta = cp["delta"]; side = cp["side"]
        n_crible = cp["n_crible"]; n_mr = cp["n_mr"]
        n_paires = cp["n_paires"]; t_cumul = cp["t_cumul"]
        print(f"\n  ♻  REPRISE : delta={delta}, a={a_depuis_position(a_centre, delta, side)}")
    else:
        delta = side = 0
        n_crible = n_mr = n_paires = 0
        t_cumul  = 0.0
        print(f"\n  Nouvelle recherche (aucun checkpoint)")

print(f"\n{'═'*65}")
print(f"  PrimeQuest v2 — p = 3m(m+1)+1,  m = 2^a·3^b−1")
print(f"{'═'*65}")
print(f"  Cible : {DIGITS_CIBLE:,} chiffres (±{TOLERANCE})")
print(f"  a_centre : {a_centre:,}   Timeout : {TIMEOUT_S/3600:.1f}h")
print(f"{'─'*65}\n")

trouve   = None
t_debut  = time.perf_counter()
t_rapport = t_debut
t_ckpt   = t_debut

_etat_global.update({"delta": delta, "side": side, "n_crible": n_crible,
                     "n_mr": n_mr, "n_paires": n_paires,
                     "t_cumul": t_cumul, "t_debut": t_debut})

while delta is not None:
    elapsed_session = time.perf_counter() - t_debut
    elapsed_total   = t_cumul + elapsed_session

    if elapsed_session >= TIMEOUT_S:
        sauver_checkpoint(delta, side, n_crible, n_mr, n_paires, elapsed_total)
        print(f"  ⏱  Timeout.  Checkpoint → {CHECKPOINT_FILE}")
        break

    a = a_depuis_position(a_centre, delta, side)

    for b in b_valides(a, DIGITS_CIBLE, TOLERANCE):
        n_paires += 1
        now = time.perf_counter()
        if now - t_ckpt >= CHECKPOINT_S:
            sauver_checkpoint(delta, side, n_crible, n_mr, n_paires,
                              t_cumul + (now - t_debut))
            t_ckpt = now
            _etat_global.update({"delta": delta, "side": side,
                                 "n_crible": n_crible, "n_mr": n_mr,
                                 "n_paires": n_paires,
                                 "t_cumul": t_cumul, "t_debut": t_debut})
        if not crible_modulaire(a, b):
            n_crible += 1
            continue
        n_mr += 1
        m = deux**a * trois**b - 1
        p = 3 * m * (m + 1) + 1
        et = t_cumul + (time.perf_counter() - t_debut)
        print(f"  [{et/3600:5.2f}h]  MR #{n_mr:4d}  a={a}, b={b} …", flush=True)
        if not gmpy2.is_prime(p, MR_TOURS):
            print("  → composé", flush=True)
            continue
        print("  *** PROBABLE — Pocklington… ***", flush=True)
        ok, temoins = pocklington(p, a, b)
        if ok:
            trouve = (a, b, m, p, int(gmpy2.num_digits(p)), temoins)
            break

    if trouve:
        break
    delta, side = position_suivante(delta, side, a_centre, 1, a_max)
    _etat_global["delta"] = delta
    _etat_global["side"]  = side

t_total = t_cumul + (time.perf_counter() - t_debut)
print(f"\n{'═'*65}")

if trouve:
    supprimer_checkpoint()
    a, b, m, p, nb, temoins = trouve
    F = deux**a * trois**(b + 1)
    print(f"  PREMIER DE {nb} CHIFFRES — PRIMALITÉ PROUVÉE ✅")
    print(f"  a={a}, b={b}  F={gmpy2.num_digits(F)} chiffres > √p  ✅")
    print(f"  q=2 → w={temoins[2]}   q=3 → w={temoins[3]}")
    nom = f"premier_{nb}chiffres.txt"
    with open(nom, "w") as f:
        f.write(f"Premier de {nb} chiffres — PrimeQuest v2\na={a}, b={b}\n")
        f.write(f"Témoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\n")
        f.write(f"Temps   : {t_total:.1f}s\n\nm =\n{m}\n\np =\n{p}\n")
    print(f"  Sauvegardé → {nom}")
else:
    print(f"  Session terminée — Checkpoint → {CHECKPOINT_FILE}")
