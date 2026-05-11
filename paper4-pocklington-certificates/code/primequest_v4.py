#!/usr/bin/env python3
"""
PrimeQuest v4 — MULTI-CŒURS + FILTRES ARITHMÉTIQUES
Famille : p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Preuve  : Théorème de Pocklington N-1

Optimisations :
  1. Crible restreint q ≡ 1 (mod 3)
  2. Filtre Théorème 2 — (a%3, b%6) ∈ FORBIDDEN_T2 ⇒ 7|p
  3. Centre corrigé — ratio b/a = 1.0
  4. Parallélisation + checkpoint
"""

import gmpy2
import math
import time
import sys
import json
import os
import signal
import multiprocessing as mp

DIGITS_CIBLE    = 30_000
TIMEOUT_S       = 14_400
TOLERANCE       = 10
MR_TOURS        = 25
TEMOINS_MAX     = 300
SIEVE_LIMIT     = 1_000_000
RAPPORT_S       = 300
CHECKPOINT_S    = 60
CHECKPOINT_FILE = f"checkpoint_{DIGITS_CIBLE}.json"
RATIO_B_SUR_A   = 1.000
N_WORKERS_PARAM = 0
DELTA_DEPART    = 0
TEMPS_CUMUL_H   = 0.0

LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)
FORBIDDEN_T2 = frozenset({(0,2),(0,3),(1,0),(1,1),(2,4),(2,5)})

def _eratosthene_filtre(lim):
    t = bytearray([1]) * (lim + 1)
    t[0] = t[1] = 0
    for i in range(2, int(lim**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    return [i for i in range(5, lim + 1) if t[i] and i % 3 == 1]

_W_PREMIERS  = None
_W_FORBIDDEN = None

def _init_worker(premiers, forbidden):
    global _W_PREMIERS, _W_FORBIDDEN
    _W_PREMIERS  = premiers
    _W_FORBIDDEN = forbidden
    signal.signal(signal.SIGINT,  signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_IGN)

def _crible_mod(a, b):
    for q in _W_PREMIERS:
        x     = (pow(2, a % (q - 1), q) * pow(3, b % (q - 1), q)) % q
        p_mod = (3 * (x - 1) * x + 1) % q
        if p_mod == 0:
            return False
    return True

def _tester_paire(args):
    a, b, mr_tours = args
    if (a % 3, b % 6) in _W_FORBIDDEN:
        return {'statut': 'theoreme2', 'a': a, 'b': b}
    if not _crible_mod(a, b):
        return {'statut': 'crible', 'a': a, 'b': b}
    deux  = gmpy2.mpz(2)
    trois = gmpy2.mpz(3)
    m = deux**a * trois**b - 1
    p = 3 * m * (m + 1) + 1
    if gmpy2.is_prime(p, mr_tours):
        return {'statut': 'probable', 'a': a, 'b': b, 'p': p, 'm': m}
    return {'statut': 'compose', 'a': a, 'b': b}

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
            if v[q] is not None: continue
            val = gmpy2.powmod(gw, p1 // gmpy2.mpz(q), p)
            if gmpy2.gcd(val - 1, p) == 1:
                v[q] = w
        if all(x is not None for x in v.values()):
            return True, v
    return False, v

def a_depuis_position(centre, delta, side):
    if delta == 0: return centre
    return centre + delta if side == 0 else centre - delta

def position_suivante(delta, side, centre, a_min, a_max):
    if delta == 0:
        if centre + 1 <= a_max: return 1, 0
        elif centre - 1 >= a_min: return 1, 1
        else: return None, None
    if side == 0:
        if centre - delta >= a_min: return delta, 1
        else: return position_suivante(delta + 1, -1, centre, a_min, a_max)
    else:
        nd = delta + 1
        if centre + nd <= a_max: return nd, 0
        elif centre - nd >= a_min: return nd, 1
        else: return None, None

def generer_batch(delta, side, centre, a_min, a_max, digits, tol, taille):
    paires = []
    d, s = delta, side
    while len(paires) < taille and d is not None:
        a = a_depuis_position(centre, d, s)
        for b in b_valides(a, digits, tol):
            paires.append((a, b, MR_TOURS))
        d, s = position_suivante(d, s, centre, a_min, a_max)
    return paires, d, s

def sauver_checkpoint(delta, side, n_t2, n_crible, n_mr, n_paires, t_cumul):
    with open(CHECKPOINT_FILE, "w", encoding="utf-8") as f:
        json.dump({"digits": DIGITS_CIBLE, "delta": delta, "side": side,
                   "n_t2": n_t2, "n_crible": n_crible, "n_mr": n_mr,
                   "n_paires": n_paires, "t_cumul": t_cumul}, f, indent=2)

def charger_checkpoint():
    if not os.path.exists(CHECKPOINT_FILE): return None
    try:
        with open(CHECKPOINT_FILE, encoding="utf-8") as f:
            data = json.load(f)
        if data.get("digits") != DIGITS_CIBLE: return None
        return data
    except Exception: return None

def supprimer_checkpoint():
    if os.path.exists(CHECKPOINT_FILE): os.remove(CHECKPOINT_FILE)

if __name__ == '__main__':
    mp.freeze_support()
    n_coeurs  = os.cpu_count() or 4
    N_WORKERS = N_WORKERS_PARAM if N_WORKERS_PARAM > 0 else max(1, n_coeurs - 1)
    BATCH     = N_WORKERS * 3

    print(f"Construction du crible jusqu'à {SIEVE_LIMIT:,}…", end=" ", flush=True)
    PREMIERS = _eratosthene_filtre(SIEVE_LIMIT)
    print(f"{len(PREMIERS):,} premiers")

    a_centre = int((DIGITS_CIBLE / 2) / (LOG10_2 + RATIO_B_SUR_A * LOG10_3))
    a_max    = int(DIGITS_CIBLE / (2 * LOG10_2)) + 10

    cp = charger_checkpoint()
    if cp:
        delta = cp["delta"]; side = cp["side"]
        n_t2 = cp.get("n_t2", 0); n_crible = cp["n_crible"]
        n_mr = cp["n_mr"]; n_paires = cp["n_paires"]
        t_cumul = cp["t_cumul"]
        print(f"\n  ♻  REPRISE : delta={delta}")
    else:
        delta = side = 0
        n_t2 = n_crible = n_mr = n_paires = 0
        t_cumul = 0.0
        print(f"\n  Nouvelle recherche v4")

    print(f"  PrimeQuest v4 — {DIGITS_CIBLE:,} ch. | {N_WORKERS} workers\n")

    trouve   = None
    t_debut  = time.perf_counter()
    t_ckpt   = t_debut

    def _handler(sig, frame):
        sauver_checkpoint(delta, side, n_t2, n_crible, n_mr, n_paires,
                          t_cumul + (time.perf_counter() - t_debut))
        sys.exit(0)
    signal.signal(signal.SIGINT, _handler)
    signal.signal(signal.SIGTERM, _handler)

    with mp.Pool(N_WORKERS, initializer=_init_worker,
                 initargs=(PREMIERS, FORBIDDEN_T2)) as pool:
        while delta is not None:
            now = time.perf_counter()
            elapsed_s = now - t_debut
            elapsed_total = t_cumul + elapsed_s
            if elapsed_s >= TIMEOUT_S:
                sauver_checkpoint(delta, side, n_t2, n_crible, n_mr, n_paires, elapsed_total)
                break
            if now - t_ckpt >= CHECKPOINT_S:
                sauver_checkpoint(delta, side, n_t2, n_crible, n_mr, n_paires, elapsed_total)
                t_ckpt = now
            batch, next_delta, next_side = generer_batch(
                delta, side, a_centre, 1, a_max, DIGITS_CIBLE, TOLERANCE, BATCH)
            if not batch: break
            resultats = pool.map(_tester_paire, batch)
            for r in resultats:
                n_paires += 1
                if r['statut'] == 'theoreme2': n_t2 += 1
                elif r['statut'] == 'crible': n_crible += 1
                elif r['statut'] == 'compose':
                    n_mr += 1
                    et = t_cumul + (time.perf_counter() - t_debut)
                    print(f"  [{et/3600:5.2f}h]  MR #{n_mr:4d}  a={r['a']}, b={r['b']}  → composé", flush=True)
                elif r['statut'] == 'probable':
                    n_mr += 1
                    print("  *** PROBABLE — Pocklington… ***", flush=True)
                    ok, temoins = pocklington(r['p'], r['a'], r['b'])
                    if ok:
                        trouve = (r['a'], r['b'], r['m'], r['p'],
                                  int(gmpy2.num_digits(r['p'])), temoins)
                        break
            if trouve: break
            delta, side = next_delta, next_side

    t_total = t_cumul + (time.perf_counter() - t_debut)
    if trouve:
        supprimer_checkpoint()
        a, b, m, p, nb, temoins = trouve
        deux  = gmpy2.mpz(2); trois = gmpy2.mpz(3)
        F = deux**a * trois**(b + 1)
        print(f"  PREMIER {nb} ch.  a={a}, b={b}  ratio={b/a:.4f}  ✅")
        nom = f"premier_{nb}chiffres.txt"
        with open(nom, "w", encoding="utf-8") as f:
            f.write(f"Premier de {nb} chiffres — PrimeQuest v4\na={a}, b={b}  (ratio={b/a:.4f})\n")
            f.write(f"Témoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\nTemps : {t_total:.1f}s\nWorkers : {N_WORKERS}\n\nm =\n{m}\n\np =\n{p}\n")
        print(f"  Sauvegardé → {nom}")
    else:
        print(f"  Session terminée — {CHECKPOINT_FILE}")
