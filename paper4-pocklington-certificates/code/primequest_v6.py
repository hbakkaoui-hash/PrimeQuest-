#!/usr/bin/env python3
"""
PrimeQuest v6 — Exploration Multi-Ratio b/a
Famille : p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Preuve  : Théorème de Pocklington N-1

Optimisations cumulées (v1 → v6) :
  1. Crible restreint aux q ≡ 1 (mod 3)
  2. Filtrage Théorème 2 (7|p gratuit)
  3. Parallélisation multi-cœurs
  4. Checkpoint 60s
  5. MR_TOURS = 3
  6. SIEVE_LIMIT = 10 000 000
  [v6] 7. RATIOS_LIST : exploration simultanée de plusieurs ratios b/a
     Chaque ratio définit un centre (a₀, b₀) distinct dans l'espace (a,b).
     Lors de chaque batch, les paires sont générées en round-robin.
     Ratios par défaut : [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15]
"""

import gmpy2
import math
import time
import sys
import json
import os
import signal
import multiprocessing as mp

DIGITS_CIBLE     = 50_000
TIMEOUT_S        = 14_400
TOLERANCE        = 10
MR_TOURS         = 3
MR_TOURS_CONFIRM = 20
TEMOINS_MAX      = 300
SIEVE_LIMIT      = 10_000_000
RAPPORT_S        = 300
CHECKPOINT_S     = 60
CHECKPOINT_FILE  = f"checkpoint_{DIGITS_CIBLE}_v6.json"

N_RATIOS      = 7
RATIO_CENTRE  = 1.00
RATIO_SPREAD  = 0.15

N_WORKERS_PARAM = 0

LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)

_RATIO_MIN_ABSOLU = 2 * LOG10_2 / DIGITS_CIBLE
_RATIO_MAX_ABSOLU = (DIGITS_CIBLE / 2 - LOG10_2) / LOG10_3

def _generer_ratios():
    step = (2 * RATIO_SPREAD / (N_RATIOS - 1)) if N_RATIOS > 1 else 0
    candidats = [
        round(RATIO_CENTRE - RATIO_SPREAD + i * step, 4)
        for i in range(N_RATIOS)
    ]
    valides = []
    for r in candidats:
        if r < _RATIO_MIN_ABSOLU or r > _RATIO_MAX_ABSOLU:
            continue
        ac = max(1, int((DIGITS_CIBLE / 2) / (LOG10_2 + r * LOG10_3)))
        bc = int(ac * r)
        if bc < 1:
            continue
        valides.append(r)
    if not valides:
        raise ValueError("Aucun ratio valide.")
    return valides

RATIOS_LIST = _generer_ratios()

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
    a, b, mr_tours, mr_confirm = args
    if (a % 3, b % 6) in _W_FORBIDDEN:
        return {'statut': 'theoreme2', 'a': a, 'b': b}
    if not _crible_mod(a, b):
        return {'statut': 'crible', 'a': a, 'b': b}
    deux  = gmpy2.mpz(2)
    trois = gmpy2.mpz(3)
    m = deux**a * trois**b - 1
    p = 3 * m * (m + 1) + 1
    if not gmpy2.is_prime(p, mr_tours):
        return {'statut': 'compose', 'a': a, 'b': b}
    if not gmpy2.is_prime(p, mr_confirm):
        return {'statut': 'compose', 'a': a, 'b': b}
    return {'statut': 'probable', 'a': a, 'b': b, 'p': p, 'm': m}

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

def generer_batch_multi(deltas, sides, centres, a_max, digits, tol, taille):
    paires     = []
    n          = len(centres)
    new_deltas = list(deltas)
    new_sides  = list(sides)
    while len(paires) < taille:
        progress = False
        for i in range(n):
            if new_deltas[i] is None:
                continue
            progress = True
            ac = centres[i][0]
            a  = a_depuis_position(ac, new_deltas[i], new_sides[i])
            for b in b_valides(a, digits, tol):
                paires.append((a, b, MR_TOURS, MR_TOURS_CONFIRM))
            nd, ns        = position_suivante(new_deltas[i], new_sides[i], ac, 1, a_max)
            new_deltas[i] = nd
            new_sides[i]  = ns
            if len(paires) >= taille:
                break
        if not progress:
            break
    return paires, new_deltas, new_sides

def sauver_checkpoint(deltas, sides, n_t2, n_crible, n_mr, n_paires, t_cumul):
    with open(CHECKPOINT_FILE, "w", encoding="utf-8") as f:
        json.dump({
            "digits":   DIGITS_CIBLE,
            "version":  6,
            "ratios":   RATIOS_LIST,
            "deltas":   [d if d is not None else -1 for d in deltas],
            "sides":    sides,
            "n_t2":     n_t2,
            "n_crible": n_crible,
            "n_mr":     n_mr,
            "n_paires": n_paires,
            "t_cumul":  t_cumul,
        }, f, indent=2)

def charger_checkpoint():
    if not os.path.exists(CHECKPOINT_FILE):
        return None
    try:
        with open(CHECKPOINT_FILE, encoding="utf-8") as f:
            data = json.load(f)
        if data.get("digits") != DIGITS_CIBLE: return None
        if data.get("version") != 6: return None
        if data.get("ratios") != RATIOS_LIST: return None
        return data
    except:
        return None

def supprimer_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        os.remove(CHECKPOINT_FILE)


if __name__ == '__main__':
    mp.freeze_support()
    n_coeurs  = os.cpu_count() or 4
    N_WORKERS = N_WORKERS_PARAM if N_WORKERS_PARAM > 0 else max(1, n_coeurs - 1)
    BATCH     = N_WORKERS * 3

    print(f"Construction du crible jusqu'à {SIEVE_LIMIT:,}…", end=" ", flush=True)
    PREMIERS = _eratosthene_filtre(SIEVE_LIMIT)
    print(f"{len(PREMIERS):,} premiers utiles")

    a_max   = int(DIGITS_CIBLE / (2 * LOG10_2)) + 10
    centres = []
    for r in RATIOS_LIST:
        ac = int((DIGITS_CIBLE / 2) / (LOG10_2 + r * LOG10_3))
        bc = int(ac * r)
        centres.append((ac, bc))

    N_RATIOS_ACTIFS = len(RATIOS_LIST)

    cp = charger_checkpoint()
    if cp:
        deltas   = [None if d == -1 else d for d in cp["deltas"]]
        sides    = cp["sides"]
        n_t2     = cp.get("n_t2", 0)
        n_crible = cp["n_crible"]
        n_mr     = cp["n_mr"]
        n_paires = cp["n_paires"]
        t_cumul  = cp["t_cumul"]
        print(f"  ♻  REPRISE v6 : MR={n_mr}  cumulé={t_cumul/3600:.2f}h")
    else:
        deltas   = [0] * N_RATIOS_ACTIFS
        sides    = [0] * N_RATIOS_ACTIFS
        n_t2 = n_crible = n_mr = n_paires = 0
        t_cumul = 0.0
        print("  Nouvelle recherche v6")

    print(f"  Ratios actifs : {RATIOS_LIST}")
    print(f"  Centres (a₀,b₀) : {centres}")

    trouve  = None
    t_debut = time.perf_counter()
    t_ckpt  = t_debut
    t_rapp  = t_debut

    def _handler(sig, frame):
        sauver_checkpoint(deltas, sides, n_t2, n_crible, n_mr, n_paires,
                          t_cumul + (time.perf_counter() - t_debut))
        sys.exit(0)

    signal.signal(signal.SIGINT,  _handler)
    signal.signal(signal.SIGTERM, _handler)

    with mp.Pool(N_WORKERS, initializer=_init_worker,
                 initargs=(PREMIERS, FORBIDDEN_T2)) as pool:

        while any(d is not None for d in deltas):
            now       = time.perf_counter()
            elapsed_s = now - t_debut
            elapsed_t = t_cumul + elapsed_s

            if elapsed_s >= TIMEOUT_S:
                sauver_checkpoint(deltas, sides, n_t2, n_crible, n_mr, n_paires, elapsed_t)
                break
            if now - t_ckpt >= CHECKPOINT_S:
                sauver_checkpoint(deltas, sides, n_t2, n_crible, n_mr, n_paires, elapsed_t)
                t_ckpt = now
            if now - t_rapp >= RAPPORT_S:
                elim_tot = n_t2 + n_crible
                print(f"  [{elapsed_t/3600:5.2f}h]  paires={n_paires:,}  MR={n_mr}  "
                      f"élim={elim_tot/max(n_paires,1)*100:.1f}%", flush=True)
                t_rapp = now

            batch, deltas, sides = generer_batch_multi(
                deltas, sides, centres, a_max, DIGITS_CIBLE, TOLERANCE, BATCH)
            if not batch: break

            resultats = pool.map(_tester_paire, batch)

            for r in resultats:
                n_paires += 1
                if r['statut'] == 'theoreme2': n_t2 += 1
                elif r['statut'] == 'crible': n_crible += 1
                elif r['statut'] == 'compose':
                    n_mr += 1
                    print(f"  MR #{n_mr:4d}  a={r['a']}, b={r['b']}  r={r['b']/r['a']:.3f}  → composé", flush=True)
                elif r['statut'] == 'probable':
                    n_mr += 1
                    print("  *** PROBABLE — Pocklington… ***", flush=True)
                    ok, temoins = pocklington(r['p'], r['a'], r['b'])
                    if ok:
                        trouve = (r['a'], r['b'], r['m'], r['p'],
                                  int(gmpy2.num_digits(r['p'])), temoins)
                        break
            if trouve: break

    if trouve:
        supprimer_checkpoint()
        a, b, m, p, nb, temoins = trouve
        sp = str(p)
        print(f"  PREMIER DE {nb} CHIFFRES ✅")
        print(f"  q=2 → w={temoins[2]}  |  q=3 → w={temoins[3]}")
        print(f"  p = {sp[:40]}…{sp[-20:]}")
        nom = f"premier_{nb}chiffres_v6.txt"
        with open(nom, "w", encoding="utf-8") as f:
            f.write(f"Premier de {nb} chiffres — PrimeQuest v6\n")
            f.write(f"a={a}, b={b}  ratio={b/a:.4f}\n")
            f.write(f"Témoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\n")
            f.write(f"Ratios : {RATIOS_LIST}\n\nm =\n{m}\n\np =\n{p}\n")
        print(f"  Sauvegardé → {nom}")
