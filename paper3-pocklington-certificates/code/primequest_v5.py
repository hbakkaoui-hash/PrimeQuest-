#!/usr/bin/env python3
"""
PrimeQuest v5 — Recherche 50 000 chiffres — OPTIMISE
Famille : p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Preuve  : Théorème de Pocklington N-1

Optimisations (en plus de v4) :
  5a. MR_TOURS 25 → 3 (initial) + 20 (confirmation) : ~8× plus rapide
      Miller-Rabin est déterministe pour les premiers (zéro faux négatifs).
      Les rares composites passant 3 tours (≤ 4^-3 ≈ 1.6%) sont arrêtés
      par les 20 tours de confirmation et le test de Pocklington.
  5b. SIEVE_LIMIT 10^6 → 10^7 : 664k premiers utiles (vs 78k), ~15% de MR en moins.
"""

import gmpy2
import math
import time
import sys
import json
import os
import signal
import multiprocessing as mp

# ═══════════════════════════════════════════════════════════════
# PARAMÈTRES
# ═══════════════════════════════════════════════════════════════
DIGITS_CIBLE      = 50_000      # nombre de chiffres souhaité
TIMEOUT_S         = 43_200      # 12h par session
TOLERANCE         = 10
MR_TOURS          = 3           # tours initiaux (rapides, déterministes pour p premier)
MR_TOURS_CONFIRM  = 20          # tours de confirmation avant Pocklington
TEMOINS_MAX       = 300
SIEVE_LIMIT       = 10_000_000  # crible étendu à 10M (664k premiers q≡1 mod3)
RAPPORT_S         = 300
CHECKPOINT_S      = 60
CHECKPOINT_FILE   = f"checkpoint_{DIGITS_CIBLE}.json"

RATIO_B_SUR_A     = 1.000

N_WORKERS_PARAM   = 0
DELTA_DEPART      = 0
TEMPS_CUMUL_H     = 0.0
# ═══════════════════════════════════════════════════════════════

LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)

F7_FORBIDDEN = frozenset({(0,2),(0,3),(1,0),(1,1),(2,4),(2,5)})


# ── Crible d'Ératosthène ─────────────────────────────────────────────────────
def _eratosthene_filtre(lim):
    t = bytearray([1]) * (lim + 1)
    t[0] = t[1] = 0
    for i in range(2, int(lim**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    return [i for i in range(5, lim + 1) if t[i] and i % 3 == 1]


# ── Fonctions worker ──────────────────────────────────────────────────────────
_W_PREMIERS = None

def _init_worker(premiers):
    global _W_PREMIERS
    _W_PREMIERS = premiers
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
    """
    Teste une paire (a, b) : filtre mod-7, crible, MR 3 tours, MR 20 tours.
    Retourne {'statut': 'filtre7'|'crible'|'compose'|'probable', ...}.
    """
    a, b, mr_tours, mr_confirm = args
    # Étape 1 : filtre mod-7 — O(1)
    if (a % 3, b % 6) in F7_FORBIDDEN:
        return {'statut': 'filtre7', 'a': a, 'b': b}
    # Étape 2 : crible (q ≡ 1 mod 3)
    if not _crible_mod(a, b):
        return {'statut': 'crible', 'a': a, 'b': b}
    # Étape 3 : Miller-Rabin initial (3 tours rapides)
    deux  = gmpy2.mpz(2)
    trois = gmpy2.mpz(3)
    m = deux**a * trois**b - 1
    p = 3 * m * (m + 1) + 1
    if not gmpy2.is_prime(p, mr_tours):
        return {'statut': 'compose', 'a': a, 'b': b}
    # Étape 4 : Miller-Rabin confirmation (20 tours)
    if not gmpy2.is_prime(p, mr_confirm):
        return {'statut': 'compose', 'a': a, 'b': b}
    return {'statut': 'probable', 'a': a, 'b': b, 'p': p, 'm': m}


# ── b valides pour un a donné ─────────────────────────────────────────────────
def b_valides(a, digits, tol):
    b_c = (digits / 2 - a * LOG10_2 - LOG10_3 / 2) / LOG10_3
    return [
        b for b in range(max(1, int(b_c) - 3), int(b_c) + 4)
        if abs(2 * (a * LOG10_2 + b * LOG10_3) + LOG10_3 - digits) <= tol + 2
    ]


# ── Pocklington ───────────────────────────────────────────────────────────────
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


# ── Zigzag depuis le centre ───────────────────────────────────────────────────
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

def generer_batch(delta, side, centre, a_min, a_max, digits, tol, taille):
    paires = []
    d, s = delta, side
    while len(paires) < taille and d is not None:
        a = a_depuis_position(centre, d, s)
        for b in b_valides(a, digits, tol):
            paires.append((a, b, MR_TOURS, MR_TOURS_CONFIRM))
        d, s = position_suivante(d, s, centre, a_min, a_max)
    return paires, d, s


# ── Checkpoint ────────────────────────────────────────────────────────────────
def sauver_checkpoint(delta, side, n_filtre7, n_crible, n_mr, n_paires, t_cumul):
    with open(CHECKPOINT_FILE, "w", encoding='utf-8') as f:
        json.dump({
            "digits":    DIGITS_CIBLE,
            "delta":     delta,
            "side":      side,
            "n_filtre7": n_filtre7,
            "n_crible":  n_crible,
            "n_mr":      n_mr,
            "n_paires":  n_paires,
            "t_cumul":   t_cumul,
        }, f, indent=2)

def charger_checkpoint():
    if not os.path.exists(CHECKPOINT_FILE):
        return None
    try:
        with open(CHECKPOINT_FILE, encoding='utf-8') as f:
            data = json.load(f)
        if data.get("digits") != DIGITS_CIBLE:
            print(f"  Checkpoint ignore (digits={data.get('digits')} != {DIGITS_CIBLE})")
            return None
        return data
    except Exception as e:
        print(f"  Checkpoint illisible : {e}")
        return None

def supprimer_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        os.remove(CHECKPOINT_FILE)


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
if __name__ == '__main__':

    mp.freeze_support()

    n_coeurs  = os.cpu_count() or 4
    N_WORKERS = N_WORKERS_PARAM if N_WORKERS_PARAM > 0 else max(1, n_coeurs - 1)
    BATCH     = N_WORKERS * 3

    print(f"Construction du crible (q ≡ 1 mod 3) jusqu'à {SIEVE_LIMIT:,}…",
          end=" ", flush=True)
    t0 = time.perf_counter()
    PREMIERS = _eratosthene_filtre(SIEVE_LIMIT)
    print(f"{len(PREMIERS):,} premiers utiles  ({time.perf_counter()-t0:.2f}s)")
    print(f"Coeurs disponibles : {n_coeurs}  →  {N_WORKERS} workers  "
          f"(batch {BATCH} paires)")

    a_centre = int((DIGITS_CIBLE - LOG10_3) / (2 * (LOG10_2 + RATIO_B_SUR_A * LOG10_3)))
    b_centre = int(RATIO_B_SUR_A * a_centre)
    a_max    = int(DIGITS_CIBLE / (2 * LOG10_2)) + 10

    if DELTA_DEPART > 0:
        delta     = DELTA_DEPART
        side      = 0
        n_filtre7 = n_crible = n_mr = n_paires = 0
        t_cumul   = TEMPS_CUMUL_H * 3600
        print(f"\n  DEPART MANUEL : delta={delta} → a={a_centre + delta}")
    else:
        cp = charger_checkpoint()
        if cp:
            delta     = cp["delta"]
            side      = cp["side"]
            n_filtre7 = cp.get("n_filtre7", 0)
            n_crible  = cp["n_crible"]
            n_mr      = cp["n_mr"]
            n_paires  = cp["n_paires"]
            t_cumul   = cp["t_cumul"]
            print(f"\n  REPRISE depuis le checkpoint :")
            print(f"     delta={delta}  MR={n_mr}  cumul={t_cumul/3600:.2f}h")
        else:
            delta = side = 0
            n_filtre7 = n_crible = n_mr = n_paires = 0
            t_cumul   = 0.0
            print(f"\n  Nouvelle recherche (aucun checkpoint)")

    print(f"\n{'='*65}")
    print(f"  PrimeQuest v5 — 50 000 chiffres  [OPTIMISE — MR×3+20, crible 10M]")
    print(f"{'='*65}")
    print(f"  Cible        : {DIGITS_CIBLE:,} chiffres (±{TOLERANCE})")
    print(f"  a_centre     : {a_centre:,}   b_centre : {b_centre:,}")
    print(f"  Workers      : {N_WORKERS}  (batch {BATCH} paires en parallele)")
    print(f"  Crible       : {len(PREMIERS):,} premiers q≡1(mod3) ≤ {SIEVE_LIMIT:,}")
    print(f"  MR tours     : {MR_TOURS} (init) + {MR_TOURS_CONFIRM} (confirm)")
    print(f"  Timeout      : {TIMEOUT_S:,}s ({TIMEOUT_S/3600:.1f}h)")
    print(f"{'─'*65}\n")

    trouve    = None
    t_debut   = time.perf_counter()
    t_rapport = t_debut
    t_ckpt    = t_debut

    def _handler(sig, frame):
        print(f"\n  Interruption — sauvegarde…", flush=True)
        sauver_checkpoint(delta, side, n_filtre7, n_crible, n_mr, n_paires,
                          t_cumul + (time.perf_counter() - t_debut))
        print(f"  Checkpoint → {CHECKPOINT_FILE}")
        sys.exit(0)

    signal.signal(signal.SIGINT,  _handler)
    signal.signal(signal.SIGTERM, _handler)

    with mp.Pool(N_WORKERS,
                 initializer=_init_worker,
                 initargs=(PREMIERS,)) as pool:

        while delta is not None:

            now           = time.perf_counter()
            elapsed_s     = now - t_debut
            elapsed_total = t_cumul + elapsed_s

            if elapsed_s >= TIMEOUT_S:
                print(f"\n  Timeout {TIMEOUT_S/3600:.1f}h atteint.")
                sauver_checkpoint(delta, side, n_filtre7, n_crible, n_mr,
                                  n_paires, elapsed_total)
                break

            if now - t_ckpt >= CHECKPOINT_S:
                sauver_checkpoint(delta, side, n_filtre7, n_crible, n_mr,
                                  n_paires, elapsed_total)
                t_ckpt = now

            if now - t_rapport >= RAPPORT_S:
                restant = max(0, TIMEOUT_S - elapsed_s)
                taux    = elapsed_s / max(n_mr, 1)
                print(
                    f"  [{elapsed_total/3600:5.2f}h]  "
                    f"a≈{a_depuis_position(a_centre,delta,side):,}  "
                    f"T2={n_filtre7}  crible={n_crible}  "
                    f"MR={n_mr}  ~{taux:.0f}s/MR  restant={restant/60:.0f}min",
                    flush=True
                )
                t_rapport = now

            batch, next_delta, next_side = generer_batch(
                delta, side, a_centre, 1, a_max,
                DIGITS_CIBLE, TOLERANCE, BATCH
            )
            if not batch:
                break

            resultats = pool.map(_tester_paire, batch)

            for r in resultats:
                n_paires += 1
                if r['statut'] == 'filtre7':
                    n_filtre7 += 1
                elif r['statut'] == 'crible':
                    n_crible += 1
                elif r['statut'] == 'compose':
                    n_mr += 1
                    now2    = time.perf_counter()
                    et      = t_cumul + (now2 - t_debut)
                    restant = max(0, TIMEOUT_S - (now2 - t_debut))
                    print(
                        f"  [{et/3600:5.2f}h]  MR #{n_mr:4d}  "
                        f"a={r['a']}, b={r['b']}  "
                        f"(restant {restant/60:.0f}min)  → compose",
                        flush=True
                    )
                elif r['statut'] == 'probable':
                    n_mr += 1
                    now2    = time.perf_counter()
                    et      = t_cumul + (now2 - t_debut)
                    restant = max(0, TIMEOUT_S - (now2 - t_debut))
                    print(
                        f"\n  [{et/3600:5.2f}h]  MR #{n_mr}  "
                        f"a={r['a']}, b={r['b']}  (restant {restant/60:.0f}min)",
                        flush=True
                    )
                    print("  *** PROBABLE — Pocklington… ***", flush=True)
                    ok, temoins = pocklington(r['p'], r['a'], r['b'])
                    if ok:
                        trouve = (r['a'], r['b'], r['m'], r['p'],
                                  int(gmpy2.num_digits(r['p'])), temoins)
                        break
                    else:
                        manque = [q for q, v in temoins.items() if v is None]
                        print(f"  Pocklington incomplet — q={manque}")

            if trouve:
                break

            delta, side = next_delta, next_side

    t_session = time.perf_counter() - t_debut
    t_total   = t_cumul + t_session

    print(f"\n{'='*65}")

    if trouve:
        supprimer_checkpoint()
        a, b, m, p, nb, temoins = trouve
        deux  = gmpy2.mpz(2)
        trois = gmpy2.mpz(3)
        F  = deux**a * trois**(b + 1)
        sp = str(p)
        sm = str(m)
        total_test = n_filtre7 + n_crible + n_mr

        print(f"  PREMIER DE {nb} CHIFFRES — PRIMALITE PROUVEE")
        print(f"{'='*65}")
        print(f"""
  a = {a},  b = {b}  (ratio b/a = {b/a:.4f})
  m = 2^{a}·3^{b}-1  ({gmpy2.num_digits(m)} chiffres)
  p = 3m(m+1)+1  ({nb} chiffres)

  F = 2^{a}·3^{{{b+1}}}  ({gmpy2.num_digits(F)} chiffres) > sqrt(p)
  q=2 → w={temoins[2]}   q=3 → w={temoins[3]}

  Temps total  : {t_total:.1f}s  ({t_total/3600:.2f}h)
  MR tests     : {n_mr}
  Filtres T2   : {n_filtre7}  ({n_filtre7/max(total_test,1)*100:.1f}%)
  Crible       : {n_crible}   ({n_crible/max(total_test,1)*100:.1f}%)
  Workers      : {N_WORKERS}

  p = {sp[:40]}...{sp[-20:]}
""")
        nom = f"premier_{nb}chiffres_v5.txt"
        with open(nom, "w", encoding='utf-8') as f:
            f.write(f"Premier de {nb} chiffres — PrimeQuest v5\n")
            f.write(f"a={a}, b={b}  (ratio b/a = {b/a:.4f})\n")
            f.write(f"Temoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\n")
            f.write(f"Temps   : {t_total:.1f}s  ({t_total/3600:.2f}h)\n")
            f.write(f"Workers : {N_WORKERS}\n")
            f.write(f"MR tests: {n_mr}  | T2: {n_filtre7} | Crible: {n_crible}\n")
            f.write(f"\nm =\n{sm}\n\np =\n{sp}\n")
        print(f"  Sauvegarde → {nom}")

    elif delta is None:
        supprimer_checkpoint()
        print("  Espace explore — aucun premier trouve. Augmenter TOLERANCE.")

    else:
        tps_mr = t_session / max(n_mr, 1)
        total_test = n_filtre7 + n_crible + n_mr
        print(f"""  COMPTE RENDU (relancer pour continuer)
  Temps session     : {t_session/3600:.2f}h  |  Cumule : {t_total/3600:.2f}h
  Filtres T2        : {n_filtre7:,}  ({n_filtre7/max(total_test,1)*100:.1f}%)
  Elimines (crible) : {n_crible:,}   ({n_crible/max(total_test,1)*100:.1f}%)
  Tests MR          : {n_mr}
  Temps/MR          : {tps_mr:.1f}s  (divise par {N_WORKERS} workers)
  Checkpoint        : {CHECKPOINT_FILE}
""")
