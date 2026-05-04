#!/usr/bin/env python3
"""
PrimeQuest v4 — Recherche et preuve de primalité MULTI-CŒURS
Famille : p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Preuve  : Théorème de Pocklington N-1

Optimisations (en plus de v3) :
  4. Filtre mod-7 (Théorème 2) — élimine 1/3 des paires au coût O(1)
     7 | p  ⟺  (a%3, b%6) ∈ F7 = {(0,2),(0,3),(1,0),(1,1),(2,4),(2,5)}

Ratio b/a :
  Résultat 29998 ch. (a=19435, b=19173) : ratio b/a = 0.9865 ≈ 1.0
  Résultat 19999 ch. (a=12228, b=13242) : ratio b/a = 1.083
  → ratio=1.0 donne un centre plus robuste (plus proche en moyenne)
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
# PARAMÈTRES — modifier ici
# ═══════════════════════════════════════════════════════════════
DIGITS_CIBLE    = 30_000      # nombre de chiffres souhaité
TIMEOUT_S       = 14_400      # durée max par session (14400 = 4h)
TOLERANCE       = 10          # tolérance ±chiffres
MR_TOURS        = 25          # tours Miller-Rabin
TEMOINS_MAX     = 300         # max témoins Pocklington
SIEVE_LIMIT     = 1_000_000   # crible jusqu'à ce nombre
RAPPORT_S       = 300         # rapport toutes les 5 min
CHECKPOINT_S    = 60          # sauvegarde checkpoint toutes les 60s
CHECKPOINT_FILE = f"checkpoint_{DIGITS_CIBLE}.json"

# Ratio b/a pour le centre du zigzag
RATIO_B_SUR_A   = 1.000       # 1.0 = centre symétrique (le plus robuste)

# ───────────────────────────────────────────────────────────────
# PARALLÉLISME — ajuster si besoin
# ───────────────────────────────────────────────────────────────
N_WORKERS_PARAM = 0           # 0 = automatique (nb_coeurs - 1)

# ───────────────────────────────────────────────────────────────
# DÉPART MANUEL — utile si le checkpoint n'existe pas
# ───────────────────────────────────────────────────────────────
DELTA_DEPART    = 0           # 0 = automatique
TEMPS_CUMUL_H   = 0.0         # heures déjà passées (sessions précédentes)
# ═══════════════════════════════════════════════════════════════

LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)

# Filtre mod-7 (Théorème 2) — paires interdites (a%3, b%6)
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
    Teste une paire (a, b) : filtre mod-7, crible, Miller-Rabin.
    Retourne {'statut': 'filtre7'|'crible'|'compose'|'probable', ...}.
    """
    a, b, mr_tours = args
    # Étape 1 : filtre mod-7 (Théorème 2) — O(1)
    if (a % 3, b % 6) in F7_FORBIDDEN:
        return {'statut': 'filtre7', 'a': a, 'b': b}
    # Étape 2 : crible (q ≡ 1 mod 3)
    if not _crible_mod(a, b):
        return {'statut': 'crible', 'a': a, 'b': b}
    # Étape 3 : Miller-Rabin
    deux  = gmpy2.mpz(2)
    trois = gmpy2.mpz(3)
    m = deux**a * trois**b - 1
    p = 3 * m * (m + 1) + 1
    if gmpy2.is_prime(p, mr_tours):
        return {'statut': 'probable', 'a': a, 'b': b, 'p': p, 'm': m}
    return {'statut': 'compose', 'a': a, 'b': b}


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
            paires.append((a, b, MR_TOURS))
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
            print(f"  Checkpoint ignoré (digits={data.get('digits')} != {DIGITS_CIBLE})")
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

    # Centre du zigzag avec ratio b/a
    a_centre = int((DIGITS_CIBLE - LOG10_3) / (2 * (LOG10_2 + RATIO_B_SUR_A * LOG10_3)))
    b_centre = int(RATIO_B_SUR_A * a_centre)
    a_max    = int(DIGITS_CIBLE / (2 * LOG10_2)) + 10

    if DELTA_DEPART > 0:
        delta    = DELTA_DEPART
        side     = 0
        n_filtre7 = n_crible = n_mr = n_paires = 0
        t_cumul  = TEMPS_CUMUL_H * 3600
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
            print(f"     delta={delta}, side={side}  "
                  f"→  a={a_depuis_position(a_centre, delta, side)}")
            print(f"     Temps cumule : {t_cumul/3600:.2f}h  |  MR : {n_mr}")
        else:
            delta = side = 0
            n_filtre7 = n_crible = n_mr = n_paires = 0
            t_cumul  = 0.0
            print(f"\n  Nouvelle recherche (aucun checkpoint)")

    print(f"\n{'='*65}")
    print(f"  PrimeQuest v4 — p = 3m(m+1)+1,  m = 2^a·3^b-1  [PARALLELE]")
    print(f"{'='*65}")
    print(f"  Cible        : {DIGITS_CIBLE:,} chiffres (±{TOLERANCE})")
    print(f"  a_centre     : {a_centre:,}   b_centre : {b_centre:,}   (ratio={RATIO_B_SUR_A})")
    print(f"  Workers      : {N_WORKERS}  (batch {BATCH} paires en parallele)")
    print(f"  Crible       : {len(PREMIERS):,} premiers q≡1(mod3) ≤ {SIEVE_LIMIT:,}")
    print(f"  Filtre mod-7 : F7 = {sorted(F7_FORBIDDEN)}  (1/3 elimines)")
    print(f"  Timeout      : {TIMEOUT_S:,}s ({TIMEOUT_S/3600:.1f}h) par session")
    print(f"  Checkpoint   : toutes les {CHECKPOINT_S}s → {CHECKPOINT_FILE}")
    print(f"{'─'*65}\n")

    trouve    = None
    t_debut   = time.perf_counter()
    t_rapport = t_debut
    t_ckpt    = t_debut

    def _handler(sig, frame):
        print(f"\n  Interruption — sauvegarde du checkpoint…", flush=True)
        sauver_checkpoint(delta, side, n_filtre7, n_crible, n_mr, n_paires,
                          t_cumul + (time.perf_counter() - t_debut))
        print(f"  Checkpoint sauvegarde → {CHECKPOINT_FILE}")
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
                print(f"  Checkpoint → {CHECKPOINT_FILE}  (relancer pour continuer)")
                break

            if now - t_ckpt >= CHECKPOINT_S:
                sauver_checkpoint(delta, side, n_filtre7, n_crible, n_mr,
                                  n_paires, elapsed_total)
                t_ckpt = now

            if now - t_rapport >= RAPPORT_S:
                restant = max(0, TIMEOUT_S - elapsed_s)
                taux    = elapsed_s / max(n_mr, 1)
                total_elim = n_filtre7 + n_crible
                print(
                    f"  [{elapsed_total/3600:5.2f}h]  "
                    f"a≈{a_depuis_position(a_centre,delta,side):,}  "
                    f"paires={n_paires:,}  T2={n_filtre7}  crible={n_crible}  "
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
                        f"a={r['a']}, b={r['b']}  "
                        f"(restant {restant/60:.0f}min)",
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
        total_elim = n_filtre7 + n_crible
        total_test = n_filtre7 + n_crible + n_mr

        print(f"  PREMIER DE {nb} CHIFFRES — PRIMALITE PROUVEE")
        print(f"{'='*65}")
        print(f"""
  Structure
  ─────────────────────────────────────────────────
  a = {a},  b = {b}  (ratio b/a = {b/a:.4f})
  m = 2^{a}·3^{b}-1         ({gmpy2.num_digits(m)} chiffres)
  p = 3m(m+1)+1            ({nb} chiffres)

  Preuve Pocklington
  ─────────────────────────────────────────────────
  p-1 = 2^{a}·3^{{{b+1}}}·m
  F   = 2^{a}·3^{{{b+1}}}  ({gmpy2.num_digits(F)} chiffres) > sqrt(p)
  q=2  →  temoin w = {temoins[2]}
  q=3  →  temoin w = {temoins[3]}

  Performance
  ─────────────────────────────────────────────────
  Temps session      : {t_session:.1f}s  ({t_session/60:.1f} min)
  Temps total cumule : {t_total:.1f}s  ({t_total/3600:.2f}h)
  Tests MR           : {n_mr}
  Filtres mod-7 (T2) : {n_filtre7}  ({n_filtre7/max(total_test,1)*100:.1f}%)
  Elimines (crible)  : {n_crible}   ({n_crible/max(total_test,1)*100:.1f}%)
  Total elimines     : {total_elim}  ({total_elim/max(total_test,1)*100:.1f}%)
  Workers paralleles : {N_WORKERS}

  p = {sp[:40]}…{sp[-20:]}
""")
        nom = f"premier_{nb}chiffres_v4.txt"
        with open(nom, "w", encoding='utf-8') as f:
            f.write(f"Premier de {nb} chiffres — PrimeQuest v4\n")
            f.write(f"a={a}, b={b}  (ratio b/a = {b/a:.4f})\n")
            f.write(f"Temoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\n")
            f.write(f"Temps   : {t_total:.1f}s  ({t_total/3600:.2f}h)\n")
            f.write(f"Workers : {N_WORKERS}\n")
            f.write(f"MR tests: {n_mr}  | Filtres T2: {n_filtre7} | Crible: {n_crible}\n")
            f.write(f"\nm =\n{sm}\n\np =\n{sp}\n")
        print(f"  Sauvegarde → {nom}")

    elif delta is None:
        supprimer_checkpoint()
        print("  Espace (a,b) entierement explore — aucun premier trouve.")
        print("  Augmenter TOLERANCE.")

    else:
        tps_mr = t_session / max(n_mr, 1)
        total_test = n_filtre7 + n_crible + n_mr
        print(f"  COMPTE RENDU — Session terminee (relancer pour continuer)")
        print(f"{'='*65}")
        print(f"""
  Temps session           : {t_session/3600:.2f}h
  Temps total cumule      : {t_total/3600:.2f}h
  a explore               : a_centre={a_centre} +/- {delta}
  Paires testees          : {n_paires:,}
  Filtres mod-7 (T2)      : {n_filtre7:,}  ({n_filtre7/max(total_test,1)*100:.1f}%)
  Eliminees (crible)      : {n_crible:,}   ({n_crible/max(total_test,1)*100:.1f}%)
  Tests MR                : {n_mr}
  Temps/MR effectif       : {tps_mr:.1f}s  (divise par {N_WORKERS} workers)
  Checkpoint              : {CHECKPOINT_FILE}
""")
