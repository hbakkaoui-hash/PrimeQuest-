#!/usr/bin/env python3
"""
PrimeQuest v6 — Exploration Multi-Ratio b/a
Famille : p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Preuve  : Théorème de Pocklington N-1

Optimisations cumulées (v1 → v6) :
  1. Crible restreint aux q ≡ 1 (mod 3) — Théorème 3 (~50% candidats en moins)
  2. Filtrage Théorème 2 — 7|p détecté sans calcul si (a%3,b%6) interdit
     → élimine ~33% des paires supplémentaires, coût nul
  3. Parallélisation multi-cœurs via multiprocessing (auto-détection)
  4. Checkpoint toutes les 60s — reprise exacte après interruption
  5. MR_TOURS = 3 au lieu de 25 (facteur 8× sur le coût dominant)
     Justification : Miller-Rabin est DÉTERMINISTE pour les premiers —
     un premier passe TOUJOURS. Les faux positifs (~1.6%) sont
     filtrés par MR_TOURS_CONFIRM + Pocklington.
  6. SIEVE_LIMIT = 10 000 000 — 664 579 primes q≡1(mod3), ~15% MR en moins.
  [v6] 7. RATIOS_LIST : exploration simultanée de plusieurs ratios b/a
     Chaque ratio définit un centre (a₀, b₀) distinct dans l'espace (a,b).
     Lors de chaque batch, les paires sont générées en round-robin sur
     tous les centres actifs.
     → Couvre len(RATIOS_LIST)× plus de territoire par session.
     → Augmente la probabilité de tomber sur un premier sans changer
       la densité théorique (le TNP ne dépend pas du ratio).
     Ratios par défaut : [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15]

Gain global v6 vs v5 : facteur ~len(RATIOS_LIST) sur la couverture spatiale.

Usage : python primequest_v6.py
        Modifier DIGITS_CIBLE, RATIOS_LIST, TIMEOUT_S ci-dessous.
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
DIGITS_CIBLE     = 50_000      # nombre de chiffres souhaité
TIMEOUT_S        = 14_400      # durée max par session (14400 = 4h)
TOLERANCE        = 10          # tolérance ±chiffres
MR_TOURS         = 3           # tours Miller-Rabin rapide
MR_TOURS_CONFIRM = 20          # confirmation avant Pocklington
TEMOINS_MAX      = 300         # max témoins Pocklington
SIEVE_LIMIT      = 10_000_000  # crible jusqu'à ce nombre
RAPPORT_S        = 300         # rapport toutes les 5 min
CHECKPOINT_S     = 60          # sauvegarde checkpoint toutes les 60s
CHECKPOINT_FILE  = f"checkpoint_{DIGITS_CIBLE}_v6.json"

# ───────────────────────────────────────────────────────────────
# [V6] RATIOS b/a explorés simultanément
# ───────────────────────────────────────────────────────────────
# Chaque ratio définit un centre indépendant dans l'espace (a, b).
# Les centres sont parcourus en round-robin : on prend une paire
# de chaque ratio à tour de rôle, couvrant tous les territoires
# en parallèle.
#
# Résultats observés (référence) :
#   9 998 ch.  (a=6212,  b=6738)  : b/a = 1.085
#   10 000 ch. (a=6213,  b=6740)  : b/a = 1.085
#   19 999 ch. (a=12228, b=13242) : b/a = 1.083
#   29 998 ch. (a=19435, b=19173) : b/a = 0.987
RATIOS_LIST = [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15]

# ───────────────────────────────────────────────────────────────
# PARALLÉLISME
# ───────────────────────────────────────────────────────────────
N_WORKERS_PARAM = 0            # 0 = auto (cpu_count - 1), >0 = valeur fixe

# ═══════════════════════════════════════════════════════════════

LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)

# ── Théorème 2 — classes interdites mod 7 ────────────────────────────────────
FORBIDDEN_T2 = frozenset({(0,2),(0,3),(1,0),(1,1),(2,4),(2,5)})


# ── Crible d'Ératosthène ─────────────────────────────────────────────────────

def _eratosthene_filtre(lim):
    t = bytearray([1]) * (lim + 1)
    t[0] = t[1] = 0
    for i in range(2, int(lim**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    return [i for i in range(5, lim + 1) if t[i] and i % 3 == 1]


# ── Fonctions worker ──────────────────────────────────────────────────────────

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
    """
    Teste une paire (a, b) :
      1. Théorème 2 : 7|p ?  (gratuit)
      2. Crible modulaire q≡1(mod3)
      3. Miller-Rabin rapide (MR_TOURS rounds)
      4. Confirmation MR (MR_TOURS_CONFIRM rounds) avant Pocklington
    """
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


# ── Zigzag ────────────────────────────────────────────────────────────────────

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


# ── [V6] Génération multi-ratio en round-robin ────────────────────────────────

def generer_batch_multi(deltas, sides, centres, a_max, digits, tol, taille):
    """
    Génère `taille` paires en round-robin sur tous les ratios actifs.
    Pour chaque cycle : on prend une valeur de a de chaque ratio à tour de rôle,
    puis tous les b valides pour cet a.

    Retourne (paires, new_deltas, new_sides).
    """
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
            nd, ns         = position_suivante(new_deltas[i], new_sides[i], ac, 1, a_max)
            new_deltas[i]  = nd
            new_sides[i]   = ns
            if len(paires) >= taille:
                break
        if not progress:
            break  # tous les ratios épuisés

    return paires, new_deltas, new_sides


# ── Checkpoint ────────────────────────────────────────────────────────────────

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
        if data.get("digits") != DIGITS_CIBLE:
            print(f"  ⚠  Checkpoint ignoré (digits={data.get('digits')} ≠ {DIGITS_CIBLE})")
            return None
        if data.get("version") != 6:
            print(f"  ⚠  Checkpoint ignoré (version {data.get('version')} ≠ 6)")
            return None
        if data.get("ratios") != RATIOS_LIST:
            print(f"  ⚠  Checkpoint ignoré (ratios différents — nouvelle recherche)")
            return None
        return data
    except Exception as e:
        print(f"  ⚠  Checkpoint illisible : {e}")
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

    # Crible
    print(f"Construction du crible (q ≡ 1 mod 3) jusqu'à {SIEVE_LIMIT:,}…",
          end=" ", flush=True)
    t0 = time.perf_counter()
    PREMIERS = _eratosthene_filtre(SIEVE_LIMIT)
    print(f"{len(PREMIERS):,} premiers utiles  ({time.perf_counter()-t0:.2f}s)")
    print(f"Cœurs disponibles : {n_coeurs}  →  {N_WORKERS} workers  (batch {BATCH})")

    # Centres — un par ratio
    a_max   = int(DIGITS_CIBLE / (2 * LOG10_2)) + 10
    centres = []
    for r in RATIOS_LIST:
        ac = int((DIGITS_CIBLE / 2) / (LOG10_2 + r * LOG10_3))
        bc = int(ac * r)
        centres.append((ac, bc))

    N_RATIOS = len(RATIOS_LIST)

    # État initial
    cp = charger_checkpoint()
    if cp:
        raw_deltas = cp["deltas"]
        deltas = [None if d == -1 else d for d in raw_deltas]
        sides  = cp["sides"]
        n_t2     = cp.get("n_t2", 0)
        n_crible = cp["n_crible"]
        n_mr     = cp["n_mr"]
        n_paires = cp["n_paires"]
        t_cumul  = cp["t_cumul"]
        print(f"\n  ♻  REPRISE v6 : {sum(1 for d in deltas if d is not None)}/{N_RATIOS} ratios actifs")
        print(f"     Cumulé : {t_cumul/3600:.2f}h  |  MR : {n_mr}")
    else:
        deltas   = [0] * N_RATIOS
        sides    = [0] * N_RATIOS
        n_t2 = n_crible = n_mr = n_paires = 0
        t_cumul = 0.0
        print(f"\n  Nouvelle recherche v6 (aucun checkpoint)")

    print(f"\n{'═'*70}")
    print(f"  PrimeQuest v6 — p = 3m(m+1)+1,  m = 2^a·3^b−1")
    print(f"{'═'*70}")
    print(f"  Cible          : {DIGITS_CIBLE:,} chiffres (±{TOLERANCE})")
    print(f"  Ratios b/a     : {RATIOS_LIST}")
    print(f"  Centres (a₀)   : {[c[0] for c in centres]}")
    print(f"  a_max          : {a_max:,}")
    print(f"  Workers        : {N_WORKERS}  (batch {BATCH} paires)")
    print(f"  Crible         : {len(PREMIERS):,} premiers q≡1(mod3) ≤ {SIEVE_LIMIT:,}")
    print(f"  Théorème 2     : {len(FORBIDDEN_T2)}/18 classes (a%3,b%6) filtrées (~33%)")
    print(f"  MR rounds      : {MR_TOURS} (filtre) + {MR_TOURS_CONFIRM} (confirmation)")
    print(f"  Timeout        : {TIMEOUT_S:,}s ({TIMEOUT_S/3600:.1f}h)")
    print(f"  Checkpoint     : toutes les {CHECKPOINT_S}s → {CHECKPOINT_FILE}")
    print(f"{'─'*70}\n")

    trouve    = None
    t_debut   = time.perf_counter()
    t_rapport = t_debut
    t_ckpt    = t_debut

    def _handler(sig, frame):
        print(f"\n  ⚡ Interruption — sauvegarde…", flush=True)
        sauver_checkpoint(deltas, sides, n_t2, n_crible, n_mr, n_paires,
                          t_cumul + (time.perf_counter() - t_debut))
        print(f"  Checkpoint → {CHECKPOINT_FILE}")
        sys.exit(0)

    signal.signal(signal.SIGINT,  _handler)
    signal.signal(signal.SIGTERM, _handler)

    with mp.Pool(N_WORKERS,
                 initializer=_init_worker,
                 initargs=(PREMIERS, FORBIDDEN_T2)) as pool:

        while any(d is not None for d in deltas):

            now           = time.perf_counter()
            elapsed_s     = now - t_debut
            elapsed_total = t_cumul + elapsed_s

            if elapsed_s >= TIMEOUT_S:
                print(f"\n  ⏱  Timeout {TIMEOUT_S/3600:.1f}h atteint.")
                sauver_checkpoint(deltas, sides, n_t2, n_crible, n_mr, n_paires,
                                  elapsed_total)
                print(f"  Checkpoint → {CHECKPOINT_FILE}")
                break

            if now - t_ckpt >= CHECKPOINT_S:
                sauver_checkpoint(deltas, sides, n_t2, n_crible, n_mr, n_paires,
                                  elapsed_total)
                t_ckpt = now

            if now - t_rapport >= RAPPORT_S:
                restant   = max(0, TIMEOUT_S - elapsed_s)
                taux      = elapsed_s / max(n_mr, 1)
                elim_tot  = n_t2 + n_crible
                actifs    = sum(1 for d in deltas if d is not None)
                print(
                    f"  [{elapsed_total/3600:5.2f}h]  "
                    f"ratios_actifs={actifs}/{N_RATIOS}  "
                    f"paires={n_paires:,}  T2={n_t2:,}  crible={n_crible:,}  "
                    f"MR={n_mr}  élim={elim_tot/max(n_paires,1)*100:.1f}%  "
                    f"~{taux:.0f}s/MR  restant={restant/60:.0f}min",
                    flush=True
                )
                t_rapport = now

            batch, deltas, sides = generer_batch_multi(
                deltas, sides, centres, a_max,
                DIGITS_CIBLE, TOLERANCE, BATCH
            )
            if not batch:
                break

            resultats = pool.map(_tester_paire, batch)

            for r in resultats:
                n_paires += 1
                if r['statut'] == 'theoreme2':
                    n_t2 += 1

                elif r['statut'] == 'crible':
                    n_crible += 1

                elif r['statut'] == 'compose':
                    n_mr += 1
                    now2    = time.perf_counter()
                    et      = t_cumul + (now2 - t_debut)
                    restant = max(0, TIMEOUT_S - (now2 - t_debut))
                    ratio_r = r['b'] / r['a']
                    print(
                        f"  [{et/3600:5.2f}h]  MR #{n_mr:4d}  "
                        f"a={r['a']}, b={r['b']}  r={ratio_r:.3f}  "
                        f"(restant {restant/60:.0f}min)  → composé",
                        flush=True
                    )

                elif r['statut'] == 'probable':
                    n_mr += 1
                    now2    = time.perf_counter()
                    et      = t_cumul + (now2 - t_debut)
                    restant = max(0, TIMEOUT_S - (now2 - t_debut))
                    ratio_r = r['b'] / r['a']
                    print(
                        f"\n  [{et/3600:5.2f}h]  MR #{n_mr}  "
                        f"a={r['a']}, b={r['b']}  r={ratio_r:.3f}  "
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
                        print(f"  ⚠  Pocklington incomplet — q={manque}")

            if trouve:
                break

    t_session = time.perf_counter() - t_debut
    t_total   = t_cumul + t_session

    print(f"\n{'═'*70}")

    if trouve:
        supprimer_checkpoint()
        a, b, m, p, nb, temoins = trouve
        deux  = gmpy2.mpz(2)
        trois = gmpy2.mpz(3)
        F  = deux**a * trois**(b + 1)
        sp = str(p)
        ratio_r = b / a

        print(f"  PREMIER DE {nb} CHIFFRES — PRIMALITÉ PROUVÉE ✅")
        print(f"{'═'*70}")
        print(f"""
  Structure
  ──────────────────────────────────────────────────────
  a = {a},  b = {b}  (ratio b/a = {ratio_r:.4f})
  m = 2^{a}·3^{b}−1   ({gmpy2.num_digits(m)} chiffres)
  p = 3m(m+1)+1      ({nb} chiffres)

  Preuve Pocklington
  ──────────────────────────────────────────────────────
  F = 2^{a}·3^{{{b+1}}}  ({gmpy2.num_digits(F)} chiffres) > √p  ✅
  q=2  →  témoin w = {temoins[2]}   ✅
  q=3  →  témoin w = {temoins[3]}   ✅

  Performance
  ──────────────────────────────────────────────────────
  Temps session      : {t_session:.1f}s  ({t_session/3600:.2f}h)
  Temps total cumulé : {t_total:.1f}s  ({t_total/3600:.2f}h)
  Tests MR           : {n_mr}
  Filtrés Théorème 2 : {n_t2}  ({n_t2/max(n_paires,1)*100:.1f}%)
  Éliminés crible    : {n_crible}  ({n_crible/max(n_paires,1)*100:.1f}%)
  Total éliminés     : {n_t2+n_crible}  ({(n_t2+n_crible)/max(n_paires,1)*100:.1f}%)
  Workers parallèles : {N_WORKERS}
  Ratios explorés    : {RATIOS_LIST}

  p = {sp[:40]}…{sp[-20:]}
""")
        nom = f"premier_{nb}chiffres_v6.txt"
        with open(nom, "w", encoding="utf-8") as f:
            f.write(f"Premier de {nb} chiffres — PrimeQuest v6\n")
            f.write(f"a={a}, b={b}  (ratio b/a={ratio_r:.4f})\n")
            f.write(f"Témoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\n")
            f.write(f"Temps   : {t_total:.1f}s  ({t_total/3600:.2f}h)\n")
            f.write(f"Workers : {N_WORKERS}\n")
            f.write(f"Ratios  : {RATIOS_LIST}\n\nm =\n{m}\n\np =\n{p}\n")
        print(f"  Sauvegardé → {nom}")

    elif not any(d is not None for d in deltas):
        supprimer_checkpoint()
        print("  Espace (a,b) entièrement exploré — augmenter TOLERANCE.")

    else:
        elim_tot = n_t2 + n_crible
        tps_mr   = t_session / max(n_mr, 1)
        actifs   = sum(1 for d in deltas if d is not None)
        print(f"  COMPTE RENDU — Session terminée (relancer pour continuer)")
        print(f"{'═'*70}")
        print(f"""
  Temps session           : {t_session/3600:.2f}h
  Temps total cumulé      : {t_total/3600:.2f}h
  Ratios actifs           : {actifs}/{N_RATIOS}
  Paires testées          : {n_paires:,}
  Filtrées Théorème 2     : {n_t2:,}  ({n_t2/max(n_paires,1)*100:.1f}%)
  Éliminées crible        : {n_crible:,}  ({n_crible/max(n_paires,1)*100:.1f}%)
  Total éliminées         : {elim_tot:,}  ({elim_tot/max(n_paires,1)*100:.1f}%)
  Tests MR                : {n_mr}
  Temps/MR effectif       : {tps_mr:.1f}s  (÷ {N_WORKERS} workers)
  Checkpoint              : {CHECKPOINT_FILE} ✅
""")
