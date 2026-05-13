#!/usr/bin/env python3
"""
PrimeQuest v7 — Spread Adaptatif + Tolérance élargie
Famille : p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Preuve  : Théorème de Pocklington N-1

Optimisations cumulées (v1 → v7) :
  1. Crible restreint aux q ≡ 1 (mod 3) — Théorème 3 (~50% candidats en moins)
  2. Filtrage Théorème 2 — 7|p détecté sans calcul si (a%3,b%6) interdit
     → élimine ~33% des paires supplémentaires, coût nul
  3. Parallélisation multi-cœurs via multiprocessing (auto-détection)
  4. Checkpoint toutes les 60s — reprise exacte après interruption
  5. MR_TOURS = 3 au lieu de 25 (facteur 8× sur le coût dominant)
  6. SIEVE_LIMIT = 10 000 000 — 664 579 primes q≡1(mod3)
  [v6] 7. Multi-ratio round-robin : plusieurs centres (a₀,b₀) explorés
     simultanément, bornes r_min/r_max calculées depuis DIGITS_CIBLE.
  [v7] 8. TOLERANCE = 50 : accepte p dans [D−50, D+50] chiffres
     → b_valides() retourne ~5× plus de candidats par valeur de a
     → densité de premiers ~5× plus élevée sans surcoût arithmétique.
  [v7] 9. Spread adaptatif : après chaque ADAPT_INTERVAL tests MR sans
     succès, N_RATIOS_ACTIFS augmente de 2 (un ratio de chaque côté).
     Démarre sur N_RATIOS_INIT ratios (plage ±SPREAD_INIT),
     s'élargit automatiquement jusqu'à N_RATIOS_MAX (plage ±SPREAD_MAX).
     Effet : la recherche reste focalisée puis s'élargit si nécessaire.

Gain v7 vs v6 :
  - TOLERANCE 10 → 50 : ~5× plus de candidats b par valeur de a
  - Spread adaptatif : couverture croissante sans intervention manuelle
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
TOLERANCE        = 50          # tolérance ±chiffres  [v7 : 10 → 50]
MR_TOURS         = 3           # tours Miller-Rabin rapide
MR_TOURS_CONFIRM = 20          # confirmation avant Pocklington
TEMOINS_MAX      = 300         # max témoins Pocklington
SIEVE_LIMIT      = 10_000_000  # crible jusqu'à ce nombre
RAPPORT_S        = 300         # rapport toutes les 5 min
CHECKPOINT_S     = 60          # sauvegarde checkpoint toutes les 60s
CHECKPOINT_FILE  = f"checkpoint_{DIGITS_CIBLE}_v7.json"

# ───────────────────────────────────────────────────────────────
# [V7] SPREAD ADAPTATIF
# ───────────────────────────────────────────────────────────────
RATIO_CENTRE    = 1.00   # centre de la plage (empirique)
N_RATIOS_INIT   = 7      # ratios actifs au départ
N_RATIOS_MAX    = 15     # ratios max après expansions complètes
SPREAD_INIT     = 0.15   # demi-amplitude initiale  → [0.85, 1.15]
SPREAD_MAX      = 0.35   # demi-amplitude maximale  → [0.65, 1.35]
ADAPT_INTERVAL  = 200    # tous les N tests MR → +2 ratios

# ───────────────────────────────────────────────────────────────
# PARALLÉLISME
# ───────────────────────────────────────────────────────────────
N_WORKERS_PARAM = 0            # 0 = auto (cpu_count - 1), >0 = valeur fixe

# ═══════════════════════════════════════════════════════════════

LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)

_RATIO_MIN_ABSOLU = 2 * LOG10_2 / DIGITS_CIBLE
_RATIO_MAX_ABSOLU = (DIGITS_CIBLE / 2 - LOG10_2) / LOG10_3

def _generer_tous_ratios():
    step = (2 * SPREAD_MAX / (N_RATIOS_MAX - 1)) if N_RATIOS_MAX > 1 else 0
    bruts = [
        round(RATIO_CENTRE - SPREAD_MAX + i * step, 4)
        for i in range(N_RATIOS_MAX)
    ]
    centre_idx = N_RATIOS_MAX // 2
    ordonnes = []
    for k in range(N_RATIOS_MAX):
        if k % 2 == 0:
            idx = centre_idx - k // 2
        else:
            idx = centre_idx + (k + 1) // 2
        if 0 <= idx < N_RATIOS_MAX:
            ordonnes.append(bruts[idx])
    valides = []
    for r in ordonnes:
        if r < _RATIO_MIN_ABSOLU or r > _RATIO_MAX_ABSOLU:
            continue
        ac = max(1, int((DIGITS_CIBLE / 2) / (LOG10_2 + r * LOG10_3)))
        if int(ac * r) < 1:
            continue
        valides.append(r)
    if len(valides) < N_RATIOS_INIT:
        raise ValueError(
            f"Seulement {len(valides)} ratios valides — "
            f"réduire N_RATIOS_INIT ou ajuster RATIO_CENTRE/SPREAD_MAX."
        )
    return valides

TOUS_RATIOS = _generer_tous_ratios()

def n_actifs_pour(n_mr_total):
    expansions = n_mr_total // ADAPT_INTERVAL
    return min(N_RATIOS_MAX, N_RATIOS_INIT + 2 * expansions)


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
        b for b in range(max(1, int(b_c) - tol - 2), int(b_c) + tol + 3)
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


def generer_batch_multi(deltas, sides, centres, a_max, digits, tol, taille, n_actifs):
    paires     = []
    new_deltas = list(deltas)
    new_sides  = list(sides)

    while len(paires) < taille:
        progress = False
        for i in range(min(n_actifs, len(centres))):
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
            "digits":        DIGITS_CIBLE,
            "version":       7,
            "tous_ratios":   TOUS_RATIOS,
            "deltas":        [d if d is not None else -1 for d in deltas],
            "sides":         sides,
            "n_t2":          n_t2,
            "n_crible":      n_crible,
            "n_mr":          n_mr,
            "n_paires":      n_paires,
            "t_cumul":       t_cumul,
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
        if data.get("version") != 7:
            print(f"  ⚠  Checkpoint ignoré (version {data.get('version')} ≠ 7)")
            return None
        if data.get("tous_ratios") != TOUS_RATIOS:
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

    print(f"Construction du crible (q ≡ 1 mod 3) jusqu'à {SIEVE_LIMIT:,}…",
          end=" ", flush=True)
    t0 = time.perf_counter()
    PREMIERS = _eratosthene_filtre(SIEVE_LIMIT)
    print(f"{len(PREMIERS):,} premiers utiles  ({time.perf_counter()-t0:.2f}s)")
    print(f"Cœurs disponibles : {n_coeurs}  →  {N_WORKERS} workers  (batch {BATCH})")

    a_max   = int(DIGITS_CIBLE / (2 * LOG10_2)) + 10
    centres = []
    for r in TOUS_RATIOS:
        ac = int((DIGITS_CIBLE / 2) / (LOG10_2 + r * LOG10_3))
        bc = int(ac * r)
        centres.append((ac, bc))

    N_TOTAL = len(TOUS_RATIOS)

    cp = charger_checkpoint()
    if cp:
        raw_deltas = cp["deltas"]
        deltas   = [None if d == -1 else d for d in raw_deltas]
        sides    = cp["sides"]
        n_t2     = cp.get("n_t2", 0)
        n_crible = cp["n_crible"]
        n_mr     = cp["n_mr"]
        n_paires = cp["n_paires"]
        t_cumul  = cp["t_cumul"]
        n_act    = n_actifs_pour(n_mr)
        print(f"\n  ♻  REPRISE v7 : {n_act}/{N_TOTAL} ratios actifs  |  "
              f"MR={n_mr}  cumulé={t_cumul/3600:.2f}h")
    else:
        deltas   = [0] * N_TOTAL
        sides    = [0] * N_TOTAL
        n_t2 = n_crible = n_mr = n_paires = 0
        t_cumul = 0.0
        n_act   = N_RATIOS_INIT
        print(f"\n  Nouvelle recherche v7 (aucun checkpoint)")

    print(f"\n{'═'*72}")
    print(f"  PrimeQuest v7 — p = 3m(m+1)+1,  m = 2^a·3^b−1")
    print(f"{'═'*72}")
    print(f"  Cible           : {DIGITS_CIBLE:,} chiffres  (tolérance ±{TOLERANCE})")
    print(f"  Bornes ratios   : [{_RATIO_MIN_ABSOLU:.6f},  {_RATIO_MAX_ABSOLU:.2f}]")
    print(f"  Spread initial  : ±{SPREAD_INIT}  ({N_RATIOS_INIT} ratios)  →"
          f"  [{RATIO_CENTRE-SPREAD_INIT:.3f}, {RATIO_CENTRE+SPREAD_INIT:.3f}]")
    print(f"  Spread max      : ±{SPREAD_MAX}  ({N_RATIOS_MAX} ratios)  →"
          f"  [{RATIO_CENTRE-SPREAD_MAX:.3f}, {RATIO_CENTRE+SPREAD_MAX:.3f}]")
    print(f"  Expansion       : +2 ratios tous les {ADAPT_INTERVAL} tests MR")
    print(f"  Ratios (tous)   : {TOUS_RATIOS}")
    print(f"  Ratios actifs   : {TOUS_RATIOS[:n_act]}  ({n_act}/{N_TOTAL})")
    print(f"  Centres (a₀)    : {[c[0] for c in centres[:n_act]]}")
    print(f"  a_max           : {a_max:,}")
    print(f"  Workers         : {N_WORKERS}  (batch {BATCH} paires)")
    print(f"  Crible          : {len(PREMIERS):,} premiers q≡1(mod3) ≤ {SIEVE_LIMIT:,}")
    print(f"  Théorème 2      : {len(FORBIDDEN_T2)}/18 classes filtrées (~33%)")
    print(f"  MR rounds       : {MR_TOURS} + {MR_TOURS_CONFIRM}")
    print(f"  Timeout         : {TIMEOUT_S:,}s ({TIMEOUT_S/3600:.1f}h)")
    print(f"  Checkpoint      : toutes les {CHECKPOINT_S}s → {CHECKPOINT_FILE}")
    print(f"{'─'*72}\n")

    trouve    = None
    t_debut   = time.perf_counter()
    t_rapport = t_debut
    t_ckpt    = t_debut
    n_act_prec = n_act

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

        while any(d is not None for d in deltas[:n_act]):

            now           = time.perf_counter()
            elapsed_s     = now - t_debut
            elapsed_total = t_cumul + elapsed_s

            n_act = n_actifs_pour(n_mr)
            if n_act > n_act_prec:
                nouveaux = TOUS_RATIOS[n_act_prec:n_act]
                print(f"\n  ↗  EXPANSION spread : {n_act_prec}→{n_act} ratios  "
                      f"(+{nouveaux})  après {n_mr} tests MR", flush=True)
                n_act_prec = n_act

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
                restant  = max(0, TIMEOUT_S - elapsed_s)
                taux     = elapsed_s / max(n_mr, 1)
                elim_tot = n_t2 + n_crible
                spread_c = SPREAD_INIT + (n_mr // ADAPT_INTERVAL) * (
                    (SPREAD_MAX - SPREAD_INIT) / max(1, (N_RATIOS_MAX - N_RATIOS_INIT) // 2)
                )
                print(
                    f"  [{elapsed_total/3600:5.2f}h]  "
                    f"ratios={n_act}/{N_TOTAL}  spread=±{min(spread_c, SPREAD_MAX):.3f}  "
                    f"paires={n_paires:,}  T2={n_t2:,}  crible={n_crible:,}  "
                    f"MR={n_mr}  élim={elim_tot/max(n_paires,1)*100:.1f}%  "
                    f"~{taux:.0f}s/MR  restant={restant/60:.0f}min",
                    flush=True
                )
                t_rapport = now

            batch, deltas, sides = generer_batch_multi(
                deltas, sides, centres, a_max,
                DIGITS_CIBLE, TOLERANCE, BATCH, n_act
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

    print(f"\n{'═'*72}")

    if trouve:
        supprimer_checkpoint()
        a, b, m, p, nb, temoins = trouve
        deux  = gmpy2.mpz(2)
        trois = gmpy2.mpz(3)
        F  = deux**a * trois**(b + 1)
        sp = str(p)
        ratio_r = b / a

        print(f"  PREMIER DE {nb} CHIFFRES — PRIMALITÉ PROUVÉE ✅")
        print(f"{'═'*72}")
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
  Tolérance          : ±{TOLERANCE} chiffres
  Ratios actifs      : {TOUS_RATIOS[:n_act]}

  p = {sp[:40]}…{sp[-20:]}
""")
        nom = f"premier_{nb}chiffres_v7.txt"
        with open(nom, "w", encoding="utf-8") as f:
            f.write(f"Premier de {nb} chiffres — PrimeQuest v7\n")
            f.write(f"a={a}, b={b}  (ratio b/a={ratio_r:.4f})\n")
            f.write(f"Témoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\n")
            f.write(f"Temps   : {t_total:.1f}s  ({t_total/3600:.2f}h)\n")
            f.write(f"Workers : {N_WORKERS}\n")
            f.write(f"Tolérance : ±{TOLERANCE} chiffres\n")
            f.write(f"Ratios  : {TOUS_RATIOS[:n_act]}\n\nm =\n{m}\n\np =\n{p}\n")
        print(f"  Sauvegardé → {nom}")

    elif not any(d is not None for d in deltas[:n_act]):
        supprimer_checkpoint()
        print("  Espace (a,b) entièrement exploré — augmenter TOLERANCE ou SPREAD_MAX.")

    else:
        elim_tot = n_t2 + n_crible
        tps_mr   = t_session / max(n_mr, 1)
        print(f"  COMPTE RENDU — Session terminée (relancer pour continuer)")
        print(f"{'═'*72}")
        print(f"""
  Temps session           : {t_session/3600:.2f}h
  Temps total cumulé      : {t_total/3600:.2f}h
  Ratios actifs           : {n_act}/{N_TOTAL}
  Paires testées          : {n_paires:,}
  Filtrées Théorème 2     : {n_t2:,}  ({n_t2/max(n_paires,1)*100:.1f}%)
  Éliminées crible        : {n_crible:,}  ({n_crible/max(n_paires,1)*100:.1f}%)
  Total éliminées         : {elim_tot:,}  ({elim_tot/max(n_paires,1)*100:.1f}%)
  Tests MR                : {n_mr}
  Temps/MR effectif       : {tps_mr:.1f}s  (÷ {N_WORKERS} workers)
  Tolérance               : ±{TOLERANCE} chiffres
  Checkpoint              : {CHECKPOINT_FILE} ✅
""")
