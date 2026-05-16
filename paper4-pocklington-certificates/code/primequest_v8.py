#!/usr/bin/env python3
"""
PrimeQuest v8 — Famille k=23, certificat Pocklington
═════════════════════════════════════════════════════
Famille : p = 23·m·(m+1)+1,   m = 2^a · 23^b − 1
Preuve  : Théorème de Pocklington N-1

Pourquoi k=23 vs k=3 ?
──────────────────────
  C(f_{23,+1}) = 5.839  vs  C(f_{3,+1}) = 3.361   [HEURISTIQUE, Bateman–Horn]
  Rapport ≈ 1.74×  — densité de premiers ~74% supérieure à taille égale.
  Source : Bakkaoui (2026), Section 5.1.2

Propriétés RIGOUREUSES et INCONDITIONNELLES :
─────────────────────────────────────────────
  1. p−1 = 2^a · 23^{b+1} · m  →  F := 2^a · 23^{b+1} divise p−1  [RIGOUREUX]
     F > √p car F²/p ≥ 23 > 1 (F/√p ≥ √23 ≈ 4.796)               [RIGOUREUX]

  2. Seuls 2 témoins Pocklington requis : q ∈ {2, 23}              [RIGOUREUX]

  3. Ces 13 premiers ne divisent JAMAIS p = 23·m·(m+1)+1 :         [RIGOUREUX]
     {2, 3, 5, 7, 11, 13, 17, 23, 29, 31, 41, 43, 59}
     Preuve : le discriminant 437 = 19·23 est un non-résidu
     quadratique modulo chacun d'eux → le polynôme 23x²+23x+1
     n'a aucune racine mod q → q ∤ p pour tout m.
     → Crible restreint aux premiers q ∉ blindspots où (437/q) = +1.

Gains v8 vs v7 (k=3) :
───────────────────────
  - Densité ×1.74 (heuristique Bateman–Horn)
  - F/√p ≥ √23 ≈ 4.8  (vs √3 ≈ 1.73 pour k=3) — plus robuste
  - 13 blindspots éliminés du crible (vs ~2 pour k=3)
  - Témoins Pocklington : toujours 2 (q=2, q=23)
  - Architecture identique à v7 (checkpoint, spread adaptatif, parallélisme)
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
K                = 23          # constante de la famille p = K·m·(m+1)+1
DIGITS_CIBLE     = 50_000      # nombre de chiffres souhaité
TIMEOUT_S        = 14_400      # durée max par session (4h)
TOLERANCE        = 50          # tolérance ±chiffres
MR_TOURS         = 3           # tours Miller-Rabin rapide
MR_TOURS_CONFIRM = 20          # confirmation avant Pocklington
TEMOINS_MAX      = 300         # max témoins Pocklington
SIEVE_LIMIT      = 10_000_000  # crible jusqu'à ce nombre
RAPPORT_S        = 300         # rapport toutes les 5 min
CHECKPOINT_S     = 60          # checkpoint toutes les 60s
CHECKPOINT_FILE  = f"checkpoint_{DIGITS_CIBLE}_v8.json"

# ── Spread adaptatif ────────────────────────────────────────────
# Ratio optimal b/a = log10(2)/log10(23) ≈ 0.221 (contributions égales)
# Données empiriques à collecter lors des premières exécutions.
RATIO_CENTRE    = 0.221   # b/a optimal pour k=23
N_RATIOS_INIT   = 7
N_RATIOS_MAX    = 15
SPREAD_INIT     = 0.08    # plage initiale [0.141, 0.301]
SPREAD_MAX      = 0.18    # plage maximale [0.041, 0.401]
ADAPT_INTERVAL  = 200     # +2 ratios tous les N tests MR

# ── Parallélisme ────────────────────────────────────────────────
N_WORKERS_PARAM = 0       # 0 = auto (cpu_count − 1)

# ═══════════════════════════════════════════════════════════════

LOG10_2 = math.log10(2)
LOG10_K = math.log10(K)   # log10(23) ≈ 1.36173

# Formule digits : D ≈ 2a·log10(2) + (2b+1)·log10(23)
_RATIO_MIN_ABSOLU = 2 * LOG10_2 / DIGITS_CIBLE
_RATIO_MAX_ABSOLU = (DIGITS_CIBLE / 2 - LOG10_2) / LOG10_K

# Blindspots : premiers qui ne divisent JAMAIS p = 23·m·(m+1)+1
# pour aucun entier m.  [RIGOUREUX, INCONDITIONNEL]
# Discriminant 437 = 19·23 est non-résidu quadratique mod chacun.
BLINDSPOTS = frozenset({2, 3, 5, 7, 11, 13, 17, 23, 29, 31, 41, 43, 59})


# ── Génération TOUS_RATIOS (centre → bords) ──────────────────────────────────

def _generer_tous_ratios():
    step = (2 * SPREAD_MAX / (N_RATIOS_MAX - 1)) if N_RATIOS_MAX > 1 else 0
    bruts = [
        round(RATIO_CENTRE - SPREAD_MAX + i * step, 4)
        for i in range(N_RATIOS_MAX)
    ]
    centre_idx = N_RATIOS_MAX // 2
    ordonnes = []
    for k_step in range(N_RATIOS_MAX):
        if k_step % 2 == 0:
            idx = centre_idx - k_step // 2
        else:
            idx = centre_idx + (k_step + 1) // 2
        if 0 <= idx < N_RATIOS_MAX:
            ordonnes.append(bruts[idx])
    valides = []
    for r in ordonnes:
        if r <= 0 or r < _RATIO_MIN_ABSOLU or r > _RATIO_MAX_ABSOLU:
            continue
        ac = max(1, int((DIGITS_CIBLE - LOG10_K) / (2 * (LOG10_2 + r * LOG10_K))))
        if int(ac * r) < 1:
            continue
        valides.append(r)
    if len(valides) < N_RATIOS_INIT:
        raise ValueError(
            f"Seulement {len(valides)} ratios valides — "
            f"ajuster RATIO_CENTRE/SPREAD_MAX/N_RATIOS_INIT."
        )
    return valides

TOUS_RATIOS = _generer_tous_ratios()

def n_actifs_pour(n_mr_total):
    expansions = n_mr_total // ADAPT_INTERVAL
    return min(N_RATIOS_MAX, N_RATIOS_INIT + 2 * expansions)


# ── Crible — primes où 437 est un carré mod q (peuvent diviser p) ─────────────
# Pour q ∉ BLINDSPOTS : la congruence 23x(x+1)+1 ≡ 0 (mod q) a des solutions
# ssi le discriminant 437 est un résidu quadratique mod q.
# Test : pow(437, (q-1)//2, q) == 1   (critère d'Euler)

def _eratosthene_filtre(lim):
    """Primes q ≤ lim, non-blindspot, où 437 est QR mod q."""
    t = bytearray([1]) * (lim + 1)
    t[0] = t[1] = 0
    for i in range(2, int(lim**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    result = []
    for q in range(3, lim + 1):
        if not t[q]:
            continue
        if q in BLINDSPOTS:
            continue
        # Critère d'Euler : 437 est QR mod q iff 437^{(q-1)/2} ≡ 1 (mod q)
        if pow(437, (q - 1) // 2, q) == 1:
            result.append(q)
    return result


# ── Worker ────────────────────────────────────────────────────────────────────

_W_PREMIERS = None

def _init_worker(premiers):
    global _W_PREMIERS
    _W_PREMIERS = premiers
    signal.signal(signal.SIGINT,  signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_IGN)

def _crible_mod(a, b):
    """Teste si p = 23·m·(m+1)+1 est divisible par un premier du crible."""
    for q in _W_PREMIERS:
        q1 = q - 1
        # m = 2^a · 23^b − 1  (mod q), par petit théorème de Fermat
        x     = (pow(2, a % q1, q) * pow(K, b % q1, q)) % q
        m_mod = (x - 1) % q
        p_mod = (K * m_mod * (m_mod + 1) + 1) % q
        if p_mod == 0:
            return False
    return True

def _tester_paire(args):
    a, b, mr_tours, mr_confirm = args

    if not _crible_mod(a, b):
        return {'statut': 'crible', 'a': a, 'b': b}

    vingt_trois = gmpy2.mpz(K)
    deux        = gmpy2.mpz(2)
    m = deux**a * vingt_trois**b - 1
    p = K * m * (m + 1) + 1

    if not gmpy2.is_prime(p, mr_tours):
        return {'statut': 'compose', 'a': a, 'b': b}

    if not gmpy2.is_prime(p, mr_confirm):
        return {'statut': 'compose', 'a': a, 'b': b}

    return {'statut': 'probable', 'a': a, 'b': b, 'p': p, 'm': m}


# ── b valides pour un a donné ─────────────────────────────────────────────────
# D ≈ 2a·log10(2) + (2b+1)·log10(23)
# b_centre = (D − 2a·log10(2) − log10(23)) / (2·log10(23))

def b_valides(a, digits, tol):
    b_c = (digits - 2 * a * LOG10_2 - LOG10_K) / (2 * LOG10_K)
    return [
        b for b in range(max(1, int(b_c) - tol - 2), int(b_c) + tol + 3)
        if abs(2 * a * LOG10_2 + (2 * b + 1) * LOG10_K - digits) <= tol + 2
    ]


# ── Pocklington — témoins q ∈ {2, 23}  [RIGOUREUX, INCONDITIONNEL] ───────────
# p−1 = 2^a · 23^{b+1} · m
# F   = 2^a · 23^{b+1} > √p  (car F²/p ≥ 23 > 1)
# Pour chaque facteur premier q de F (q=2 et q=23) :
#   trouver w tel que w^{p−1} ≡ 1 (mod p)  ET  pgcd(w^{(p−1)/q}−1, p) = 1

def pocklington(p, a, b):
    p1 = p - 1
    v  = {2: None, 23: None}
    for w in range(2, 2 + TEMOINS_MAX):
        gw = gmpy2.mpz(w)
        if gmpy2.powmod(gw, p1, p) != 1:
            continue
        for q in [2, 23]:
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


# ── Génération multi-ratio round-robin ────────────────────────────────────────

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


# ── Checkpoint ────────────────────────────────────────────────────────────────

def sauver_checkpoint(deltas, sides, n_crible, n_mr, n_paires, t_cumul):
    with open(CHECKPOINT_FILE, "w", encoding="utf-8") as f:
        json.dump({
            "digits":      DIGITS_CIBLE,
            "version":     8,
            "k":           K,
            "tous_ratios": TOUS_RATIOS,
            "deltas":      [d if d is not None else -1 for d in deltas],
            "sides":       sides,
            "n_crible":    n_crible,
            "n_mr":        n_mr,
            "n_paires":    n_paires,
            "t_cumul":     t_cumul,
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
        if data.get("version") != 8:
            print(f"  ⚠  Checkpoint ignoré (version {data.get('version')} ≠ 8)")
            return None
        if data.get("k") != K:
            print(f"  ⚠  Checkpoint ignoré (k={data.get('k')} ≠ {K})")
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

    # Crible — primes où 437 est QR (peuvent diviser p), hors blindspots
    print(f"Construction du crible (437 QR mod q, hors {len(BLINDSPOTS)} blindspots) "
          f"jusqu'à {SIEVE_LIMIT:,}…", end=" ", flush=True)
    t0 = time.perf_counter()
    PREMIERS = _eratosthene_filtre(SIEVE_LIMIT)
    print(f"{len(PREMIERS):,} premiers utiles  ({time.perf_counter()-t0:.2f}s)")
    print(f"Cœurs disponibles : {n_coeurs}  →  {N_WORKERS} workers  (batch {BATCH})")

    # Centres — un par ratio
    # a_centre = (D − log10(23)) / (2·(log10(2) + r·log10(23)))
    a_max   = int((DIGITS_CIBLE - LOG10_K) / (2 * LOG10_2)) + 10
    centres = []
    for r in TOUS_RATIOS:
        ac = max(1, int((DIGITS_CIBLE - LOG10_K) / (2 * (LOG10_2 + r * LOG10_K))))
        bc = max(1, int(ac * r))
        centres.append((ac, bc))

    N_TOTAL = len(TOUS_RATIOS)

    # État initial
    cp = charger_checkpoint()
    if cp:
        raw_deltas = cp["deltas"]
        deltas   = [None if d == -1 else d for d in raw_deltas]
        sides    = cp["sides"]
        n_crible = cp["n_crible"]
        n_mr     = cp["n_mr"]
        n_paires = cp["n_paires"]
        t_cumul  = cp["t_cumul"]
        n_act    = n_actifs_pour(n_mr)
        print(f"\n  ♻  REPRISE v8 : {n_act}/{N_TOTAL} ratios actifs  |  "
              f"MR={n_mr}  cumulé={t_cumul/3600:.2f}h")
    else:
        deltas   = [0] * N_TOTAL
        sides    = [0] * N_TOTAL
        n_crible = n_mr = n_paires = 0
        t_cumul  = 0.0
        n_act    = N_RATIOS_INIT
        print(f"\n  Nouvelle recherche v8 (aucun checkpoint)")

    print(f"\n{'═'*72}")
    print(f"  PrimeQuest v8 — p = 23·m·(m+1)+1,  m = 2^a·23^b−1")
    print(f"{'═'*72}")
    print(f"  Cible            : {DIGITS_CIBLE:,} chiffres  (tolérance ±{TOLERANCE})")
    print(f"  Densité Bateman  : C=5.839 (vs 3.361 pour k=3) ≈ ×1.74  [HEURISTIQUE]")
    print(f"  F > √p           : F/√p ≥ √23 ≈ 4.796  [RIGOUREUX, INCONDITIONNEL]")
    print(f"  Blindspots       : {sorted(BLINDSPOTS)}")
    print(f"  Bornes ratios    : [{_RATIO_MIN_ABSOLU:.6f},  {_RATIO_MAX_ABSOLU:.2f}]")
    print(f"  Ratio optimal    : b/a = log10(2)/log10(23) ≈ 0.2211")
    print(f"  Spread initial   : ±{SPREAD_INIT}  ({N_RATIOS_INIT} ratios)  →"
          f"  [{RATIO_CENTRE-SPREAD_INIT:.3f}, {RATIO_CENTRE+SPREAD_INIT:.3f}]")
    print(f"  Spread max       : ±{SPREAD_MAX}  ({N_RATIOS_MAX} ratios)  →"
          f"  [{RATIO_CENTRE-SPREAD_MAX:.3f}, {RATIO_CENTRE+SPREAD_MAX:.3f}]")
    print(f"  Expansion        : +2 ratios tous les {ADAPT_INTERVAL} tests MR")
    print(f"  Ratios (tous)    : {TOUS_RATIOS}")
    print(f"  Ratios actifs    : {TOUS_RATIOS[:n_act]}  ({n_act}/{N_TOTAL})")
    print(f"  Centres (a₀,b₀) : {centres[:n_act]}")
    print(f"  a_max            : {a_max:,}")
    print(f"  Workers          : {N_WORKERS}  (batch {BATCH} paires)")
    print(f"  Crible           : {len(PREMIERS):,} premiers utiles ≤ {SIEVE_LIMIT:,}")
    print(f"  Témoins Pockl.   : q ∈ {{2, 23}}  [RIGOUREUX]")
    print(f"  MR rounds        : {MR_TOURS} + {MR_TOURS_CONFIRM}")
    print(f"  Timeout          : {TIMEOUT_S:,}s ({TIMEOUT_S/3600:.1f}h)")
    print(f"  Checkpoint       : toutes les {CHECKPOINT_S}s → {CHECKPOINT_FILE}")
    print(f"{'─'*72}\n")

    trouve    = None
    t_debut   = time.perf_counter()
    t_rapport = t_debut
    t_ckpt    = t_debut
    n_act_prec = n_act

    def _handler(sig, frame):
        print(f"\n  ⚡ Interruption — sauvegarde…", flush=True)
        sauver_checkpoint(deltas, sides, n_crible, n_mr, n_paires,
                          t_cumul + (time.perf_counter() - t_debut))
        print(f"  Checkpoint → {CHECKPOINT_FILE}")
        sys.exit(0)

    signal.signal(signal.SIGINT,  _handler)
    signal.signal(signal.SIGTERM, _handler)

    with mp.Pool(N_WORKERS,
                 initializer=_init_worker,
                 initargs=(PREMIERS,)) as pool:

        while any(d is not None for d in deltas[:n_act]):

            now           = time.perf_counter()
            elapsed_s     = now - t_debut
            elapsed_total = t_cumul + elapsed_s

            # Expansion spread
            n_act = n_actifs_pour(n_mr)
            if n_act > n_act_prec:
                nouveaux = TOUS_RATIOS[n_act_prec:n_act]
                print(f"\n  ↗  EXPANSION spread : {n_act_prec}→{n_act} ratios  "
                      f"(+{nouveaux})  après {n_mr} tests MR", flush=True)
                n_act_prec = n_act

            if elapsed_s >= TIMEOUT_S:
                print(f"\n  ⏱  Timeout {TIMEOUT_S/3600:.1f}h atteint.")
                sauver_checkpoint(deltas, sides, n_crible, n_mr, n_paires, elapsed_total)
                print(f"  Checkpoint → {CHECKPOINT_FILE}")
                break

            if now - t_ckpt >= CHECKPOINT_S:
                sauver_checkpoint(deltas, sides, n_crible, n_mr, n_paires, elapsed_total)
                t_ckpt = now

            if now - t_rapport >= RAPPORT_S:
                restant  = max(0, TIMEOUT_S - elapsed_s)
                taux     = elapsed_s / max(n_mr, 1)
                spread_c = SPREAD_INIT + (n_mr // ADAPT_INTERVAL) * (
                    (SPREAD_MAX - SPREAD_INIT) / max(1, (N_RATIOS_MAX - N_RATIOS_INIT) // 2)
                )
                print(
                    f"  [{elapsed_total/3600:5.2f}h]  "
                    f"ratios={n_act}/{N_TOTAL}  spread=±{min(spread_c, SPREAD_MAX):.3f}  "
                    f"paires={n_paires:,}  crible={n_crible:,}  "
                    f"MR={n_mr}  élim={n_crible/max(n_paires,1)*100:.1f}%  "
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
                if r['statut'] == 'crible':
                    n_crible += 1
                elif r['statut'] == 'compose':
                    n_mr += 1
                    now2    = time.perf_counter()
                    et      = t_cumul + (now2 - t_debut)
                    restant = max(0, TIMEOUT_S - (now2 - t_debut))
                    ratio_r = r['b'] / r['a']
                    print(
                        f"  [{et/3600:5.2f}h]  MR #{n_mr:4d}  "
                        f"a={r['a']}, b={r['b']}  r={ratio_r:.4f}  "
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
                        f"a={r['a']}, b={r['b']}  r={ratio_r:.4f}  "
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
        deux        = gmpy2.mpz(2)
        vingt_trois = gmpy2.mpz(K)
        F       = deux**a * vingt_trois**(b + 1)
        sp      = str(p)
        ratio_r = b / a

        print(f"  PREMIER DE {nb} CHIFFRES — PRIMALITÉ PROUVÉE ✅")
        print(f"{'═'*72}")
        print(f"""
  Structure
  ──────────────────────────────────────────────────────
  Famille  : p = 23·m·(m+1)+1,  m = 2^a·23^b−1
  a = {a},  b = {b}  (ratio b/a = {ratio_r:.4f})
  m = 2^{a}·23^{b}−1   ({gmpy2.num_digits(m)} chiffres)
  p = 23·m·(m+1)+1    ({nb} chiffres)

  Preuve Pocklington  [RIGOUREUSE, INCONDITIONNELLE]
  ──────────────────────────────────────────────────────
  F = 2^{a}·23^{{{b+1}}}  ({gmpy2.num_digits(F)} chiffres)
  F > √p : F²/p ≥ 23 > 1  ✅
  q=2   →  témoin w = {temoins[2]}   ✅
  q=23  →  témoin w = {temoins[23]}  ✅

  Performance
  ──────────────────────────────────────────────────────
  Temps session      : {t_session:.1f}s  ({t_session/3600:.2f}h)
  Temps total cumulé : {t_total:.1f}s  ({t_total/3600:.2f}h)
  Tests MR           : {n_mr}
  Éliminés crible    : {n_crible}  ({n_crible/max(n_paires,1)*100:.1f}%)
  Workers parallèles : {N_WORKERS}
  Tolérance          : ±{TOLERANCE} chiffres
  Ratios actifs      : {TOUS_RATIOS[:n_act]}

  p = {sp[:40]}…{sp[-20:]}
""")
        nom = f"premier_{nb}chiffres_v8.txt"
        with open(nom, "w", encoding="utf-8") as f:
            f.write(f"Premier de {nb} chiffres — PrimeQuest v8 (k=23)\n")
            f.write(f"Famille : p = 23·m·(m+1)+1,  m = 2^a·23^b−1\n")
            f.write(f"a={a}, b={b}  (ratio b/a={ratio_r:.4f})\n")
            f.write(f"Témoins : q=2→w={temoins[2]}, q=23→w={temoins[23]}\n")
            f.write(f"Temps   : {t_total:.1f}s  ({t_total/3600:.2f}h)\n")
            f.write(f"Workers : {N_WORKERS}\n")
            f.write(f"Tolérance : ±{TOLERANCE} chiffres\n")
            f.write(f"Ratios  : {TOUS_RATIOS[:n_act]}\n\nm =\n{m}\n\np =\n{p}\n")
        print(f"  Sauvegardé → {nom}")

    elif not any(d is not None for d in deltas[:n_act]):
        supprimer_checkpoint()
        print("  Espace (a,b) entièrement exploré — augmenter TOLERANCE ou SPREAD_MAX.")

    else:
        print(f"  COMPTE RENDU — Session terminée (relancer pour continuer)")
        print(f"{'═'*72}")
        print(f"""
  Temps session           : {t_session/3600:.2f}h
  Temps total cumulé      : {t_total/3600:.2f}h
  Ratios actifs           : {n_act}/{N_TOTAL}
  Paires testées          : {n_paires:,}
  Éliminées crible        : {n_crible:,}  ({n_crible/max(n_paires,1)*100:.1f}%)
  Tests MR                : {n_mr}
  Temps/MR effectif       : {t_session/max(n_mr,1):.1f}s  (÷ {N_WORKERS} workers)
  Tolérance               : ±{TOLERANCE} chiffres
  Checkpoint              : {CHECKPOINT_FILE} ✅
""")
