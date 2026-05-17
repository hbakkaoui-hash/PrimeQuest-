#!/usr/bin/env python3
"""
PrimeQuest v9 — Famille k=23, trois optimisations cœur
═══════════════════════════════════════════════════════
Famille : p = 23·m·(m+1)+1,   m = 2^a · 23^b − 1
Preuve  : Théorème de Pocklington N-1  (témoins q ∈ {2, 23})

Optimisations v9 vs v8 (mathématiquement équivalent) :
───────────────────────────────────────────────────────
  [v9-1] SINGLE MR CALL
         v8 faisait is_prime(p,3) puis is_prime(p,20) : le second refait
         entièrement BPSW+20 rounds même si le premier avait déjà répondu.
         Pour les composites (>99 % des cas après crible), les deux calls
         échouent au BPSW interne → coût identique. Pour les pseudo-premiers
         (cas rarissime), le double appel était redondant.
         Fix : un seul is_prime(p, MR_TOURS_CONFIRM).
         Gain attendu : faible (BPSW domine, pas les rounds supplémentaires),
         mais élimine une incohérence et allège le code.

  [v9-2] CRIBLE PAR GRILLE (grid sieve) pour q ≤ GRID_LIMIT
         Observation : m = 2^a · 23^b − 1 est périodique en a mod ord_q(2)
         et en b mod ord_q(23). On peut donc pré-calculer UNE FOIS, pour
         chaque q, l'ensemble des paires interdites
             FORBIDDEN_q = {(i,j) : p ≡ 0 mod q  avec  i=a%r2, j=b%r23}
         puis tester chaque paire (a,b) par simple lookup O(1) au lieu de
         deux appels pow() Python.
         Équivalence : 2^a ≡ 2^(a mod r2) mod q par définition de l'ordre
         multiplicatif (r2 = ord_q(2) divise q-1 ; Fermat donne aussi
         2^(a%(q-1)) ≡ 2^a mod q et ces deux expressions coïncident mod q).
         Pour q > GRID_LIMIT : on conserve le pow() classique (tables trop
         grandes à pré-calculer).
         Gain attendu : ×5-10 sur les q ≤ GRID_LIMIT ; négligeable sur le
         temps total (MR domine à grande taille), mais prépare une réduction
         future de SIEVE_LIMIT guidée par l'instrumentation.

  [v9-3] MÉMOÏSATION DE 2^a PAR WORKER
         b_valides() retourne plusieurs b pour chaque a ; le worker
         recalculait 2^a à chaque appel. On cache les derniers résultats
         dans un dict process-local (taille ≤ CACHE_POW2A_SIZE).
         Équivalence : triviale (aucune logique modifiée).
         Gain attendu : faible pour les grands a (gmpy2 exponentiation
         rapide), mais mesurable quand TOLERANCE est grande (5+ b par a).

  [v9-4] INSTRUMENTATION t_crible / t_mr
         Chaque worker renvoie ses temps de crible et MR. Le main agrège
         et affiche le ratio pour guider le choix optimal de SIEVE_LIMIT.

Propriétés mathématiques inchangées (contraintes absolues respectées) :
  · Blindspots {2,3,5,7,11,13,17,23,29,31,41,43,59} : intacts
  · Discriminant 437 = 19·23 : non modifié
  · Certificat Pocklington q∈{2,23}, F>√p : non modifié
  · Checkpointing v9 : format identique à v8 sauf "version":9
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
K                = 23
DIGITS_CIBLE     = 50_000
TIMEOUT_S        = 14_400
TOLERANCE        = 50
MR_TOURS_CONFIRM = 15          # [v9-1] un seul appel (était 3+20 en v8)
TEMOINS_MAX      = 300
SIEVE_LIMIT      = 10_000_000  # tunable via instrumentation t_crible/t_mr
GRID_LIMIT       = 1_000       # [v9-2] grid lookup pour q ≤ GRID_LIMIT
CACHE_POW2A_SIZE = 8           # [v9-3] nb max d'entrées dans le cache 2^a
RAPPORT_S        = 300
CHECKPOINT_S     = 60
CHECKPOINT_FILE  = f"checkpoint_{DIGITS_CIBLE}_v9.json"

RATIO_CENTRE    = 0.221
N_RATIOS_INIT   = 7
N_RATIOS_MAX    = 15
SPREAD_INIT     = 0.08
SPREAD_MAX      = 0.18
ADAPT_INTERVAL  = 200

N_WORKERS_PARAM = 0

# ═══════════════════════════════════════════════════════════════

LOG10_2 = math.log10(2)
LOG10_K = math.log10(K)

_RATIO_MIN_ABSOLU = 2 * LOG10_2 / DIGITS_CIBLE
_RATIO_MAX_ABSOLU = (DIGITS_CIBLE / 2 - LOG10_2) / LOG10_K

BLINDSPOTS = frozenset({2, 3, 5, 7, 11, 13, 17, 23, 29, 31, 41, 43, 59})


# ── Génération TOUS_RATIOS ────────────────────────────────────────────────────

def _generer_tous_ratios():
    step = (2 * SPREAD_MAX / (N_RATIOS_MAX - 1)) if N_RATIOS_MAX > 1 else 0
    bruts = [round(RATIO_CENTRE - SPREAD_MAX + i * step, 4) for i in range(N_RATIOS_MAX)]
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
        raise ValueError(f"Seulement {len(valides)} ratios valides.")
    return valides

TOUS_RATIOS = _generer_tous_ratios()

def n_actifs_pour(n_mr_total):
    return min(N_RATIOS_MAX, N_RATIOS_INIT + 2 * (n_mr_total // ADAPT_INTERVAL))


# ── [v9-2] Construction du crible ────────────────────────────────────────────

def _eratosthene(lim):
    """Primes q ≤ lim, non-blindspot, où 437 est QR mod q."""
    t = bytearray([1]) * (lim + 1)
    t[0] = t[1] = 0
    for i in range(2, int(lim**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    result = []
    for q in range(3, lim + 1):
        if t[q] and q not in BLINDSPOTS and pow(437, (q - 1) // 2, q) == 1:
            result.append(q)
    return result

def _ordre_multiplicatif(base, q):
    """Ordre multiplicatif de base mod q (q premier, pgcd(base,q)=1)."""
    cur = base % q
    for k in range(1, q):
        if cur == 1:
            return k
        cur = cur * (base % q) % q
    return q - 1

def _build_grille(premiers_grille):
    """
    [v9-2] Pour chaque q ≤ GRID_LIMIT, pré-calcule :
      r2  = ord_q(2),  r23 = ord_q(23)
      FORBIDDEN_q = {(i,j) : p(2^i·23^j−1) ≡ 0 mod q}
    Retourne une liste de tuples (r2, r23, frozenset).

    Équivalence avec le crible pow() classique :
      2^a ≡ 2^(a mod r2) mod q  (définition de l'ordre multiplicatif)
      23^b ≡ 23^(b mod r23) mod q  (idem)
    Donc (a mod r2, b mod r23) ∈ FORBIDDEN_q  ⟺  q | p(a,b).
    """
    grille = []
    for q in premiers_grille:
        r2  = _ordre_multiplicatif(2, q)
        r23 = _ordre_multiplicatif(K, q)
        pow2  = [pow(2, i, q) for i in range(r2)]
        pow23 = [pow(K, j, q) for j in range(r23)]
        forbidden = frozenset(
            (i, j)
            for i in range(r2)
            for j in range(r23)
            if (K * ((pow2[i] * pow23[j] % q) - 1) % q
                * ((pow2[i] * pow23[j] % q)) % q + 1) % q == 0
        )
        if forbidden:
            grille.append((r2, r23, forbidden))
    return grille


# ── Worker ────────────────────────────────────────────────────────────────────

_W_GRILLE         = None   # [(r2, r23, frozenset), …] pour q ≤ GRID_LIMIT
_W_PREMIERS_LARGE = None   # [q, …] pour GRID_LIMIT < q ≤ SIEVE_LIMIT
_W_CACHE_POW2A    = {}     # [v9-3] cache 2^a par worker

def _init_worker(grille, premiers_large):
    global _W_GRILLE, _W_PREMIERS_LARGE, _W_CACHE_POW2A
    _W_GRILLE         = grille
    _W_PREMIERS_LARGE = premiers_large
    _W_CACHE_POW2A    = {}
    signal.signal(signal.SIGINT,  signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_IGN)

def _crible(a, b):
    """Retourne False si un premier du crible divise p(a,b)."""
    # Tier 1 : grid O(1) par q (q ≤ GRID_LIMIT)
    for r2, r23, forbidden in _W_GRILLE:
        if (a % r2, b % r23) in forbidden:
            return False
    # Tier 2 : pow classique (q > GRID_LIMIT)
    for q in _W_PREMIERS_LARGE:
        q1    = q - 1
        x     = pow(2, a % q1, q) * pow(K, b % q1, q) % q
        m_mod = (x - 1) % q
        if (K * m_mod * (m_mod + 1) + 1) % q == 0:
            return False
    return True

def _get_pow2a(a):
    """[v9-3] Retourne gmpy2.mpz(2)**a avec cache process-local."""
    if a not in _W_CACHE_POW2A:
        if len(_W_CACHE_POW2A) >= CACHE_POW2A_SIZE:
            # Éviction FIFO : supprime la plus ancienne clé
            _W_CACHE_POW2A.pop(next(iter(_W_CACHE_POW2A)))
        _W_CACHE_POW2A[a] = gmpy2.mpz(2) ** a
    return _W_CACHE_POW2A[a]

def _tester_paire(args):
    """[v9-1] Un seul appel is_prime (MR_TOURS_CONFIRM rounds)."""
    a, b = args

    t0 = time.perf_counter()
    if not _crible(a, b):
        return {'statut': 'crible', 'a': a, 'b': b,
                't_crible': time.perf_counter() - t0, 't_mr': 0.0}
    t_crible = time.perf_counter() - t0

    t1 = time.perf_counter()
    pow2a = _get_pow2a(a)
    m = pow2a * gmpy2.mpz(K) ** b - 1
    p = K * m * (m + 1) + 1

    if not gmpy2.is_prime(p, MR_TOURS_CONFIRM):
        return {'statut': 'compose', 'a': a, 'b': b,
                't_crible': t_crible, 't_mr': time.perf_counter() - t1}

    return {'statut': 'probable', 'a': a, 'b': b, 'p': p, 'm': m,
            't_crible': t_crible, 't_mr': time.perf_counter() - t1}


# ── b valides ─────────────────────────────────────────────────────────────────

def b_valides(a, digits, tol):
    b_c = (digits - 2 * a * LOG10_2 - LOG10_K) / (2 * LOG10_K)
    return [
        b for b in range(max(1, int(b_c) - tol - 2), int(b_c) + tol + 3)
        if abs(2 * a * LOG10_2 + (2 * b + 1) * LOG10_K - digits) <= tol + 2
    ]


# ── Pocklington ───────────────────────────────────────────────────────────────

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
            if gmpy2.gcd(gmpy2.powmod(gw, p1 // gmpy2.mpz(q), p) - 1, p) == 1:
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
        if centre + 1 <= a_max: return 1, 0
        if centre - 1 >= a_min: return 1, 1
        return None, None
    if side == 0:
        if centre - delta >= a_min: return delta, 1
        return position_suivante(delta + 1, -1, centre, a_min, a_max)
    nd = delta + 1
    if centre + nd <= a_max: return nd, 0
    if centre - nd >= a_min: return nd, 1
    return None, None


# ── Multi-ratio round-robin ───────────────────────────────────────────────────

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
                paires.append((a, b))
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
            "version":     9,
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
        if data.get("version") != 9:
            print(f"  ⚠  Checkpoint ignoré (version {data.get('version')} ≠ 9)")
            return None
        if data.get("k") != K:
            print(f"  ⚠  Checkpoint ignoré (k={data.get('k')} ≠ {K})")
            return None
        if data.get("tous_ratios") != TOUS_RATIOS:
            print(f"  ⚠  Checkpoint ignoré (ratios différents)")
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
    BATCH     = N_WORKERS * 4   # légèrement agrandi pour amortir le overhead

    # ── Construction du crible ────────────────────────────────────────────────
    print(f"Construction du crible (blindspots={len(BLINDSPOTS)}, "
          f"GRID_LIMIT={GRID_LIMIT:,}, SIEVE_LIMIT={SIEVE_LIMIT:,})…")

    t0 = time.perf_counter()
    tous_premiers = _eratosthene(SIEVE_LIMIT)
    premiers_grille = [q for q in tous_premiers if q <= GRID_LIMIT]
    premiers_large  = [q for q in tous_premiers if q > GRID_LIMIT]
    t_sieve = time.perf_counter() - t0

    print(f"  Crible Ératosthène : {len(tous_premiers):,} premiers utiles  ({t_sieve:.2f}s)")

    t0 = time.perf_counter()
    GRILLE = _build_grille(premiers_grille)
    t_grid = time.perf_counter() - t0
    print(f"  Grille [v9-2]      : {len(GRILLE)} entrées (q ≤ {GRID_LIMIT}) "
          f" en {t_grid:.2f}s")
    print(f"  Pow classique      : {len(premiers_large):,} premiers "
          f"({GRID_LIMIT} < q ≤ {SIEVE_LIMIT:,})")
    print(f"Cœurs : {n_coeurs}  →  {N_WORKERS} workers  (batch {BATCH})")

    # ── Centres ───────────────────────────────────────────────────────────────
    a_max   = int((DIGITS_CIBLE - LOG10_K) / (2 * LOG10_2)) + 10
    centres = []
    for r in TOUS_RATIOS:
        ac = max(1, int((DIGITS_CIBLE - LOG10_K) / (2 * (LOG10_2 + r * LOG10_K))))
        bc = max(1, int(ac * r))
        centres.append((ac, bc))
    N_TOTAL = len(TOUS_RATIOS)

    # ── Reprise checkpoint ─────────────────────────────────────────────────────
    cp = charger_checkpoint()
    if cp:
        deltas   = [None if d == -1 else d for d in cp["deltas"]]
        sides    = cp["sides"]
        n_crible = cp["n_crible"]
        n_mr     = cp["n_mr"]
        n_paires = cp["n_paires"]
        t_cumul  = cp["t_cumul"]
        n_act    = n_actifs_pour(n_mr)
        print(f"\n  ♻  REPRISE v9 : {n_act}/{N_TOTAL} ratios  |  "
              f"MR={n_mr}  cumulé={t_cumul/3600:.2f}h")
    else:
        deltas   = [0] * N_TOTAL
        sides    = [0] * N_TOTAL
        n_crible = n_mr = n_paires = 0
        t_cumul  = 0.0
        n_act    = N_RATIOS_INIT
        print(f"\n  Nouvelle recherche v9 (aucun checkpoint)")

    # Accumulateurs instrumentation [v9-4]
    acc_t_crible = 0.0
    acc_t_mr     = 0.0

    print(f"\n{'═'*72}")
    print(f"  PrimeQuest v9 — p = 23·m·(m+1)+1,  m = 2^a·23^b−1")
    print(f"{'═'*72}")
    print(f"  Cible           : {DIGITS_CIBLE:,} chiffres  (±{TOLERANCE})")
    print(f"  MR rounds       : {MR_TOURS_CONFIRM}  (single call)  [v9-1]")
    print(f"  Grid sieve      : {len(GRILLE)} primes ≤ {GRID_LIMIT} → O(1)/q  [v9-2]")
    print(f"  Cache 2^a       : {CACHE_POW2A_SIZE} entrées/worker  [v9-3]")
    print(f"  Instrum.        : t_crible/t_mr affichés  [v9-4]")
    print(f"  Spread          : ±{SPREAD_INIT}→±{SPREAD_MAX}  "
          f"({N_RATIOS_INIT}→{N_RATIOS_MAX} ratios, +2 / {ADAPT_INTERVAL} MR)")
    print(f"  Ratios actifs   : {TOUS_RATIOS[:n_act]}")
    print(f"  Centres (a₀)    : {[c[0] for c in centres[:n_act]]}")
    print(f"  Timeout         : {TIMEOUT_S/3600:.1f}h  |  ckpt toutes {CHECKPOINT_S}s")
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
                 initargs=(GRILLE, premiers_large)) as pool:

        while any(d is not None for d in deltas[:n_act]):

            now           = time.perf_counter()
            elapsed_s     = now - t_debut
            elapsed_total = t_cumul + elapsed_s

            n_act = n_actifs_pour(n_mr)
            if n_act > n_act_prec:
                nouveaux = TOUS_RATIOS[n_act_prec:n_act]
                print(f"\n  ↗  EXPANSION : {n_act_prec}→{n_act} ratios  "
                      f"(+{nouveaux})  après {n_mr} MR", flush=True)
                n_act_prec = n_act

            if elapsed_s >= TIMEOUT_S:
                print(f"\n  ⏱  Timeout {TIMEOUT_S/3600:.1f}h.")
                sauver_checkpoint(deltas, sides, n_crible, n_mr, n_paires, elapsed_total)
                break

            if now - t_ckpt >= CHECKPOINT_S:
                sauver_checkpoint(deltas, sides, n_crible, n_mr, n_paires, elapsed_total)
                t_ckpt = now

            if now - t_rapport >= RAPPORT_S:
                restant  = max(0, TIMEOUT_S - elapsed_s)
                taux_mr  = elapsed_s / max(n_mr, 1)
                elim_pct = n_crible / max(n_paires, 1) * 100
                ratio_tc = (acc_t_crible / max(acc_t_mr, 1e-9))
                spread_c = SPREAD_INIT + (n_mr // ADAPT_INTERVAL) * (
                    (SPREAD_MAX - SPREAD_INIT) / max(1, (N_RATIOS_MAX - N_RATIOS_INIT) // 2))
                print(
                    f"  [{elapsed_total/3600:5.2f}h]  "
                    f"ratios={n_act}/{N_TOTAL}  spread=±{min(spread_c,SPREAD_MAX):.3f}  "
                    f"paires={n_paires:,}  crible={n_crible:,}({elim_pct:.1f}%)  "
                    f"MR={n_mr}  ~{taux_mr:.0f}s/MR  "
                    f"t_crible/t_mr={ratio_tc:.3f}  restant={restant/60:.0f}min",
                    flush=True
                )
                # [v9-4] Conseil SIEVE_LIMIT
                if ratio_tc < 0.05 and n_mr >= 10:
                    print(f"  💡 t_crible/t_mr={ratio_tc:.3f} < 0.05 : "
                          f"MR domine, SIEVE_LIMIT optimal ou réductible.", flush=True)
                elif ratio_tc > 0.3 and n_mr >= 10:
                    print(f"  💡 t_crible/t_mr={ratio_tc:.3f} > 0.3 : "
                          f"crible coûteux, envisager SIEVE_LIMIT plus bas.", flush=True)
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
                acc_t_crible += r.get('t_crible', 0.0)
                acc_t_mr     += r.get('t_mr', 0.0)

                if r['statut'] == 'crible':
                    n_crible += 1
                elif r['statut'] == 'compose':
                    n_mr += 1
                    now2    = time.perf_counter()
                    et      = t_cumul + (now2 - t_debut)
                    restant = max(0, TIMEOUT_S - (now2 - t_debut))
                    print(
                        f"  [{et/3600:5.2f}h]  MR #{n_mr:4d}  "
                        f"a={r['a']}, b={r['b']}  r={r['b']/r['a']:.4f}  "
                        f"(restant {restant/60:.0f}min)  → composé",
                        flush=True
                    )
                elif r['statut'] == 'probable':
                    n_mr += 1
                    now2    = time.perf_counter()
                    et      = t_cumul + (now2 - t_debut)
                    restant = max(0, TIMEOUT_S - (now2 - t_debut))
                    print(
                        f"\n  [{et/3600:5.2f}h]  MR #{n_mr}  "
                        f"a={r['a']}, b={r['b']}  r={r['b']/r['a']:.4f}  "
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
        F       = gmpy2.mpz(2)**a * gmpy2.mpz(K)**(b + 1)
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
  t_crible/t_mr      : {acc_t_crible/max(acc_t_mr,1e-9):.3f}
  Workers            : {N_WORKERS}
  Tolérance          : ±{TOLERANCE} chiffres
  Ratios actifs      : {TOUS_RATIOS[:n_act]}

  p = {sp[:40]}…{sp[-20:]}
""")
        nom = f"premier_{nb}chiffres_v9.txt"
        with open(nom, "w", encoding="utf-8") as f:
            f.write(f"Premier de {nb} chiffres — PrimeQuest v9 (k=23)\n")
            f.write(f"Famille : p = 23·m·(m+1)+1,  m = 2^a·23^b−1\n")
            f.write(f"a={a}, b={b}  (ratio b/a={ratio_r:.4f})\n")
            f.write(f"Témoins : q=2→w={temoins[2]}, q=23→w={temoins[23]}\n")
            f.write(f"Temps   : {t_total:.1f}s  ({t_total/3600:.2f}h)\n")
            f.write(f"Workers : {N_WORKERS}\n")
            f.write(f"MR rounds : {MR_TOURS_CONFIRM}  (single call, v9)\n")
            f.write(f"t_crible/t_mr : {acc_t_crible/max(acc_t_mr,1e-9):.3f}\n")
            f.write(f"Tolérance : ±{TOLERANCE} chiffres\n")
            f.write(f"Ratios  : {TOUS_RATIOS[:n_act]}\n\nm =\n{m}\n\np =\n{p}\n")
        print(f"  Sauvegardé → {nom}")

    elif not any(d is not None for d in deltas[:n_act]):
        supprimer_checkpoint()
        print("  Espace exploré entièrement — augmenter TOLERANCE ou SPREAD_MAX.")

    else:
        elim_pct = n_crible / max(n_paires, 1) * 100
        print(f"  COMPTE RENDU — Session terminée (relancer pour continuer)")
        print(f"{'═'*72}")
        print(f"""
  Temps session           : {t_session/3600:.2f}h
  Temps total cumulé      : {t_total/3600:.2f}h
  Ratios actifs           : {n_act}/{N_TOTAL}
  Paires testées          : {n_paires:,}
  Éliminées crible        : {n_crible:,}  ({elim_pct:.1f}%)
  Tests MR                : {n_mr}
  Temps/MR effectif       : {t_session/max(n_mr,1):.1f}s  (÷{N_WORKERS} workers)
  t_crible / t_mr         : {acc_t_crible/max(acc_t_mr,1e-9):.3f}
  → Si ratio < 0.05 : MR domine, SIEVE_LIMIT bien réglé ou réductible.
  → Si ratio > 0.30 : crible trop lourd, réduire SIEVE_LIMIT.
  Checkpoint              : {CHECKPOINT_FILE} ✅
""")
