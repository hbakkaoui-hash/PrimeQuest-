#!/usr/bin/env python3
"""
PrimeQuest v2 — Recherche et preuve de primalité AVEC REPRISE
Famille : p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Preuve  : Théorème de Pocklington N-1

Optimisations intégrées (issues de l'analyse arithmétique) :
  1. Crible modulaire restreint aux premiers q ≡ 1 (mod 3)
     Théorème prouvé : q | p impossible si q ≡ 2 (mod 3)
     → ~50% des premiers du crible supprimés, taux élimination identique
  2. Zigzag depuis a_centre = (D/2)/log10(6) — zone la plus fertile
  3. Checkpoint automatique — reprise exacte après interruption

Usage : python3 primequest_final.py
        Modifier DIGITS_CIBLE, TIMEOUT_S et DELTA_DEPART ci-dessous.

Dépendance : pip install gmpy2
"""

import gmpy2
import math
import time
import sys
import json
import os
import signal

# ═══════════════════════════════════════════════════════════════
# PARAMÈTRES — modifier ici
# ═══════════════════════════════════════════════════════════════
DIGITS_CIBLE    = 20_000      # nombre de chiffres souhaité
TIMEOUT_S       = 14_400      # durée max par session (14400 = 4h)
TOLERANCE       = 10          # tolérance ±chiffres
MR_TOURS        = 25          # tours Miller-Rabin
TEMOINS_MAX     = 300         # max témoins Pocklington
SIEVE_LIMIT     = 1_000_000   # crible jusqu'à ce nombre
RAPPORT_S       = 300         # rapport toutes les 5 min
CHECKPOINT_S    = 60          # sauvegarde checkpoint toutes les 60s
CHECKPOINT_FILE = f"checkpoint_{DIGITS_CIBLE}.json"

# ───────────────────────────────────────────────────────────────
# DÉPART MANUEL — utile si le checkpoint n'existe pas
# ───────────────────────────────────────────────────────────────
# DELTA_DEPART > 0 : démarre à a_centre + DELTA_DEPART (ignore checkpoint)
# DELTA_DEPART = 0 : comportement automatique (checkpoint ou début)
DELTA_DEPART    = 0           # ← CHANGER ICI
TEMPS_CUMUL_H   = 0.0         # ← heures déjà passées (sessions précédentes)
# ═══════════════════════════════════════════════════════════════

deux    = gmpy2.mpz(2)
trois   = gmpy2.mpz(3)
LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)


# ── Crible d'Ératosthène — uniquement q ≡ 1 (mod 3) ────────────────────────
#
# Théorème (prouvé) : q | p est impossible si q ≡ 2 (mod 3).
# Preuve : 3x²-3x+1 ≡ 0 (mod q) a des solutions ssi Δ=-3 est un carré mod q,
#          ssi (-3/q)=1, ssi q ≡ 1 (mod 3).
# Conséquence : on ne conserve que les q ≡ 1 (mod 3), soit ~50% du crible.

def _eratosthene_filtre(lim):
    """Crible d'Ératosthène — retourne uniquement les premiers q ≡ 1 (mod 3)."""
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


# ── Crible modulaire (sans calculer p) ───────────────────────────────────────
#
# Pour q ≡ 1 (mod 3) : teste si q | p en O(1) sans calculer p.
# x = 2^a · 3^b mod q  (via petit théorème de Fermat)
# p mod q = (3(x-1)x + 1) mod q

def crible_modulaire(a, b):
    for q in PREMIERS:
        x     = (pow(2, a % (q - 1), q) * pow(3, b % (q - 1), q)) % q
        p_mod = (3 * (x - 1) * x + 1) % q
        if p_mod == 0:
            return False
    return True


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


# ── Zigzag depuis le centre (reprennable) ────────────────────────────────────
#
# Centre optimal : a_centre = (D/2) / log10(6)
# Séquence : centre, centre+1, centre-1, centre+2, centre-2, …
# Position encodée par (delta, side) pour reprise exacte via checkpoint.

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


# ── Checkpoint ────────────────────────────────────────────────────────────────

def sauver_checkpoint(delta, side, n_crible, n_mr, n_paires, t_cumul):
    data = {
        "digits":   DIGITS_CIBLE,
        "delta":    delta,
        "side":     side,
        "n_crible": n_crible,
        "n_mr":     n_mr,
        "n_paires": n_paires,
        "t_cumul":  t_cumul,
    }
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump(data, f, indent=2)

def charger_checkpoint():
    if not os.path.exists(CHECKPOINT_FILE):
        return None
    try:
        with open(CHECKPOINT_FILE) as f:
            data = json.load(f)
        if data.get("digits") != DIGITS_CIBLE:
            print(f"  ⚠  Checkpoint ignoré (digits={data.get('digits')} ≠ {DIGITS_CIBLE})")
            return None
        return data
    except Exception as e:
        print(f"  ⚠  Checkpoint illisible : {e}")
        return None

def supprimer_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        os.remove(CHECKPOINT_FILE)


# ── Gestionnaire de signal (Ctrl+C / SIGTERM) ────────────────────────────────

_etat_global = {}

def handler_signal(sig, frame):
    print(f"\n  ⚡ Interruption — sauvegarde du checkpoint…", flush=True)
    if _etat_global:
        sauver_checkpoint(
            _etat_global["delta"],
            _etat_global["side"],
            _etat_global["n_crible"],
            _etat_global["n_mr"],
            _etat_global["n_paires"],
            _etat_global["t_cumul"] + (time.perf_counter() - _etat_global["t_debut"])
        )
        print(f"  Checkpoint sauvegardé → {CHECKPOINT_FILE}")
        print(f"  Relancer le script pour reprendre.")
    sys.exit(0)

signal.signal(signal.SIGINT,  handler_signal)
signal.signal(signal.SIGTERM, handler_signal)


# ── Initialisation ────────────────────────────────────────────────────────────

a_centre = int((DIGITS_CIBLE / 2) / math.log10(6))
a_max    = int(DIGITS_CIBLE / (2 * LOG10_2)) + 10

if DELTA_DEPART > 0:
    delta    = DELTA_DEPART
    side     = 0
    n_crible = 0
    n_mr     = 0
    n_paires = 0
    t_cumul  = TEMPS_CUMUL_H * 3600
    print(f"\n  ▶  DÉPART MANUEL : delta={delta} → a={a_centre + delta}")
    print(f"     Temps cumulé : {TEMPS_CUMUL_H}h")
else:
    cp = charger_checkpoint()
    if cp:
        delta    = cp["delta"]
        side     = cp["side"]
        n_crible = cp["n_crible"]
        n_mr     = cp["n_mr"]
        n_paires = cp["n_paires"]
        t_cumul  = cp["t_cumul"]
        print(f"\n  ♻  REPRISE depuis le checkpoint :")
        print(f"     delta={delta}, side={side}  →  a={a_depuis_position(a_centre, delta, side)}")
        print(f"     Temps cumulé  : {t_cumul/3600:.2f}h")
        print(f"     MR effectués  : {n_mr}")
    else:
        delta    = 0
        side     = 0
        n_crible = 0
        n_mr     = 0
        n_paires = 0
        t_cumul  = 0.0
        print(f"\n  Nouvelle recherche (aucun checkpoint)")

print(f"\n{'═'*65}")
print(f"  PrimeQuest v2 — p = 3m(m+1)+1,  m = 2^a·3^b−1")
print(f"{'═'*65}")
print(f"  Cible         : {DIGITS_CIBLE:,} chiffres (±{TOLERANCE})")
print(f"  a_centre      : {a_centre:,}   a_max : {a_max:,}")
print(f"  Crible        : {len(PREMIERS):,} premiers q ≡ 1(mod3) ≤ {SIEVE_LIMIT:,}")
print(f"  Timeout       : {TIMEOUT_S:,}s ({TIMEOUT_S/3600:.1f}h) par session")
print(f"  Checkpoint    : toutes les {CHECKPOINT_S}s → {CHECKPOINT_FILE}")
print(f"{'─'*65}\n")


# ── Boucle principale ─────────────────────────────────────────────────────────

trouve    = None
t_debut   = time.perf_counter()
t_rapport = t_debut
t_ckpt    = t_debut

_etat_global.update({
    "delta": delta, "side": side,
    "n_crible": n_crible, "n_mr": n_mr, "n_paires": n_paires,
    "t_cumul": t_cumul, "t_debut": t_debut
})

while delta is not None:

    elapsed_session = time.perf_counter() - t_debut
    elapsed_total   = t_cumul + elapsed_session

    if elapsed_session >= TIMEOUT_S:
        print(f"\n  ⏱  Timeout {TIMEOUT_S/3600:.1f}h atteint.")
        sauver_checkpoint(delta, side, n_crible, n_mr, n_paires, elapsed_total)
        print(f"  Checkpoint → {CHECKPOINT_FILE}  (relancer pour continuer)")
        break

    a = a_depuis_position(a_centre, delta, side)

    for b in b_valides(a, DIGITS_CIBLE, TOLERANCE):
        n_paires += 1

        now = time.perf_counter()

        # ── Checkpoint périodique ─────────────────────────────────────────
        if now - t_ckpt >= CHECKPOINT_S:
            sauver_checkpoint(delta, side, n_crible, n_mr, n_paires,
                              t_cumul + (now - t_debut))
            t_ckpt = now
            _etat_global.update({
                "delta": delta, "side": side,
                "n_crible": n_crible, "n_mr": n_mr,
                "n_paires": n_paires, "t_cumul": t_cumul, "t_debut": t_debut
            })

        # ── Rapport de progression ────────────────────────────────────────
        if now - t_rapport >= RAPPORT_S:
            elapsed_s   = now - t_debut
            elapsed_tot = t_cumul + elapsed_s
            restant     = max(0, TIMEOUT_S - elapsed_s)
            taux        = elapsed_s / max(n_mr, 1)
            print(
                f"  [{elapsed_tot/3600:5.2f}h]  a={a:7,}  "
                f"paires={n_paires:9,}  crible={n_crible:8,}  "
                f"MR={n_mr:4,}  ~{taux:.0f}s/MR  "
                f"restant={restant/60:.0f}min",
                flush=True
            )
            t_rapport = now

        # ── Crible modulaire — q ≡ 1 (mod 3) uniquement ──────────────────
        if not crible_modulaire(a, b):
            n_crible += 1
            continue

        # ── Miller-Rabin ──────────────────────────────────────────────────
        n_mr    += 1
        now      = time.perf_counter()
        elapsed_tot = t_cumul + (now - t_debut)
        restant  = max(0, TIMEOUT_S - (now - t_debut))
        print(
            f"  [{elapsed_tot/3600:5.2f}h]  MR #{n_mr:4d}  "
            f"a={a}, b={b}  (restant {restant/60:.0f}min) …",
            flush=True
        )

        m = deux**a * trois**b - 1
        p = 3 * m * (m + 1) + 1

        if not gmpy2.is_prime(p, MR_TOURS):
            print("  → composé", flush=True)
            continue

        # ── Pocklington ───────────────────────────────────────────────────
        print("  *** PROBABLE — Pocklington… ***", flush=True)
        ok, temoins = pocklington(p, a, b)
        if ok:
            trouve = (a, b, m, p, int(gmpy2.num_digits(p)), temoins)
            break
        else:
            manque = [q for q, v in temoins.items() if v is None]
            print(f"  ⚠  Pocklington incomplet — témoin manquant q={manque}")

    if trouve:
        break

    delta, side = position_suivante(delta, side, a_centre, 1, a_max)
    _etat_global["delta"] = delta
    _etat_global["side"]  = side

t_session = time.perf_counter() - t_debut
t_total   = t_cumul + t_session


# ── Résultat ──────────────────────────────────────────────────────────────────

print(f"\n{'═'*65}")

if trouve:
    supprimer_checkpoint()
    a, b, m, p, nb, temoins = trouve
    F  = deux**a * trois**(b + 1)
    sp = str(p)

    print(f"  PREMIER DE {nb} CHIFFRES — PRIMALITÉ PROUVÉE ✅")
    print(f"{'═'*65}")
    print(f"""
  Structure
  ─────────────────────────────────────────────────
  a = {a},  b = {b}
  m = 2^{a}·3^{b}−1         ({gmpy2.num_digits(m)} chiffres)
  p = 3m(m+1)+1            ({nb} chiffres)

  Preuve Pocklington
  ─────────────────────────────────────────────────
  F = 2^{a}·3^{{{b+1}}}  ({gmpy2.num_digits(F)} chiffres) > √p  ✅
  q=2  →  témoin w = {temoins[2]}   ✅
  q=3  →  témoin w = {temoins[3]}   ✅

  Performance
  ─────────────────────────────────────────────────
  Temps session     : {t_session:.1f}s  ({t_session/60:.1f} min)
  Temps total cumulé: {t_total:.1f}s  ({t_total/3600:.2f}h)
  Tests MR          : {n_mr}
  Éliminés (crible) : {n_crible}

  p = {sp[:40]}…{sp[-20:]}
""")

    nom = f"premier_{nb}chiffres.txt"
    with open(nom, "w") as f:
        f.write(f"Premier de {nb} chiffres — PrimeQuest v2\n")
        f.write(f"a={a}, b={b}\n")
        f.write(f"Témoins : q=2→w={temoins[2]}, q=3→w={temoins[3]}\n")
        f.write(f"Temps   : {t_total:.1f}s\n\nm =\n{m}\n\np =\n{p}\n")
    print(f"  Sauvegardé → {nom}")

elif delta is None:
    supprimer_checkpoint()
    print("  Espace (a,b) entièrement exploré — aucun premier trouvé.")
    print("  Augmenter TOLERANCE.")

else:
    tps_par_mr = t_session / max(n_mr, 1)
    print(f"  COMPTE RENDU — Session terminée")
    print(f"{'═'*65}")
    print(f"""
  Temps session           : {t_session/3600:.2f}h
  Temps total cumulé      : {t_total/3600:.2f}h
  a exploré               : a_centre={a_centre} ± {delta}
  Paires testées          : {n_paires:,}
  Éliminées (crible)      : {n_crible:,}  ({n_crible/max(n_paires,1)*100:.1f}%)
  Tests MR                : {n_mr}
  Temps par MR            : {tps_par_mr:.1f}s
  Résultat                : aucun premier — relancer pour continuer

  Checkpoint              : {CHECKPOINT_FILE} ✅
""")
