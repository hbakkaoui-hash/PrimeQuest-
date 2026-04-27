#!/usr/bin/env python3
"""
PrimeQuest-HPC — Première ébauche
Recherche de premiers certifiés de l'ordre de plusieurs millions de chiffres
Conçu pour supercalculateur (MPI) ou cluster multi-cœurs (multiprocessing)

Famille cible : p = 3m(m+1)+1,  m = 2^a · 3^b − 1
Preuve        : Pocklington N-1  (p−1 = 2^a · 3^(b+1) · m, complètement factorisé)

Architecture :
  ┌─────────────────────────────────────────────────────────┐
  │  Rang 0 (coordinateur)                                  │
  │    - distribue les plages (a,b) aux workers             │
  │    - reçoit les candidats probables                     │
  │    - applique la preuve Pocklington                     │
  │    - écrit le résultat final                            │
  ├─────────────────────────────────────────────────────────┤
  │  Rangs 1..N (workers)                                   │
  │    - crible modulaire (µs/candidat, sans calculer p)    │
  │    - Miller-Rabin sur les survivants (min à h/candidat) │
  │    - envoie les probables au rang 0                     │
  └─────────────────────────────────────────────────────────┘

Usage (MPI)          : mpirun -np 1024 python3 primequest_hpc_draft.py --digits 1000000
Usage (local, test)  : python3 primequest_hpc_draft.py --digits 10000 --mode local
"""

import gmpy2
import math
import sys
import time
import random
import argparse
import os
from typing import Optional

# ── Import MPI (optionnel — fallback multiprocessing) ────────────────────────
try:
    from mpi4py import MPI
    MPI_OK = True
except ImportError:
    MPI_OK = False

# ── Fallback multiprocessing pour tests locaux ────────────────────────────────
import multiprocessing as mp

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — adapter selon la machine cible
# ═══════════════════════════════════════════════════════════════════════════════

DIGITS_CIBLE       = 1_000_000   # objectif (chiffres décimaux)
TOLERANCE          = 30          # fenêtre ±digits (plus large pour grands nombres)
MR_TOURS           = 20          # tours Miller-Rabin (20 = sécurité 2^-40)
TEMOINS_MAX        = 500         # max témoins cherchés pour Pocklington
SIEVE_LIMIT        = 10_000_000  # crible modulaire jusqu'à 10^7 (~664 000 premiers)
CHECKPOINT_ITER    = 500         # log tous les N candidats passés au crible

# ── Constantes gmpy2 ──────────────────────────────────────────────────────────
deux  = gmpy2.mpz(2)
trois = gmpy2.mpz(3)

# ═══════════════════════════════════════════════════════════════════════════════
# BLOC 1 : CRIBLE D'ÉRATOSTHÈNE (exécuté une seule fois au démarrage)
# ═══════════════════════════════════════════════════════════════════════════════

def _eratosthene(limite: int) -> list:
    t = bytearray([1]) * (limite + 1)
    t[0] = t[1] = 0
    for i in range(2, int(limite**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    return [i for i in range(2, limite + 1) if t[i]]

print(f"[init] Construction du crible jusqu'à {SIEVE_LIMIT:,}…", flush=True)
t_s = time.perf_counter()
PREMIERS_PETITS = _eratosthene(SIEVE_LIMIT)
print(f"[init] {len(PREMIERS_PETITS):,} premiers ({time.perf_counter()-t_s:.1f}s)", flush=True)

# ═══════════════════════════════════════════════════════════════════════════════
# BLOC 2 : CRIBLE MODULAIRE  ← innovation clé pour les très grands nombres
#
# Pour p de plusieurs millions de chiffres, calculer p explicitement rien que
# pour le cribler prendrait des minutes. On exploite la structure algébrique :
#
#   p = 3·(2^a·3^b − 1)·(2^a·3^b) + 1
#   p mod q = (3·(2^a mod q·3^b mod q − 1)·(2^a mod q·3^b mod q) + 1) mod q
#
# Par le petit théorème de Fermat : 2^a mod q = 2^(a mod (q−1)) mod q
# Coût : O(log q) par premier q — négligeable, même pour 664 000 premiers.
# Gain : élimine >99% des composés AVANT de calculer p.
# ═══════════════════════════════════════════════════════════════════════════════

def crible_modulaire(a: int, b: int) -> bool:
    """
    Teste p = 3m(m+1)+1 contre tous les petits premiers
    sans jamais calculer p entièrement.
    Retourne False (composé certain) ou True (passe le crible).
    """
    for q in PREMIERS_PETITS:
        # p ≡ 1 (mod 2) toujours (p impair)
        # p ≡ 1 (mod 3) toujours (par construction)
        if q <= 3:
            continue

        # Petit théorème de Fermat : ordre de 2 mod q divise q-1
        # → 2^a ≡ 2^(a mod (q-1)) (mod q)
        a_red = a % (q - 1)
        b_red = b % (q - 1)
        pow2a = pow(2, a_red, q)   # pow() natif Python : rapide
        pow3b = pow(3, b_red, q)

        m_mod = (pow2a * pow3b - 1) % q
        p_mod = (3 * m_mod * (m_mod + 1) + 1) % q

        if p_mod == 0:
            return False   # composé certain

    return True   # survivant — passe au Miller-Rabin

# ═══════════════════════════════════════════════════════════════════════════════
# BLOC 3 : CALCUL DE b POUR UN a DONNÉ  (sans calcul en grande précision)
#
# Pour un exposant a fixé, la valeur de b telle que p ait environ D chiffres :
#   log10(p) ≈ 2·(a·log10(2) + b·log10(3)) + log10(3)
#   → b ≈ (D/2 − a·log10(2) − log10(3)/2) / log10(3)
#
# On ne calcule que quelques entiers autour de cette valeur approchée,
# ce qui permet de parcourir le domaine (a,b) en O(1) par valeur de a.
# ═══════════════════════════════════════════════════════════════════════════════

_LOG10_2 = math.log10(2)
_LOG10_3 = math.log10(3)

def b_approxime(a: int, digits: int) -> list:
    """
    Retourne la liste des b (typiquement 2-4 valeurs) tels que
    p = 3m(m+1)+1 a ≈ digits ± TOLERANCE chiffres.
    Calcul entièrement en virgule flottante — pas de grande précision.
    """
    b_centre = (digits / 2 - a * _LOG10_2 - _LOG10_3 / 2) / _LOG10_3
    candidats = []
    for b in range(max(1, int(b_centre) - 4), int(b_centre) + 5):
        nb_approx = 2 * (a * _LOG10_2 + b * _LOG10_3) + _LOG10_3
        if abs(nb_approx - digits) <= TOLERANCE + 5:
            candidats.append(b)
    return candidats

# ═══════════════════════════════════════════════════════════════════════════════
# BLOC 4 : MILLER-RABIN SUR LE NOMBRE COMPLET
#
# C'est l'étape coûteuse. Pour 1 million de chiffres (3,3 Mbits) :
#   - 1 exponentiation modulaire ≈ quelques minutes sur 1 cœur
#   - 20 tours MR ≈ quelques heures sur 1 cœur
#
# Sur supercalculateur, plusieurs stratégies :
#   (A) Chaque nœud teste un candidat différent — embarrassingly parallel ✅
#   (B) Paralléliser les tours MR entre threads du même nœud
#   (C) Paralléliser la multiplication elle-même (FFT distribuée) — complexe
#
# L'approche (A) est adoptée ici car elle maximise le débit global.
# ═══════════════════════════════════════════════════════════════════════════════

def miller_rabin_complet(a: int, b: int):
    """
    Calcule p entièrement (coûteux !) puis applique Miller-Rabin.
    Retourne (True, p) ou (False, None).
    """
    m = deux**a * trois**b - 1
    p = 3 * m * (m + 1) + 1
    if gmpy2.is_prime(p, MR_TOURS):
        return True, p
    return False, None

# ═══════════════════════════════════════════════════════════════════════════════
# BLOC 5 : PREUVE DE POCKLINGTON
#
# Une fois un probable premier trouvé, la preuve est rapide :
#   p−1 = 2^a · 3^(b+1) · m,  F = 2^a · 3^(b+1) > √p (vrai par construction)
#   On cherche un témoin w pour q=2 et un témoin w pour q=3.
#   Typiquement trouvés en < 20 essais — quasi-instantané même pour 1M chiffres.
# ═══════════════════════════════════════════════════════════════════════════════

def pocklington(p: gmpy2.mpz, a: int, b: int) -> tuple:
    """
    Théorème de Pocklington :
    Pour chaque facteur premier q de F={2,3} :
      ∃w : w^(p-1) ≡ 1 (mod p)  et  gcd(w^((p-1)/q) − 1, p) = 1
    → p est premier.
    Retourne (True, {2: w2, 3: w3}) ou (False, état_partiel).
    """
    p1     = p - 1
    valide = {2: None, 3: None}

    for w in range(2, 2 + TEMOINS_MAX):
        gw = gmpy2.mpz(w)
        if gmpy2.powmod(gw, p1, p) != 1:
            continue
        for q in [2, 3]:
            if valide[q] is not None:
                continue
            val = gmpy2.powmod(gw, p1 // gmpy2.mpz(q), p)
            if gmpy2.gcd(val - 1, p) == 1:
                valide[q] = w
        if all(v is not None for v in valide.values()):
            return True, valide

    return False, valide

# ═══════════════════════════════════════════════════════════════════════════════
# BLOC 6 : BOUCLE DE RECHERCHE D'UN WORKER
#
# Chaque worker (rang MPI ou process) reçoit une plage [a_debut, a_fin)
# et cherche de façon autonome.  Aucune communication pendant la recherche
# (modèle "embarrassingly parallel").
# ═══════════════════════════════════════════════════════════════════════════════

def chercher(a_debut: int, a_fin: int, digits: int,
             worker_id: int, result_queue=None) -> Optional[dict]:
    """
    Recherche dans la plage a ∈ [a_debut, a_fin).
    result_queue : multiprocessing.Queue (mode local) ou None (mode MPI).
    Retourne un dict résultat ou None.
    """
    log_path = f"log_worker_{worker_id:05d}.txt"
    n_crible = 0
    n_mr     = 0
    n_paires = 0
    t0       = time.perf_counter()

    with open(log_path, "w", buffering=1) as log:
        log.write(f"Worker {worker_id} : a ∈ [{a_debut}, {a_fin})\n")
        log.write(f"Cible : {digits} chiffres (±{TOLERANCE})\n\n")

        for a in range(a_debut, a_fin):
            for b in b_approxime(a, digits):
                n_paires += 1

                # ── Étape 1 : Crible modulaire (µs) ──────────────────────
                if not crible_modulaire(a, b):
                    n_crible += 1
                    continue

                # ── Checkpoint ───────────────────────────────────────────
                if n_paires % CHECKPOINT_ITER == 0:
                    elapsed = time.perf_counter() - t0
                    log.write(
                        f"  [t={elapsed:8.1f}s] paires={n_paires} "
                        f"crible={n_crible} MR={n_mr}\n"
                    )

                # ── Étape 2 : Miller-Rabin (coûteux) ─────────────────────
                n_mr += 1
                elapsed = time.perf_counter() - t0
                log.write(
                    f"  MR #{n_mr:4d}  a={a}, b={b}  t={elapsed:.1f}s\n"
                )

                probable, p = miller_rabin_complet(a, b)
                if not probable:
                    continue

                # ── Étape 3 : Pocklington ─────────────────────────────────
                log.write(f"  *** PROBABLE PREMIER  a={a}, b={b} ***\n")
                prouve, temoins = pocklington(p, a, b)

                if prouve:
                    nb = int(gmpy2.num_digits(p))
                    log.write(f"  ✅ PROUVÉ — {nb} chiffres\n")
                    res = {
                        "a": a, "b": b, "p": p, "nb": nb,
                        "temoins": temoins,
                        "worker": worker_id,
                        "n_mr": n_mr, "n_crible": n_crible,
                        "temps": time.perf_counter() - t0
                    }
                    if result_queue is not None:
                        result_queue.put(res)
                    return res

    return None

# ═══════════════════════════════════════════════════════════════════════════════
# BLOC 7 : ORCHESTRATION
# ═══════════════════════════════════════════════════════════════════════════════

def partitionner(digits: int, n_workers: int) -> list:
    """
    Divise [1, a_max] en n_workers tranches équilibrées.
    a_max : valeur maximale de a pour laquelle p peut avoir ~digits chiffres.
    """
    # p ≈ 3·(2^a)^2 quand b=1 → log10(p) ≈ 2a·log10(2)
    # → a_max ≈ digits / (2·log10(2))
    a_max  = int(digits / (2 * _LOG10_2)) + 50
    chunk  = max(1, a_max // n_workers)
    plages = []
    for i in range(n_workers):
        debut = i * chunk + 1
        fin   = debut + chunk if i < n_workers - 1 else a_max + 1
        plages.append((debut, fin))
    return plages

def run_local(digits: int, n_workers: int):
    """Mode multiprocessing pour tests sur machine locale."""
    print(f"\n[local] {n_workers} workers — cible {digits} chiffres")
    plages = partitionner(digits, n_workers)
    queue  = mp.Queue()

    procs = []
    for wid, (a0, a1) in enumerate(plages):
        p = mp.Process(
            target=chercher,
            args=(a0, a1, digits, wid, queue)
        )
        p.start()
        procs.append(p)
        print(f"  Worker {wid} lancé : a ∈ [{a0}, {a1})")

    # Attendre le premier résultat
    result = queue.get()

    # Arrêter tous les workers
    for p in procs:
        p.terminate()
    for p in procs:
        p.join()

    return result

def run_mpi(digits: int):
    """Mode MPI pour supercalculateur."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    plages = partitionner(digits, size)
    a0, a1 = plages[rank]

    if rank == 0:
        print(f"[MPI] {size} rangs — cible {digits} chiffres")
        print(f"[MPI] Domaine total a ∈ [1, {plages[-1][1]}]")

    # Chaque rang cherche dans sa plage
    local_result = chercher(a0, a1, digits, rank)

    # Collecte sur rang 0
    all_results = comm.gather(local_result, root=0)

    if rank == 0:
        for res in all_results:
            if res is not None:
                return res
    return None

# ═══════════════════════════════════════════════════════════════════════════════
# BLOC 8 : SAUVEGARDE DU RÉSULTAT
# ═══════════════════════════════════════════════════════════════════════════════

def sauvegarder(res: dict, t_total: float):
    p  = res["p"]
    nb = res["nb"]
    a, b   = res["a"], res["b"]
    t      = res["temoins"]
    F      = deux**a * trois**(b + 1)
    m      = deux**a * trois**b - 1
    sp     = str(p)

    print(f"\n{'═'*70}")
    print(f"  PREMIER DE {nb} CHIFFRES — PRIMALITÉ PROUVÉE")
    print(f"{'═'*70}")
    print(f"  Formule : p = 3m(m+1)+1,  m = 2^{a}·3^{b}−1")
    print(f"  F = 2^{a}·3^{{{b+1}}}  ({gmpy2.num_digits(F)} chiffres) > √p  ✅")
    print(f"  Témoins Pocklington :")
    print(f"    q=2  → w={t[2]}  ✅")
    print(f"    q=3  → w={t[3]}  ✅")
    print(f"  Worker   : {res['worker']}")
    print(f"  Temps    : {t_total:.1f}s  ({t_total/3600:.2f}h)")
    print(f"  Crible   : {res['n_crible']} candidats éliminés")
    print(f"  MR tests : {res['n_mr']}")
    print(f"\n  p = {sp[:60]}…")

    nom = f"premier_{nb}ch_HPC.txt"
    with open(nom, "w") as f:
        f.write(f"Premier de {nb} chiffres — PrimeQuest-HPC\n")
        f.write(f"Méthode  : Pocklington, p=3m(m+1)+1\n")
        f.write(f"a={a}, b={b}\n")
        f.write(f"Témoins  : q=2→w={t[2]}, q=3→w={t[3]}\n")
        f.write(f"Temps    : {t_total:.1f}s\n\n")
        f.write(f"m =\n{m}\n\np =\n{p}\n")
    print(f"  Sauvegardé → {nom}")

# ═══════════════════════════════════════════════════════════════════════════════
# POINT D'ENTRÉE
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PrimeQuest-HPC")
    parser.add_argument("--digits",   type=int, default=DIGITS_CIBLE,
                        help="Nombre de chiffres cible")
    parser.add_argument("--mode",     choices=["mpi", "local"], default="mpi",
                        help="Mode d'exécution")
    parser.add_argument("--workers",  type=int, default=mp.cpu_count(),
                        help="Nombre de workers (mode local uniquement)")
    args = parser.parse_args()

    t0 = time.perf_counter()

    if args.mode == "mpi" and MPI_OK:
        res = run_mpi(args.digits)
    else:
        if args.mode == "mpi" and not MPI_OK:
            print("[warn] mpi4py non disponible — passage en mode local")
        res = run_local(args.digits, args.workers)

    t_total = time.perf_counter() - t0

    # Seul le rang 0 (ou le process principal) affiche le résultat
    rank = MPI.COMM_WORLD.Get_rank() if MPI_OK and args.mode == "mpi" else 0
    if rank == 0:
        if res:
            sauvegarder(res, t_total)
        else:
            print(f"\n  Aucun premier trouvé dans le domaine exploré.")
            print(f"  Augmenter TOLERANCE ou le nombre de workers.")
            sys.exit(1)
