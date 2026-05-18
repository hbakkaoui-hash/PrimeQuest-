#!/usr/bin/env python3
# =============================================================================
# Auteur   : Hassane BAKKAOUI — Chercheur indépendant
# Contact  : bakkahassa@hotmail.com
# Dépôt    : https://github.com/hbakkaoui-hash/PrimeQuest-
# Licence  : MIT License — voir LICENSE dans le dépôt
# Date     : Mai 2026
# Version  : 1.0.0
#
# Références (Papers I–IV) :
#   Paper I  — "Une famille paramétrique de nombres premiers : distribution
#               heuristique du décalage minimal et validation numérique"
#   Paper II — "Familles paramétriques et progressions arithmétiques :
#               densités de Bateman-Horn pour k ∈ {3,15,21,23}"
#   Paper III— "Connexion GUE / zéros de Riemann via la décomposition
#               Q(p_N) ≈ α√p + β·D_N et le terme de Weyl fractionnaire"
#   Paper IV — "Certificats de primalité de Pocklington pour la famille
#               p = k·m(m+1) + ε + 2k·q jusqu'à 50 000 chiffres"
#   arXiv    : math.NT [à compléter après soumission]
# =============================================================================
"""
generate_database.py
====================
Générateur de la base de données unifiée de la famille paramétrique de premiers.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FORMULE CENTRALE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    p = k · m·(m+1) + ε + 2k · q

où :
  p  — nombre premier étudié
  k  — paramètre de famille ∈ {3, 15, 21, 23}
       • k=3  : 2k=6,  φ(6)=2  → couvre TOUS les premiers > 3 (p ≡ ±1 mod 6)
       • k=15 : 2k=30, φ(30)=8 → premiers ≡ ±1 (mod 30)
       • k=21 : 2k=42, φ(42)=12→ premiers ≡ ±1 (mod 42)
       • k=23 : 2k=46, φ(46)=22→ premiers ≡ ±1 (mod 46)
  ε  — signe ∈ {+1, −1}, déterminé uniquement par p mod 2k
       • p ≡ 1  (mod 2k) → ε = +1
       • p ≡ −1 (mod 2k) → ε = −1
       • autrement         → p n'appartient pas à la k-famille
  m  — indice de couche (entier ≥ 0), choisi canoniquement (voir ci-dessous)
  q  — décalage signé (entier ∈ ℤ), = q_min par construction canonique

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CONVENTION CANONIQUE POUR (m, q) = q_min
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Pour un premier p et un paramètre k donnés :

  1.  ε est fixé par p mod 2k (unique valeur qui rend q entier).

  2.  Poser N = (p − ε) / k  [entier exact, toujours pair].
      Relation : N = m·(m+1) + 2q

  3.  Calculer m_low = ⌊(−1 + √(1 + 4N)) / 2⌋
      → le plus grand m tel que m·(m+1) ≤ N.

  4.  En déduire :
        q_low  = (N − m_low·(m_low+1)) / 2  ≥ 0
        q_high = q_low − (m_low + 1)          < 0
      [q_high correspond au choix m_high = m_low + 1]

  5.  Choix canonique minimisant |q| :
        |q_high| < |q_low|  →  m = m_low+1,  q = q_high  (q négatif)
        |q_high| ≥ |q_low|  →  m = m_low,    q = q_low   (q positif ; priorité
                                                            en cas d'égalité)

  →  q_min ≡ q  (même valeur, même convention dans tous les articles)
  →  Borne inconditionnelle : |q| ≤ m  (Proposition 2.3, Paper I)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PARAMÈTRES CALCULÉS — 23 COLONNES PAR LIGNE (une ligne par premier par k)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

GROUPE A — Paramètres fondamentaux de la formule
  p      (int64)  Le nombre premier.
  k      (int64)  Paramètre de famille ∈ {3, 15, 21, 23}.
  m      (int64)  Indice de couche canonique (voir convention ci-dessus).
  eps    (int64)  Signe ε ∈ {+1, −1}.
  q      (int64)  Décalage signé = q_min (convention canonique).
  q_min  (int64)  Identique à q (rappel explicite pour les articles).
  H_m    (int64)  Base hexagonale généralisée = k·m·(m+1).
                  Interprétation : "centre" de la couche m.
                  Vérification : p = H_m + ε + 2k·q  (exacte).

GROUPE B — Structure de couche
  r      (int64)  Rang de p dans la k-famille, ordre croissant de p, 1-indexé.
  n_m    (int64)  Nombre de premiers dans la couche m pour ce k
                  (cardinal de {p ∈ k-famille : couche canonique = m}).

GROUPE C — Fonctionnelles cumulatives (ordre croissant de p dans la k-famille)
  Q_r    (float64) Q_r = Σᵢ₌₁ʳ qᵢ
                   Somme cumulative des décalages signés.
                   Loi des grands nombres : Q_r / r → 0 (centrage asymptotique).

  X_r    (float64) X_r = Σᵢ₌₁ʳ (qᵢ + εᵢ) = Σᵢ₌₁ʳ cᵢ
                   Fonctionnelle Q du papier (notation Paper III).
                   Décomposition : Q_r ≈ α·√p_r + β·D_r  (Théorème 5.1).

  D_r    (float64) D_r = Σᵢ₌₁ʳ ({ √(pᵢ′/k) } − 1/2)  avec pᵢ′ = pᵢ − εᵢ
                   Terme de Weyl / équidistribution.
                   { x } = partie fractionnaire de x.
                   Teste si √(p/k) mod 1 est équidistribué dans [0,1)
                   (loi de Weyl pour les polynômes, Paper III Prop. 5).
                   Corrélation r(X_r, D_r) ≥ 0.84 (Table 6, Paper III).

GROUPE D — Distribution modulaire
  q_mod_2  (int64)  q mod 2  — parité du décalage
  q_mod_3  (int64)  q mod 3
  q_mod_5  (int64)  q mod 5
  q_mod_7  (int64)  q mod 7
  q_mod_11 (int64)  q mod 11
  q_mod_13 (int64)  q mod 13
  → Test d'équidistribution de q dans les classes résiduelles (Paper II §4).

GROUPE E — Grandeurs analytiques par couche
  ln_p    (float64) log(p)  — logarithme naturel.
  sqrt_p  (float64) √p.
  ln_m    (float64) log(m) si m > 0, sinon 0.
                    Loi logarithmique : E[|q_min|(m)] ≈ ρ_k · log(m)
                    avec ρ_k = pente × C_k (Observation 5.1, Paper I).
  nm_hat  (float64) Espérance théorique du nombre de premiers dans la couche m :
                      n̂_m = 4k(m+1) / (φ(2k) · log(k·m·(m+1)))
                    Dérivée de la densité de Dirichlet dans la progression
                    ≡ ±1 (mod 2k) sur l'intervalle [H_m − 1, H_{m+1} − 1].
  T_m     (float64) Statistique de Poisson normalisée :
                      T_m = (n_m − n̂_m) / √n̂_m
                    Sous H₀ (indépendance) : T_m ∼ N(0,1).
                    Surdispersion si T_m >> 1 sur plusieurs couches.
                    ⚠ La dernière couche (m_max) est tronquée → T_m non fiable.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
VALIDATIONS INTÉGRÉES (par chunk, loggées dans generation.log)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  1. Reconstruction exacte   : p = H_m + ε + 2k·q  pour chaque ligne
  2. Borne inconditionnelle  : |q| ≤ m  (Prop. 2.3, Paper I)
  3. Annulation arithmétique : 2k·|Σq| / Σp < 5 %
     → reflète p ≈ H_m pour grands m ; précision croissante avec N
     → 99.9964 % à N = 10^4 (Table 8, Paper III)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ARCHITECTURE TECHNIQUE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Phase 1 — Crible segmenté numpy (segments 2^23 ≈ 8 Mo)
            → Parquet Snappy brut par k (colonnes Groupes A, D, E partiels)
            → multiprocessing.Pool (4 workers, un par k) sur compute_params
            → checkpoint automatique toutes les 10^8 nombres (--resume)

  Phase 2 — Enrichissement en streaming (500 k lignes/batch)
            → comptage n_m par couche (passe rapide sur raw Parquet)
            → cumulatifs Q_r, X_r, D_r avec état courant inter-batches
            → Parquet final partitionné : FINAL_DIR/k={k}/part_XXXXXX.parquet

  Phase 3 — Export CSV gzip (10^5 lignes/fichier) + SQLite léger
            → tables : primes_parametric_sample, layers, modular_stats

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
USAGE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    python generate_database.py                      # N = 10^12 (défaut)
    python generate_database.py --n 1e9              # borne personnalisée
    python generate_database.py --n 1e9 --workers 4  # 4 workers parallèles
    python generate_database.py --resume             # reprise checkpoint
    python generate_database.py --quick              # N = 10^6, test rapide

Dépendances : numpy, pandas, pyarrow, tqdm
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import multiprocessing as mp
import os
import queue
import sqlite3
import tempfile
import threading
import time
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

K_VALUES: list[int] = [3, 15, 21, 23]

# φ(2k) = φ(6)=2, φ(30)=8, φ(42)=12, φ(46)=22
PHI_2K: Dict[int, int] = {3: 2, 15: 8, 21: 12, 23: 22}

N_MAX: int = 10**12         # borne supérieure par défaut (extensible)
SEGMENT_SIZE: int = 2**23   # ~8 Mo de booléens — meilleur débit pour grands N

# Nombre de lignes par fichier Parquet final (par k)
PARQUET_ROWS: int = 5 * 10**6

# Nombre de premiers par fichier CSV compressé (toutes familles confondues)
CSV_PRIMES: int = 10**5

# Moduli pour le Groupe D
L_MODULI: list[int] = [2, 3, 5, 7, 11, 13]

OUTPUT_DIR = Path("prime_database")
RAW_DIR    = OUTPUT_DIR / "raw"
FINAL_DIR  = OUTPUT_DIR / "parquet"
CSV_DIR    = OUTPUT_DIR / "csv"
DB_FILE    = OUTPUT_DIR / "primes_parametric.sqlite"
LOG_FILE   = OUTPUT_DIR / "generation.log"
CKPT_FILE  = OUTPUT_DIR / "checkpoint.json"

# ─────────────────────────────────────────────────────────────────────────────
# SCHÉMA PYARROW
# ─────────────────────────────────────────────────────────────────────────────

_I = pa.int64()
_F = pa.float64()

# Colonnes produites en Phase 1 (brut)
_RAW_FIELDS = [
    pa.field("p",      _I),
    pa.field("k",      _I),
    pa.field("m",      _I),
    pa.field("eps",    _I),
    pa.field("q",      _I),
    pa.field("q_min",  _I),
    pa.field("H_m",    _I),
    pa.field("ln_p",   _F),
    pa.field("sqrt_p", _F),
    pa.field("ln_m",   _F),
] + [pa.field(f"q_mod_{l}", _I) for l in L_MODULI]

SCHEMA_RAW = pa.schema(_RAW_FIELDS)

# Colonnes ajoutées en Phase 2 (enrichi)
_ENRICH_FIELDS = [
    pa.field("r",      _I),
    pa.field("n_m",    _I),
    pa.field("Q_r",    _F),
    pa.field("X_r",    _F),
    pa.field("D_r",    _F),
    pa.field("nm_hat", _F),
    pa.field("T_m",    _F),
]

SCHEMA_FULL = pa.schema(_RAW_FIELDS + _ENRICH_FIELDS)
COL_ORDER   = [f.name for f in _RAW_FIELDS + _ENRICH_FIELDS]

# ─────────────────────────────────────────────────────────────────────────────
# JOURNALISATION
# ─────────────────────────────────────────────────────────────────────────────

def setup_logging() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.FileHandler(LOG_FILE),
            logging.StreamHandler(),
        ],
    )


# ─────────────────────────────────────────────────────────────────────────────
# CRIBLE SEGMENTÉ
# ─────────────────────────────────────────────────────────────────────────────

def _small_primes(limit: int) -> np.ndarray:
    """Crible d'Ératosthène pour les premiers ≤ limit."""
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.nonzero(sieve)[0]


def prime_segments(n_max: int, seg: int = SEGMENT_SIZE) -> Iterator[np.ndarray]:
    """
    Génère des tableaux numpy de premiers en ordre croissant,
    un tableau par segment [low, high], couvrant [2, n_max].

    Crible unifié : chaque segment est traité de la même manière
    (y compris le premier [2, seg]) en évitant le biais du premier segment
    qui ne renvoyait que les petits premiers jusqu'à √N.
    """
    sqrt_n = int(n_max**0.5) + 1
    small  = _small_primes(sqrt_n)   # petits premiers pour le crible segmenté

    low = 2
    while low <= n_max:
        high  = min(low + seg - 1, n_max)
        flags = np.ones(high - low + 1, dtype=bool)
        for p in small:
            # Premier multiple de p dans [low, high] (≥ 2p pour éviter p lui-même)
            start = math.ceil(low / p) * p
            if start == p:
                start += p   # sauter p (il est premier) → commencer à 2p
            if start <= high:
                flags[start - low::p] = False
        yield np.nonzero(flags)[0] + low
        low += seg


# ─────────────────────────────────────────────────────────────────────────────
# CALCUL DES PARAMÈTRES (VECTORISÉ)
# ─────────────────────────────────────────────────────────────────────────────

def compute_params(primes: np.ndarray, k: int) -> Optional[pd.DataFrame]:
    """
    Pour un tableau de premiers et un paramètre k, calcule tous les champs
    bruts (Groupes A, D, E) pour les premiers de la k-famille.

    Convention canonique q_min :
      – ε est déterminé par p mod 2k (seule valeur possible qui rende q entier)
      – On calcule q_low ≥ 0 et q_high = q_low − (m_low+1) < 0
      – On retient le m qui minimise |q| ; en cas d'égalité, q positif (m_low)

    Retourne None si aucun premier de ce segment n'appartient à la k-famille.
    """
    two_k = 2 * k
    res   = primes % two_k

    mask_pos = res == 1            # ε = +1
    mask_neg = res == two_k - 1   # ε = −1
    mask     = mask_pos | mask_neg
    if not mask.any():
        return None

    p   = primes[mask].astype(np.int64)
    eps = np.where(mask_pos[mask], np.int64(1), np.int64(-1))

    # N = (p − ε) / k  [division entière exacte, N toujours pair]
    N = (p - eps) // np.int64(k)

    # m_low = ⌊(−1 + √(1+4N)) / 2⌋ — le plus grand m avec m(m+1) ≤ N
    m_low = (0.5 * (-1.0 + np.sqrt(1.0 + 4.0 * N.astype(np.float64)))).astype(np.int64)

    # Correction flottante (rare pour très grands N)
    m_low = np.where((m_low + 1) * (m_low + 2) <= N, m_low + 1, m_low)
    m_low = np.where(m_low * (m_low + 1) > N,        m_low - 1, m_low)

    # q_low ≥ 0,  q_high = q_low − (m_low+1) < 0
    q_low  = (N - m_low * (m_low + 1)) // 2
    q_high = q_low - (m_low + 1)

    # Choix canonique : |q_high| < |q_low| → on remonte d'une couche
    use_high = np.abs(q_high) < q_low   # q_low = |q_low| car q_low ≥ 0
    m  = np.where(use_high, m_low + 1, m_low).astype(np.int64)
    q  = np.where(use_high, q_high,    q_low).astype(np.int64)

    H_m  = (k * m * (m + 1)).astype(np.int64)
    p_f  = p.astype(np.float64)
    m_f  = m.astype(np.float64)

    data: dict = {
        "p":      p,
        "k":      np.full(len(p), k,  dtype=np.int64),
        "m":      m,
        "eps":    eps,
        "q":      q,
        "q_min":  q,                    # identique par construction canonique
        "H_m":    H_m,
        "ln_p":   np.log(p_f),
        "sqrt_p": np.sqrt(p_f),
        "ln_m":   np.where(m_f > 0.0, np.log(m_f), 0.0),
    }
    for l in L_MODULI:
        data[f"q_mod_{l}"] = (q % l).astype(np.int64)

    return pd.DataFrame(data)


# ─────────────────────────────────────────────────────────────────────────────
# VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

def validate_raw_chunk(df: pd.DataFrame, k: int) -> bool:
    """
    Trois vérifications systématiques sur un chunk brut :
      1. Reconstruction : p = H_m + ε + 2k·q
      2. Borne inconditionnelle : |q| ≤ m
      3. Annulation arithmétique : 2k·|Σq| / Σp < seuil
    Retourne True si tout est valide.
    """
    two_k = 2 * k
    ok    = True

    # 1. Reconstruction
    recon = df["H_m"] + df["eps"] + two_k * df["q"]
    bad   = int((recon != df["p"]).sum())
    if bad:
        logging.error("k=%d RECONSTRUCTION ÉCHOUÉE — %d lignes", k, bad)
        ok = False

    # 2. Borne |q| ≤ m
    bad2 = int((df["q"].abs() > df["m"]).sum())
    if bad2:
        logging.error("k=%d: %d lignes violent |q| ≤ m", k, bad2)
        ok = False

    # 3. Annulation arithmétique  2k·|Q_r| / Σp < 5 %
    if len(df) >= 100:
        sum_p  = float(df["p"].sum())
        cumq   = float(df["q"].sum())
        ratio  = two_k * abs(cumq) / sum_p if sum_p > 0 else 0.0
        if ratio > 0.05:
            logging.warning(
                "k=%d: annulation 2k·|Σq|/Σp = %.4f%% (seuil 5%%)",
                k, ratio * 100,
            )
    return ok


# ─────────────────────────────────────────────────────────────────────────────
# PHASE 1 — CRIBLE → PARQUET BRUT  (avec checkpoint par segment)
# ─────────────────────────────────────────────────────────────────────────────

def _compute_k_worker(args: Tuple[np.ndarray, int]) -> Optional[bytes]:
    """
    Worker multiprocessing : calcule les paramètres pour (segment, k) et
    sérialise en bytes Arrow IPC pour retour sans pickle numpy.
    Retourne None si aucun premier de la famille.
    """
    seg_primes, k = args
    df = compute_params(seg_primes, k)
    if df is None or df.empty:
        return None
    tbl = pa.Table.from_pandas(df, schema=SCHEMA_RAW, preserve_index=False)
    buf = pa.BufferOutputStream()
    with pa.ipc.new_stream(buf, tbl.schema) as writer:
        writer.write_table(tbl)
    return buf.getvalue().to_pybytes()


def _bytes_to_df(data: bytes) -> pd.DataFrame:
    reader = pa.ipc.open_stream(data)
    return reader.read_all().to_pandas()


def phase1_generate(
    n_max:       int,
    resume_from: int = 0,
    n_workers:   int = 1,
) -> Dict[int, int]:
    """
    Crible les premiers jusqu'à n_max, calcule les colonnes brutes et écrit
    un fichier Parquet par k dans RAW_DIR.

    – resume_from : dernier p traité (pour reprise de checkpoint).
    – n_workers   : nombre de workers parallèles pour le calcul des params
                    (un worker par k, max 4 utile). Défaut = min(4, cpu_count).

    Checkpoint automatique tous les CKPT_INTERVAL segments dans CKPT_FILE.
    Retourne le dictionnaire {k: nb_lignes}.
    """
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    writers:     Dict[int, pq.ParquetWriter] = {}
    row_counts:  Dict[int, int]              = {}
    ckpt_segs:   int                         = max(1, 100_000_000 // SEGMENT_SIZE)

    for k in K_VALUES:
        path = RAW_DIR / f"raw_k{k}.parquet"
        if path.exists() and resume_from == 0:
            path.unlink()
        writers[k]    = pq.ParquetWriter(path, SCHEMA_RAW, compression="snappy")
        row_counts[k] = 0

    total_segs = math.ceil(n_max / SEGMENT_SIZE)
    pool: Optional[mp.Pool] = None
    if n_workers > 1:
        pool = mp.Pool(processes=min(n_workers, len(K_VALUES)))

    pbar = tqdm(
        enumerate(prime_segments(n_max, seg=SEGMENT_SIZE)),
        total=total_segs,
        desc="Phase 1 — crible",
        unit="seg",
        ncols=100,
    )

    try:
        for seg_idx, seg_primes in pbar:
            if len(seg_primes) == 0:
                continue
            if seg_primes[-1] <= resume_from:
                continue
            if seg_primes[0] <= resume_from:
                seg_primes = seg_primes[seg_primes > resume_from]
            if len(seg_primes) == 0:
                continue

            if pool is not None:
                # Calcul parallèle sur les 4 valeurs de k
                results = pool.map(_compute_k_worker,
                                   [(seg_primes, k) for k in K_VALUES])
                for k, data in zip(K_VALUES, results):
                    if data is None:
                        continue
                    df = _bytes_to_df(data)
                    validate_raw_chunk(df, k)
                    tbl = pa.Table.from_pandas(df, schema=SCHEMA_RAW, preserve_index=False)
                    writers[k].write_table(tbl)
                    row_counts[k] += len(df)
            else:
                for k in K_VALUES:
                    df = compute_params(seg_primes, k)
                    if df is None or df.empty:
                        continue
                    validate_raw_chunk(df, k)
                    tbl = pa.Table.from_pandas(df, schema=SCHEMA_RAW, preserve_index=False)
                    writers[k].write_table(tbl)
                    row_counts[k] += len(df)

            # Checkpoint toutes les ckpt_segs segments
            if (seg_idx + 1) % ckpt_segs == 0:
                last_p = int(seg_primes[-1]) if len(seg_primes) else resume_from
                _update_ckpt_p1(last_p, row_counts)

            pbar.set_postfix({
                f"k{k}": f"{row_counts[k] / 1e6:.2f}M" for k in K_VALUES
            })
    finally:
        if pool:
            pool.close()
            pool.join()
        for w in writers.values():
            w.close()

    logging.info("Phase 1 terminée. Lignes : %s", row_counts)
    return row_counts


def _update_ckpt_p1(last_p: int, row_counts: Dict[int, int]) -> None:
    """Met à jour le champ phase1 du checkpoint sans écraser les autres clés."""
    ckpt = {}
    if CKPT_FILE.exists():
        try:
            with open(CKPT_FILE) as f:
                ckpt = json.load(f)
        except Exception:
            pass
    ckpt["phase1_last_p"]    = last_p
    ckpt["phase1_row_counts"] = row_counts
    with open(CKPT_FILE, "w") as f:
        json.dump(ckpt, f, indent=2)


# ─────────────────────────────────────────────────────────────────────────────
# ESPÉRANCE THÉORIQUE n̂_m
# ─────────────────────────────────────────────────────────────────────────────

def nm_hat(k: int, m_arr: np.ndarray) -> np.ndarray:
    """
    n̂_m = 4k(m+1) / (φ(2k) · log(k·m·(m+1)))

    Densité de Dirichlet dans la progression ≡ ±1 (mod 2k), multipliée par
    le nombre de candidats 2(m+1) dans la couche m.
    Garde un plancher log = 1 pour m = 0.
    """
    phi  = PHI_2K[k]
    H    = (k * m_arr * (m_arr + 1)).astype(np.float64)
    logH = np.where(H > math.e, np.log(H), 1.0)
    return 4.0 * k * (m_arr + 1).astype(np.float64) / (phi * logH)


# ─────────────────────────────────────────────────────────────────────────────
# PHASE 2 — ENRICHISSEMENT PAR k
# ─────────────────────────────────────────────────────────────────────────────

def _count_layers(k: int) -> Dict[int, int]:
    """Passe rapide sur raw_k{k}.parquet pour compter n_m par couche."""
    pf    = pq.ParquetFile(RAW_DIR / f"raw_k{k}.parquet")
    counts: Dict[int, int] = {}
    for batch in pf.iter_batches(batch_size=500_000, columns=["m"]):
        for mv in batch.column("m").to_pylist():
            counts[mv] = counts.get(mv, 0) + 1
    return counts


def phase2_enrich(k: int) -> None:
    """
    Lit raw_k{k}.parquet (par chunks pour économiser la RAM), calcule les
    colonnes d'enrichissement et écrit les fichiers Parquet finaux dans
    FINAL_DIR/k={k}/.

    Colonnes calculées ici :
      r      : rang dans la k-famille (ordre croissant de p)
      n_m    : nombre de premiers dans la couche m
      Q_r    : Σᵢ qᵢ  (somme cumulative des décalages signés)
      X_r    : Σᵢ (qᵢ + εᵢ)  (= Σ cᵢ ; terme Q du papier)
      D_r    : Σᵢ ({√(pᵢ'/k)} − 1/2)  (terme de Weyl / équidistribution)
               avec pᵢ' = pᵢ − εᵢ
      nm_hat : espérance théorique de n_m (Bateman-Horn / Dirichlet)
      T_m    : (n_m − nm_hat) / √nm_hat  (statistique de Poisson normalisée)
    """
    logging.info("k=%d : décompte des couches …", k)
    layer_nm = _count_layers(k)

    m_keys    = np.array(sorted(layer_nm.keys()), dtype=np.int64)
    nm_hat_v  = nm_hat(k, m_keys)
    nm_hat_map: Dict[int, float] = dict(zip(m_keys.tolist(), nm_hat_v.tolist()))

    out_dir = FINAL_DIR / f"k={k}"
    out_dir.mkdir(parents=True, exist_ok=True)

    pf         = pq.ParquetFile(RAW_DIR / f"raw_k{k}.parquet")
    rank       = 0
    Q_cum      = 0.0
    X_cum      = 0.0
    D_cum      = 0.0
    file_index = 0

    total_rows = sum(layer_nm.values())
    pbar = tqdm(
        pf.iter_batches(batch_size=500_000),
        total=math.ceil(total_rows / 500_000),
        desc=f"Phase 2 — k={k}",
        unit="lot",
        ncols=100,
    )

    for batch in pbar:
        df = batch.to_pandas()
        n  = len(df)

        # ── Rang r ────────────────────────────────────────────────────────────
        df["r"] = np.arange(rank + 1, rank + n + 1, dtype=np.int64)
        rank   += n

        # ── n_m, nm_hat, T_m ──────────────────────────────────────────────────
        df["n_m"]    = df["m"].map(layer_nm).astype(np.int64)
        df["nm_hat"] = df["m"].map(nm_hat_map).astype(np.float64)
        nm_arr = df["nm_hat"].to_numpy()
        nm_obs = df["n_m"].to_numpy(dtype=np.float64)
        df["T_m"] = np.where(
            nm_arr > 0.0,
            (nm_obs - nm_arr) / np.sqrt(nm_arr),
            0.0,
        )

        # ── Cumulatifs ────────────────────────────────────────────────────────
        q_arr   = df["q"].to_numpy(dtype=np.float64)
        eps_arr = df["eps"].to_numpy(dtype=np.float64)
        p_pr    = (df["p"] - df["eps"]).to_numpy(dtype=np.float64)  # p' = p − ε

        # Q_r = Σ qᵢ
        Q_loc    = np.cumsum(q_arr)
        df["Q_r"] = Q_loc + Q_cum
        Q_cum    = float(Q_loc[-1]) + Q_cum

        # X_r = Σ (qᵢ + εᵢ) = Σ cᵢ
        X_loc    = np.cumsum(q_arr + eps_arr)
        df["X_r"] = X_loc + X_cum
        X_cum    = float(X_loc[-1]) + X_cum

        # D_r = Σ ({√(p'/k)} − 1/2)  — terme de Weyl (Prop. 5 du papier)
        frac     = np.modf(np.sqrt(p_pr / k))[0]   # partie fractionnaire
        D_loc    = np.cumsum(frac - 0.5)
        df["D_r"] = D_loc + D_cum
        D_cum    = float(D_loc[-1]) + D_cum

        # ── Ordonnancement des colonnes + écriture ────────────────────────────
        df = df[COL_ORDER]
        tbl = pa.Table.from_pandas(df, schema=SCHEMA_FULL, preserve_index=False)
        pq.write_table(
            tbl,
            out_dir / f"part_{file_index:06d}.parquet",
            compression="snappy",
        )
        file_index += 1

    logging.info("k=%d : %d fichiers Parquet écrits dans %s", k, file_index, out_dir)


# ─────────────────────────────────────────────────────────────────────────────
# PHASE 3 — EXPORT CSV + SQLITE
# ─────────────────────────────────────────────────────────────────────────────

def export_csv_k(k: int) -> None:
    """Exporte les fichiers Parquet finaux en CSV compressés (CSV_PRIMES lignes/fichier)."""
    CSV_DIR.mkdir(parents=True, exist_ok=True)
    src_files = sorted((FINAL_DIR / f"k={k}").glob("*.parquet"))
    if not src_files:
        logging.warning("k=%d : aucun fichier Parquet final trouvé — CSV ignoré", k)
        return

    buf:     list[pd.DataFrame] = []
    buf_len: int                = 0
    idx:     int                = 0

    for src in tqdm(src_files, desc=f"CSV k={k}", unit="file", ncols=80):
        chunk = pd.read_parquet(src)
        buf.append(chunk)
        buf_len += len(chunk)

        while buf_len >= CSV_PRIMES:
            merged = pd.concat(buf, ignore_index=True)
            out    = merged.iloc[:CSV_PRIMES]
            rest   = merged.iloc[CSV_PRIMES:]
            path   = CSV_DIR / f"k{k}_{idx:06d}.csv.gz"
            out.to_csv(path, index=False, compression="gzip")
            idx     += 1
            buf      = [rest] if len(rest) > 0 else []
            buf_len  = len(rest)

    if buf_len > 0:
        merged = pd.concat(buf, ignore_index=True)
        path   = CSV_DIR / f"k{k}_{idx:06d}.csv.gz"
        merged.to_csv(path, index=False, compression="gzip")

    logging.info("k=%d : %d fichiers CSV écrits", k, idx + 1)


def build_sqlite() -> None:
    """
    Construit la base SQLite légère avec trois tables :
      primes_parametric_sample : ~50 000 lignes par k (head + tail)
      layers                   : une ligne par (k, m) avec agrégats
      modular_stats            : distribution de q mod l par blocs de CSV_PRIMES
    """
    logging.info("Construction SQLite …")
    if DB_FILE.exists():
        DB_FILE.unlink()
    con = sqlite3.connect(DB_FILE)
    con.execute("PRAGMA journal_mode = WAL")
    con.execute("PRAGMA synchronous  = NORMAL")

    sample_parts: list[pd.DataFrame] = []
    layers_rows:  list[dict]         = []
    mod_rows:     list[dict]         = []

    for k in K_VALUES:
        src_files = sorted((FINAL_DIR / f"k={k}").glob("*.parquet"))
        if not src_files:
            continue

        # ── Échantillon head/tail (dédupliqué sur p) ──────────────────────────
        head_df = pd.read_parquet(src_files[0]).head(25_000)
        tail_df = pd.read_parquet(src_files[-1]).tail(25_000)
        combined = pd.concat([head_df, tail_df], ignore_index=True).drop_duplicates(subset="p")
        sample_parts.append(combined)

        # ── Table layers ──────────────────────────────────────────────────────
        for src in tqdm(src_files, desc=f"layers k={k}", unit="file", ncols=80, leave=False):
            chunk = pd.read_parquet(src, columns=["m", "p", "n_m", "nm_hat", "T_m", "Q_r", "ln_p"])
            for m_val, grp in chunk.groupby("m", sort=False):
                layers_rows.append({
                    "k":        int(k),
                    "m":        int(m_val),
                    "H_m":      int(k * int(m_val) * (int(m_val) + 1)),
                    "n_m":      int(grp["n_m"].iloc[0]),
                    "nm_hat":   float(grp["nm_hat"].iloc[0]),
                    "T_m":      float(grp["T_m"].iloc[0]),
                    "sum_q_r":  float(grp["Q_r"].iloc[-1]),
                    "avg_ln_p": float(grp["ln_p"].mean()),
                })

        # ── Table modular_stats ────────────────────────────────────────────────
        r_offset = 0
        for src in src_files:
            chunk  = pd.read_parquet(src, columns=["r"] + [f"q_mod_{l}" for l in L_MODULI])
            n_blks = math.ceil(len(chunk) / CSV_PRIMES)
            for b in range(n_blks):
                bl = chunk.iloc[b * CSV_PRIMES:(b + 1) * CSV_PRIMES]
                row: dict = {
                    "k":     int(k),
                    "block": r_offset + b,
                    "r_min": int(bl["r"].iloc[0]),
                    "r_max": int(bl["r"].iloc[-1]),
                }
                for l in L_MODULI:
                    vc = bl[f"q_mod_{l}"].value_counts().sort_index()
                    for v in range(l):
                        row[f"cnt_mod{l}_{v}"] = int(vc.get(v, 0))
                mod_rows.append(row)
            r_offset += n_blks

    # Écriture SQLite
    pd.concat(sample_parts, ignore_index=True).to_sql(
        "primes_parametric_sample", con, if_exists="replace", index=False
    )
    con.execute(
        "CREATE INDEX IF NOT EXISTS idx_pps ON primes_parametric_sample (p, k)"
    )

    # Dédupliquer layers (plusieurs chunks peuvent contribuer la même (k,m))
    layers_df = (
        pd.DataFrame(layers_rows)
        .groupby(["k", "m"], as_index=False)
        .agg({"H_m": "first", "n_m": "first", "nm_hat": "first",
              "T_m": "first", "sum_q_r": "last", "avg_ln_p": "mean"})
    )
    layers_df.to_sql("layers", con, if_exists="replace", index=False)
    con.execute("CREATE INDEX IF NOT EXISTS idx_layers ON layers (k, m)")

    pd.DataFrame(mod_rows).to_sql(
        "modular_stats", con, if_exists="replace", index=False
    )

    con.commit()
    con.close()
    logging.info("SQLite écrit : %s", DB_FILE)


# ─────────────────────────────────────────────────────────────────────────────
# CHECKPOINT
# ─────────────────────────────────────────────────────────────────────────────

def load_checkpoint() -> dict:
    if CKPT_FILE.exists():
        with open(CKPT_FILE) as f:
            return json.load(f)
    return {"phase": 0, "last_p": 0}


def save_checkpoint(state: dict) -> None:
    with open(CKPT_FILE, "w") as f:
        json.dump(state, f, indent=2)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Générateur de la base de données paramétrique")
    parser.add_argument("--n",       type=float, default=N_MAX,
                        help="Borne supérieure (défaut 10^12)")
    parser.add_argument("--quick",   action="store_true",
                        help="Mode test : N=10^6")
    parser.add_argument("--resume",  action="store_true",
                        help="Reprendre depuis le checkpoint")
    parser.add_argument("--workers", type=int,
                        default=min(4, mp.cpu_count()),
                        help="Workers parallèles pour Phase 1 (défaut = min(4, cœurs))")
    args = parser.parse_args()

    n_max    = 10**6 if args.quick else int(args.n)
    n_workers = args.workers
    setup_logging()
    logging.info("N_MAX = %d, k ∈ %s, workers = %d", n_max, K_VALUES, n_workers)

    ckpt = load_checkpoint() if args.resume else {"phase": 0, "last_p": 0}
    t0   = time.time()

    # ── Phase 1 ──────────────────────────────────────────────────────────────
    if ckpt.get("phase", 0) < 1:
        resume_p   = ckpt.get("phase1_last_p", 0)
        row_counts = phase1_generate(n_max, resume_from=resume_p, n_workers=n_workers)
        ckpt.update({"phase": 1, "last_p": n_max, "row_counts": row_counts})
        save_checkpoint(ckpt)
    else:
        logging.info("Phase 1 ignorée (checkpoint)")

    # ── Phase 2 ──────────────────────────────────────────────────────────────
    for k in K_VALUES:
        phase_key = f"phase2_k{k}_done"
        if ckpt.get(phase_key):
            logging.info("Phase 2 k=%d ignorée (checkpoint)", k)
            continue
        phase2_enrich(k)
        ckpt[phase_key] = True
        save_checkpoint(ckpt)

    # ── Phase 3 ──────────────────────────────────────────────────────────────
    for k in K_VALUES:
        export_csv_k(k)

    build_sqlite()

    elapsed = time.time() - t0
    logging.info("Terminé en %.1f s (%.1f min)", elapsed, elapsed / 60)
    save_checkpoint({**ckpt, "phase": 3, "done": True})

    print(f"\n✓ Base générée en {elapsed:.1f}s")
    print(f"  Répertoire : {OUTPUT_DIR.resolve()}")
    print(f"  SQLite     : {DB_FILE}")
    print(f"  Parquet    : {FINAL_DIR}/k={{3,15,21,23}}/")
    print(f"  CSV        : {CSV_DIR}/")


if __name__ == "__main__":
    main()
