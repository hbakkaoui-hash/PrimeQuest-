#!/usr/bin/env python3
# =============================================================================
# Auteur   : Hassane BAKKAOUI — Chercheur indépendant
# Contact  : bakkahassa@hotmail.com
# Dépôt    : https://github.com/hbakkaoui-hash/PrimeQuest-
# Licence  : MIT License — voir LICENSE dans le dépôt
# Date     : Mai 2026 — Version 1.1.0
#
# Références :
#   Paper I  — Distribution du décalage minimal q_min et validation numérique
#   Paper II — Densités de Bateman-Horn pour k ∈ {3,15,21,23}
#   Paper III— Connexion GUE / zéros de Riemann, terme de Weyl D_r
# =============================================================================
"""
primes_parametric_generator.py
================================
Générateur haute performance de la table primes_parametric pour k = 3.
Famille paramétrique : p = 3·m(m+1) + ε + 6·q  avec N = 10^8 premiers.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FORMULE CENTRALE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    p = k · m·(m+1) + ε + 2k · q     k = 3, ε ∈ {+1,−1}, m ≥ 0, q ∈ ℤ

CONVENTION CANONIQUE q_min
━━━━━━━━━━━━━━━━━━━━━━━━━━
    1. ε : p ≡ 1 (mod 6) → +1 ; p ≡ 5 (mod 6) → −1
    2. N = (p − ε) / k                    (entier pair)
    3. m_low = ⌊(−1 + √(1 + 4N)) / 2⌋
    4. q_low  = (N − m_low·(m_low+1)) / 2   ≥ 0
       q_high = q_low − (m_low + 1)          < 0
    5. |q_high| < q_low  →  m = m_low+1, q = q_high
       sinon              →  m = m_low,   q = q_low   (priorité signe+ à égalité)

    Borne inconditionnelle : |q| ≤ m  (Proposition 2.3, Paper I)

COLONNES GÉNÉRÉES
━━━━━━━━━━━━━━━━━
    r       rang dans la famille (1-indexé)
    p       nombre premier
    k       paramètre de famille (= 3)
    m       indice de couche canonique
    eps     ε ∈ {+1, −1}
    q_min   décalage signé minimal
    H_m     base hexagonale = k·m·(m+1)
    ln_p    log(p)
    sqrt_p  √p
    Q_r     Σᵢ₌₁ʳ qᵢ                        (somme cumulative des décalages)
    X_r     Σᵢ₌₁ʳ (qᵢ + εᵢ)                 (fonctionnelle Q du papier)
    D_r     Σᵢ₌₁ʳ ({√(pᵢ'/k)} − ½)          (terme de Weyl, pᵢ' = pᵢ − εᵢ)

ARCHITECTURE
━━━━━━━━━━━━
    • Crible segmenté numpy (segments de 2^23 ≈ 8 Mo)
    • Traitement vectorisé par batch (numpy, ~100× plus rapide que boucle Python)
    • Écriture Parquet partitionnée par m_block = m // 1000
    • Checkpoints JSON (reprise après interruption en quelques secondes)
    • Validation par chunk : reconstruction, borne |q| ≤ m, annulation arithmétique

Dépendances : numpy, pandas, pyarrow, tqdm
"""

import json
import logging
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

K             = 3
N_TARGET      = 10**8           # Nombre de premiers à collecter
SEGMENT_SIZE  = 2**23           # Taille du segment crible (~8 Mo de RAM)
BATCH_SIZE    = 10**6           # Premiers traités par batch numpy
OUTPUT_DIR    = Path("primes_k3")
CHECKPOINT    = Path("checkpoint_k3.json")
VIOLATIONS_LOG= Path("violations_k3.log")
COMPRESSION   = "snappy"

# ═══════════════════════════════════════════════════════════════════════════════
# LOGGING DES VIOLATIONS
# ═══════════════════════════════════════════════════════════════════════════════

logging.basicConfig(
    filename=VIOLATIONS_LOG,
    level=logging.CRITICAL,
    format="%(asctime)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def log_violation(msg: str) -> None:
    logging.critical("VIOLATION: %s", msg)
    print(f"\n[VIOLATION CRITIQUE] {msg}")
    print(f"  Détails dans : {VIOLATIONS_LOG}")
    sys.exit(1)


# ═══════════════════════════════════════════════════════════════════════════════
# CRIBLE SEGMENTÉ (numpy)
# ═══════════════════════════════════════════════════════════════════════════════

def _small_primes(limit: int) -> np.ndarray:
    """Petits premiers jusqu'à limit via crible d'Ératosthène classique."""
    if limit < 2:
        return np.array([], dtype=np.int64)
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.nonzero(sieve)[0].astype(np.int64)


def prime_segments(n_max: int, seg: int = SEGMENT_SIZE):
    """
    Générateur de segments de premiers (np.ndarray) couvrant [2, n_max].

    Pour chaque segment [low, high], retourne les premiers dans ce segment
    sous forme de tableau numpy int64.
    """
    sqrt_n  = int(n_max**0.5) + 1
    small   = _small_primes(sqrt_n)
    low     = 2

    while low <= n_max:
        high  = min(low + seg - 1, n_max)
        flags = np.ones(high - low + 1, dtype=bool)
        for p in small:
            start = math.ceil(low / p) * p
            if start == p:
                start += p
            if start <= high:
                flags[start - low::p] = False
        yield np.nonzero(flags)[0] + low
        low += seg


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENTION CANONIQUE — VECTORISÉE (numpy)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_params(primes: np.ndarray):
    """
    Calcule les paramètres canoniques (m, ε, q, H_m) pour un tableau de premiers.

    Seuls les premiers p ≡ ±1 (mod 6) sont conservés (= toute la k=3 famille).

    Retourne (p, eps, m, q, H_m) — tous np.ndarray int64.
    """
    two_k = 2 * K  # 6
    res   = primes % two_k
    mask  = (res == 1) | (res == two_k - 1)

    p   = primes[mask].astype(np.int64)
    eps = np.where(res[mask] == 1, np.int64(1), np.int64(-1))

    # N = (p − ε) / k  (entier exact, toujours pair)
    N = (p - eps) // np.int64(K)

    # m_low = ⌊(−1 + √(1 + 4N)) / 2⌋  avec correction des erreurs flottantes
    m_low = ((-1.0 + np.sqrt(1.0 + 4.0 * N.astype(np.float64))) / 2.0).astype(np.int64)
    m_low = np.where((m_low + 1) * (m_low + 2) <= N, m_low + 1, m_low)
    m_low = np.where(m_low * (m_low + 1) > N,        m_low - 1, m_low)

    # q_low et q_high
    q_low  = (N - m_low * (m_low + 1)) // 2   # ≥ 0
    q_high = q_low - (m_low + 1)               # < 0

    # Choix canonique : |q_high| < q_low → prendre m_high
    use_high = np.abs(q_high) < q_low
    m   = np.where(use_high, m_low + 1, m_low).astype(np.int64)
    q   = np.where(use_high, q_high,    q_low ).astype(np.int64)
    H_m = (np.int64(K) * m * (m + 1)).astype(np.int64)

    return p, eps, m, q, H_m


# ═══════════════════════════════════════════════════════════════════════════════
# VALIDATION PAR CHUNK
# ═══════════════════════════════════════════════════════════════════════════════

def validate_chunk(p: np.ndarray, eps: np.ndarray, m: np.ndarray,
                   q: np.ndarray, H_m: np.ndarray) -> None:
    """
    Trois vérifications sur le chunk courant (arrêt immédiat sur violation).

    1. Reconstruction : p == H_m + ε + 2k·q
    2. Borne : |q| ≤ m  (Proposition 2.3, Paper I)
    3. H_m = k·m·(m+1)
    """
    # 1. Reconstruction
    recon = H_m + eps + 2 * np.int64(K) * q
    bad   = np.where(recon != p)[0]
    if bad.size:
        i = bad[0]
        log_violation(
            f"Reconstruction échouée : p={p[i]}, H_m+ε+2kq={recon[i]}, "
            f"m={m[i]}, eps={eps[i]}, q={q[i]}"
        )

    # 2. Borne |q| ≤ m
    bad = np.where(np.abs(q) > m)[0]
    if bad.size:
        i = bad[0]
        log_violation(
            f"Borne violée : |q|={abs(q[i])} > m={m[i]} pour p={p[i]}"
        )

    # 3. H_m correct
    expected_H = np.int64(K) * m * (m + 1)
    bad = np.where(H_m != expected_H)[0]
    if bad.size:
        i = bad[0]
        log_violation(
            f"H_m incorrect : H_m={H_m[i]}, attendu={expected_H[i]} pour p={p[i]}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# CHECKPOINT (JSON)
# ═══════════════════════════════════════════════════════════════════════════════

def save_checkpoint(r: int, p_last: int, Q_r: int, X_r: int, D_r: float) -> None:
    data = {"r": r, "p_last": p_last, "Q_r": Q_r, "X_r": X_r, "D_r": D_r}
    with open(CHECKPOINT, "w") as f:
        json.dump(data, f)


def load_checkpoint() -> dict | None:
    if not CHECKPOINT.exists():
        return None
    try:
        with open(CHECKPOINT) as f:
            return json.load(f)
    except Exception:
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# ÉCRITURE PARQUET PARTITIONNÉE (par m_block = m // 1000)
# ═══════════════════════════════════════════════════════════════════════════════

SCHEMA = pa.schema([
    ("r",      pa.int64()),
    ("p",      pa.int64()),
    ("k",      pa.int16()),
    ("m",      pa.int64()),
    ("eps",    pa.int8()),
    ("q_min",  pa.int64()),
    ("H_m",    pa.int64()),
    ("ln_p",   pa.float64()),
    ("sqrt_p", pa.float64()),
    ("Q_r",    pa.int64()),
    ("X_r",    pa.int64()),
    ("D_r",    pa.float64()),
])


def flush_batch(df: pd.DataFrame, out_dir: Path) -> None:
    """Écrit le DataFrame partitionné par m_block = m // 1000."""
    for m_block, grp in df.groupby(df["m"] // 1000):
        part_dir = out_dir / f"m_block={m_block:04d}"
        part_dir.mkdir(parents=True, exist_ok=True)
        r_first = grp["r"].iloc[0]
        fname   = part_dir / f"part_{r_first:010d}.parquet"
        tbl = pa.Table.from_pandas(
            grp[list(SCHEMA.names)], schema=SCHEMA, preserve_index=False
        )
        pq.write_table(tbl, fname, compression=COMPRESSION)


# ═══════════════════════════════════════════════════════════════════════════════
# FONCTION PRINCIPALE
# ═══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    t0 = time.time()
    OUTPUT_DIR.mkdir(exist_ok=True)

    # ── Estimation de la limite du crible (N-ième premier) ──────────────────
    # Formule de Rosser-Schoenfeld + marge 15 %
    n     = float(N_TARGET)
    ln_n  = math.log(n)
    p_est = int(n * (ln_n + math.log(ln_n) - 1.0) * 1.15)

    print("=" * 70)
    print(f"  GÉNÉRATEUR primes_parametric — k={K}, N={N_TARGET:,}")
    print(f"  Limite crible estimée : {p_est:,}")
    print("=" * 70)

    # ── Reprise depuis checkpoint ────────────────────────────────────────────
    ckpt    = load_checkpoint()
    r_off   = 0
    Q_r_cum = 0
    X_r_cum = 0
    D_r_cum = 0.0
    p_resume = 0

    if ckpt:
        r_off    = ckpt["r"]
        p_resume = ckpt["p_last"]
        Q_r_cum  = ckpt["Q_r"]
        X_r_cum  = ckpt["X_r"]
        D_r_cum  = ckpt["D_r"]
        print(f"  Reprise : r={r_off:,}, p_last={p_resume:,}")

    # ── Boucle principale ────────────────────────────────────────────────────
    buffer: list[pd.DataFrame] = []
    buf_rows = 0
    r_total  = r_off

    pbar = tqdm(
        total=N_TARGET, initial=r_off,
        desc="Premiers", unit="p", ncols=80,
    )

    for seg in prime_segments(p_est):
        # Reprise : ignorer les segments entièrement déjà traités
        if p_resume and seg[-1] <= p_resume:
            continue
        if p_resume and seg[0] <= p_resume:
            seg = seg[seg > p_resume]
        if seg.size == 0:
            continue

        # Paramètres canoniques (vectorisé)
        p, eps, m, q, H_m = compute_params(seg)
        if p.size == 0:
            continue

        # Tronquer si on dépasse N_TARGET
        remaining = N_TARGET - r_total
        if p.size > remaining:
            p, eps, m, q, H_m = (
                p[:remaining], eps[:remaining], m[:remaining],
                q[:remaining], H_m[:remaining],
            )

        validate_chunk(p, eps, m, q, H_m)

        n   = p.size
        r   = np.arange(r_total + 1, r_total + n + 1, dtype=np.int64)

        # Cumulatifs locaux puis mise à jour globale
        Q_loc = np.cumsum(q)
        X_loc = np.cumsum(q + eps)
        p_prime = p - eps                                          # p' = p − ε
        frac    = np.sqrt(p_prime.astype(np.float64) / K) % 1.0   # {√(p'/k)}
        D_loc   = np.cumsum(frac - 0.5)

        df = pd.DataFrame({
            "r":      r,
            "p":      p,
            "k":      np.full(n, K, dtype=np.int16),
            "m":      m,
            "eps":    eps.astype(np.int8),
            "q_min":  q,
            "H_m":    H_m,
            "ln_p":   np.log(p.astype(np.float64)),
            "sqrt_p": np.sqrt(p.astype(np.float64)),
            "Q_r":    Q_loc + Q_r_cum,
            "X_r":    X_loc + X_r_cum,
            "D_r":    D_loc + D_r_cum,
        })

        Q_r_cum += int(Q_loc[-1])
        X_r_cum += int(X_loc[-1])
        D_r_cum += float(D_loc[-1])
        r_total += n

        buffer.append(df)
        buf_rows += n
        pbar.update(n)

        if buf_rows >= BATCH_SIZE:
            flush_batch(pd.concat(buffer, ignore_index=True), OUTPUT_DIR)
            save_checkpoint(r_total, int(p[-1]), Q_r_cum, X_r_cum, D_r_cum)
            buffer   = []
            buf_rows = 0

        if r_total >= N_TARGET:
            break

    pbar.close()

    # Flush final
    if buffer:
        flush_batch(pd.concat(buffer, ignore_index=True), OUTPUT_DIR)

    if CHECKPOINT.exists():
        CHECKPOINT.unlink()

    # ── Validation finale (streaming, sans tout charger en RAM) ─────────────
    elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print("  VALIDATION FINALE")
    print("=" * 70)

    all_files = sorted(OUTPUT_DIR.rglob("*.parquet"))
    if not all_files:
        print("  Aucun fichier Parquet trouvé.")
        return

    # 5 premières lignes
    first_df = pq.read_table(all_files[0]).to_pandas().sort_values("r").head(5)
    print("\n  5 premières lignes :")
    print(first_df[["r", "p", "m", "eps", "q_min", "Q_r", "X_r", "D_r"]].to_string(index=False))

    # 5 dernières lignes (lire uniquement le dernier fichier)
    last_df = pq.read_table(all_files[-1]).to_pandas().sort_values("r").tail(5)
    print("\n  5 dernières lignes :")
    print(last_df[["r", "p", "m", "eps", "q_min", "Q_r", "X_r", "D_r"]].to_string(index=False))

    # Statistiques globales (streaming)
    total_rows   = 0
    p_max        = 0
    m_max        = 0
    sum_q_first100 = None

    first100_rows: list[pd.DataFrame] = []
    first100_done = False

    for f in all_files:
        tbl = pq.read_table(f, columns=["r", "p", "m", "q_min"])
        df  = tbl.to_pandas()
        total_rows += len(df)
        p_max = max(p_max, int(df["p"].max()))
        m_max = max(m_max, int(df["m"].max()))
        if not first100_done:
            first100_rows.append(df[df["r"] <= 100])
            combined = pd.concat(first100_rows, ignore_index=True)
            if len(combined) >= 100:
                sum_q_first100 = int(combined[combined["r"] <= 100]["q_min"].sum())
                first100_done = True

    total_size = sum(f.stat().st_size for f in all_files)
    n_parts    = len(sorted(OUTPUT_DIR.glob("m_block=*")))

    print(f"\n  Lignes générées : {total_rows:,} / {N_TARGET:,}")
    m_max_expected = int(math.sqrt(p_max / K))
    print(f"  p_max atteint   : {p_max:,}  (attendu ≈ 2.04×10^9)")
    print(f"  m_max           : {m_max}  (attendu ≈ {m_max_expected:,}  = √(p_max/k))")
    if sum_q_first100 is not None:
        print(f"  Σ q_min [1..100]: {sum_q_first100}")
    print(f"  Partitions      : {n_parts}  (m_block)")
    print(f"  Taille Parquet  : {total_size / 1024**3:.2f} Go")
    print(f"  Durée totale    : {elapsed:.1f}s  ({elapsed/60:.1f} min)")
    print("\n" + "=" * 70)
    print("  GENERATION TERMINEE")
    print("=" * 70)


if __name__ == "__main__":
    main()
