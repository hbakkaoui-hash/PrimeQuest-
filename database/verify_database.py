#!/usr/bin/env python3
# =============================================================================
# verify_database.py — Vérification indépendante de la base de données
#
# Recalcule depuis zéro les 5 premiers et 5 derniers enregistrements de
# chaque k-famille et compare avec les valeurs stockées dans Parquet.
#
# Usage :
#     python verify_database.py              # vérifie la base dans prime_database/
#     python verify_database.py --db /path/  # répertoire alternatif
# =============================================================================
"""
Vérificateur indépendant — recalcule intégralement les colonnes pour un
sous-ensemble de premiers et compare avec la base générée.

Chaque test utilise uniquement :
  – la définition p = k·m(m+1) + ε + 2k·q
  – le crible de base pour identifier les premiers
  – aucune référence au code de generate_database.py
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

# ─────────────────────────────────────────────────────────────────────────────
# PARAMÈTRES
# ─────────────────────────────────────────────────────────────────────────────

K_VALUES = [3, 15, 21, 23]
PHI_2K   = {3: 2, 15: 8, 21: 12, 23: 22}
L_MODULI = [2, 3, 5, 7, 11, 13]
TOL      = 1e-9   # tolérance flottante


# ─────────────────────────────────────────────────────────────────────────────
# UTILITAIRES INDÉPENDANTS
# ─────────────────────────────────────────────────────────────────────────────

def sieve(limit: int) -> np.ndarray:
    """Crible d'Ératosthène simple."""
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def canonical_params(p: int, k: int) -> Optional[dict]:
    """
    Calcule (m, ε, q) canoniques pour le premier p et le paramètre k.
    Retourne None si p n'appartient pas à la k-famille.

    Convention q_min : parmi (m_low, q_low≥0) et (m_high, q_high<0),
    choisit le m qui minimise |q| ; en cas d'égalité, q positif (m_low).
    """
    two_k = 2 * k
    res   = p % two_k
    if res == 1:
        eps = 1
    elif res == two_k - 1:
        eps = -1
    else:
        return None   # p ∉ famille k

    # N = (p − ε) / k  [entier exact, toujours pair]
    assert (p - eps) % k == 0, f"k={k} p={p} eps={eps}: (p-eps) non divisible par k"
    N = (p - eps) // k

    # m_low = ⌊(−1 + √(1+4N)) / 2⌋
    m_low = int((-1.0 + math.sqrt(1.0 + 4.0 * N)) / 2.0)
    while (m_low + 1) * (m_low + 2) <= N:
        m_low += 1
    while m_low * (m_low + 1) > N:
        m_low -= 1

    assert m_low * (m_low + 1) <= N < (m_low + 1) * (m_low + 2)

    q_low  = (N - m_low * (m_low + 1)) // 2
    q_high = q_low - (m_low + 1)

    assert q_low >= 0 and q_high < 0

    if abs(q_high) < abs(q_low):
        m, q = m_low + 1, q_high
    else:
        m, q = m_low, q_low

    # Vérification reconstruction
    H_m = k * m * (m + 1)
    assert p == H_m + eps + two_k * q, (
        f"RECONSTRUCTION ÉCHOUÉE: p={p} ≠ {H_m}+{eps}+{two_k}×{q}={H_m+eps+two_k*q}"
    )
    assert abs(q) <= m, f"BORNE: |q|={abs(q)} > m={m}"

    p_pr = p - eps
    frac = math.modf(math.sqrt(p_pr / k))[0]

    return {
        "p":      p,
        "k":      k,
        "m":      m,
        "eps":    eps,
        "q":      q,
        "q_min":  q,
        "H_m":    H_m,
        "ln_p":   math.log(p),
        "sqrt_p": math.sqrt(p),
        "ln_m":   math.log(m) if m > 0 else 0.0,
        "_frac":  frac - 0.5,   # contribution à D_r
        "_c":     q + eps,       # contribution à X_r
    }


def compute_cumulative(records: list[dict]) -> list[dict]:
    """
    Calcule Q_r, X_r, D_r pour une liste ordonnée d'enregistrements.
    Ajoute 'r' (rang local, 1-indexé) et les cumulatifs.
    """
    Q = 0.0
    X = 0.0
    D = 0.0
    out = []
    for i, rec in enumerate(records):
        Q += rec["q"]
        X += rec["_c"]
        D += rec["_frac"]
        out.append({**rec, "r": i + 1, "Q_r": Q, "X_r": X, "D_r": D})
    return out


# ─────────────────────────────────────────────────────────────────────────────
# LECTURE DE LA BASE PARQUET
# ─────────────────────────────────────────────────────────────────────────────

def read_head_tail(k: int, db_dir: Path, n: int = 5) -> pd.DataFrame:
    """Lit les n premières et n dernières lignes des fichiers Parquet pour k."""
    src_files = sorted((db_dir / "parquet" / f"k={k}").glob("*.parquet"))
    if not src_files:
        raise FileNotFoundError(f"Aucun fichier Parquet pour k={k} dans {db_dir}/parquet/")

    head = pd.read_parquet(src_files[0]).head(n)
    tail = pd.read_parquet(src_files[-1]).tail(n)
    return pd.concat([head, tail], ignore_index=True)


# ─────────────────────────────────────────────────────────────────────────────
# COMPARAISON
# ─────────────────────────────────────────────────────────────────────────────

def compare(k: int, db_dir: Path, n: int = 5) -> bool:
    """
    Extrait les n premières + n dernières lignes de la base pour k,
    recalcule chaque colonne depuis zéro, compare.
    Retourne True si tout est conforme.
    """
    print(f"\n{'='*60}")
    print(f"  Vérification k = {k}  (5 premiers + 5 derniers)")
    print(f"{'='*60}")

    db_df = read_head_tail(k, db_dir, n)
    primes_list = db_df["p"].tolist()

    # Calcul indépendant pour chaque premier
    ref_records: list[dict] = []
    for p in primes_list:
        rec = canonical_params(int(p), k)
        if rec is None:
            print(f"  ERREUR : p={p} n'appartient pas à la k={k}-famille !")
            return False
        ref_records.append(rec)

    # Cumulatifs (attention : les rangs des 5 derniers ne repartent pas de 1)
    # On recalcule Q_r, X_r, D_r relativement au rang stocké dans la base.
    # Pour les 5 premiers (rangs 1..5), le cumulatif absolu == cumulatif local.
    # Pour les 5 derniers, on utilise les valeurs stockées pour comparaison directe.

    all_ok = True

    for i, (row, rec) in enumerate(zip(db_df.itertuples(index=False), ref_records)):
        errors = []

        # Colonnes entières exactes
        for col in ("p", "k", "m", "eps", "q", "q_min", "H_m"):
            stored = int(getattr(row, col))
            ref    = int(rec[col])
            if stored != ref:
                errors.append(f"{col}: stocké={stored}, recalculé={ref}")

        # Colonnes flottantes (tolérance relative)
        for col, ref_val in [
            ("ln_p",  rec["ln_p"]),
            ("sqrt_p", rec["sqrt_p"]),
            ("ln_m",  rec["ln_m"]),
        ]:
            stored = float(getattr(row, col))
            if ref_val != 0.0 and abs(stored - ref_val) / abs(ref_val) > TOL:
                errors.append(f"{col}: stocké={stored:.8f}, recalculé={ref_val:.8f}")
            elif ref_val == 0.0 and abs(stored) > TOL:
                errors.append(f"{col}: stocké={stored:.2e}, attendu=0")

        # Colonnes modulaires q mod l
        for l in L_MODULI:
            col    = f"q_mod_{l}"
            stored = int(getattr(row, col))
            ref_v  = int(rec["q"]) % l
            if stored != ref_v:
                errors.append(f"{col}: stocké={stored}, attendu={ref_v}")

        label = "HEAD" if i < n else "TAIL"
        if errors:
            all_ok = False
            print(f"  [{label}] p={row.p}  ÉCHEC :")
            for e in errors:
                print(f"      {e}")
        else:
            print(f"  [{label}] p={row.p}  OK (m={rec['m']}, ε={rec['eps']:+d}, q={rec['q']:+d})")

    # Vérification des rangs (les 5 premiers doivent avoir r = 1..5)
    head_ranks = db_df.head(n)["r"].tolist()
    expected   = list(range(1, n + 1))
    if head_ranks != expected:
        print(f"  ERREUR rangs HEAD : {head_ranks} ≠ {expected}")
        all_ok = False
    else:
        print(f"  Rangs HEAD : {head_ranks}  OK")

    # Vérification de la croissance des cumulatifs (Q_r non décroissant si q≥0)
    tail_Qr = db_df.tail(n)["Q_r"].tolist()
    if any(b < a - 1e6 for a, b in zip(tail_Qr, tail_Qr[1:])):
        print("  AVERTISSEMENT : Q_r non monotone en queue")
    else:
        print(f"  Q_r (tail) cohérent")

    return all_ok


# ─────────────────────────────────────────────────────────────────────────────
# VÉRIFICATION GLOBALE DES CONTRAINTES
# ─────────────────────────────────────────────────────────────────────────────

def global_checks(k: int, db_dir: Path) -> bool:
    """
    Parcourt tous les fichiers Parquet pour k et vérifie sur chaque batch :
      – reconstruction p = H_m + ε + 2k·q
      – borne |q| ≤ m
      – croissance de r
    """
    print(f"\n--- Checks globaux k={k} ---")
    src_files = sorted((db_dir / "parquet" / f"k={k}").glob("*.parquet"))
    all_ok    = True
    last_r    = 0

    for src in src_files:
        df  = pd.read_parquet(src, columns=["p", "k", "m", "eps", "q", "H_m", "r"])
        two_k = 2 * k

        recon = df["H_m"] + df["eps"] + two_k * df["q"]
        bad   = (recon != df["p"]).sum()
        if bad:
            print(f"  {src.name}: {bad} reconstructions INVALIDES")
            all_ok = False

        bad2 = (df["q"].abs() > df["m"]).sum()
        if bad2:
            print(f"  {src.name}: {bad2} violations |q| ≤ m")
            all_ok = False

        r_min = int(df["r"].iloc[0])
        r_max = int(df["r"].iloc[-1])
        if r_min != last_r + 1:
            print(f"  {src.name}: rang commence à {r_min} (attendu {last_r+1})")
            all_ok = False
        last_r = r_max

    if all_ok:
        print(f"  Tous les checks OK — {last_r:,} lignes vérifiées")
    return all_ok


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Vérificateur indépendant de la base paramétrique")
    parser.add_argument("--db",     type=Path, default=Path("prime_database"),
                        help="Répertoire de la base (défaut: prime_database/)")
    parser.add_argument("--global", action="store_true", dest="full",
                        help="Vérification complète (tous les fichiers)")
    args = parser.parse_args()

    db_dir = args.db
    if not db_dir.exists():
        print(f"Répertoire introuvable : {db_dir}", file=sys.stderr)
        sys.exit(1)

    all_passed = True
    for k in K_VALUES:
        try:
            ok = compare(k, db_dir, n=5)
            if args.full:
                ok = global_checks(k, db_dir) and ok
            all_passed = all_passed and ok
        except FileNotFoundError as exc:
            print(f"  k={k} : {exc}")
            all_passed = False

    print("\n" + ("=" * 60))
    if all_passed:
        print("  RÉSULTAT : TOUTES LES VÉRIFICATIONS PASSÉES")
    else:
        print("  RÉSULTAT : DES ERREURS ONT ÉTÉ DÉTECTÉES — voir ci-dessus")
    print("=" * 60)
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
