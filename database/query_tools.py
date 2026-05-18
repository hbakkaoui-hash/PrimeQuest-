#!/usr/bin/env python3
# =============================================================================
# query_tools.py — Fonctions de requête pour la base de données paramétrique
#
# Formule : p = k·m(m+1) + ε + 2k·q
# =============================================================================
"""
Fonctions utilitaires pour interroger la base Parquet générée par
generate_database.py.

Fonctions exportées :
  query_layer(m, k)              → DataFrame de tous les premiers de la couche m
  compute_c(m_min, m_max, k)     → statistique χ² cumulée sur [m_min, m_max]

Toutes les fonctions lisent directement les fichiers Parquet partitionnés
et ne nécessitent pas que SQLite soit ouvert.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.dataset as ds

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

DEFAULT_DB_DIR = Path("prime_database")
K_VALUES       = [3, 15, 21, 23]
PHI_2K         = {3: 2, 15: 8, 21: 12, 23: 22}


def _parquet_dir(k: int, db_dir: Path = DEFAULT_DB_DIR) -> Path:
    return db_dir / "parquet" / f"k={k}"


def _files_for_k(k: int, db_dir: Path = DEFAULT_DB_DIR) -> list[Path]:
    d = _parquet_dir(k, db_dir)
    if not d.exists():
        raise FileNotFoundError(
            f"Aucun répertoire Parquet pour k={k}: {d}\n"
            "Avez-vous lancé generate_database.py ?"
        )
    return sorted(d.glob("*.parquet"))


# ─────────────────────────────────────────────────────────────────────────────
# query_layer
# ─────────────────────────────────────────────────────────────────────────────

def query_layer(
    m: int,
    k: int,
    db_dir: Path = DEFAULT_DB_DIR,
    columns: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    Retourne un DataFrame contenant tous les premiers de la couche m
    pour le paramètre k, dans l'ordre croissant de p.

    Paramètres
    ----------
    m        : indice de couche (entier ≥ 0)
    k        : paramètre de famille ∈ {3, 15, 21, 23}
    db_dir   : répertoire racine de la base (défaut : prime_database/)
    columns  : liste de colonnes à retourner (défaut : toutes)

    Retourne
    --------
    pd.DataFrame  (vide si la couche est vide dans la base générée)

    Exemple
    -------
    >>> df = query_layer(m=100, k=3)
    >>> print(df[["p", "eps", "q", "Q_r"]].head())
    """
    if k not in K_VALUES:
        raise ValueError(f"k={k} invalide ; valeurs autorisées : {K_VALUES}")

    files = _files_for_k(k, db_dir)

    # Filtre PyArrow poussé vers le bas (pushdown predicate)
    filt  = ds.field("m") == pa.scalar(m, type=pa.int64())
    parts: list[pd.DataFrame] = []

    for f in files:
        # Read with pandas and filter; avoids schema-merge issues across files
        chunk = pd.read_parquet(f, columns=columns)
        chunk = chunk[chunk["m"] == m]
        if not chunk.empty:
            parts.append(chunk)

    if not parts:
        return pd.DataFrame()

    result = pd.concat(parts, ignore_index=True).sort_values("p")
    return result.reset_index(drop=True)


# ─────────────────────────────────────────────────────────────────────────────
# compute_c
# ─────────────────────────────────────────────────────────────────────────────

def compute_c(
    m_min: int,
    m_max: int,
    k: int,
    db_dir: Path = DEFAULT_DB_DIR,
) -> Tuple[float, pd.DataFrame]:
    """
    Calcule la statistique de chi-deux cumulée sur la plage de couches
    [m_min, m_max] pour le paramètre k.

    c(m_min, m_max) = (1 / (m_max − m_min + 1)) · Σ_{m=m_min}^{m_max} T_m²

    où T_m = (n_m − n̂_m) / √n̂_m est la statistique de Poisson normalisée.

    Sous l'hypothèse nulle d'indépendance, T_m ~ N(0,1) et c ≈ 1.
    Des valeurs c >> 1 indiquent une surdispersion.

    Paramètres
    ----------
    m_min, m_max : bornes inclusives de la plage de couches
    k            : paramètre de famille ∈ {3, 15, 21, 23}
    db_dir       : répertoire racine de la base

    Retourne
    --------
    (c, detail_df) où :
      c         : float — statistique chi-deux normalisée
      detail_df : DataFrame avec colonnes [m, n_m, nm_hat, T_m, T_m_sq]
                  une ligne par couche dans [m_min, m_max]

    Exemple
    -------
    >>> c_val, detail = compute_c(m_min=10, m_max=1000, k=3)
    >>> print(f"c = {c_val:.4f}")
    """
    if k not in K_VALUES:
        raise ValueError(f"k={k} invalide ; valeurs autorisées : {K_VALUES}")
    if m_min > m_max:
        raise ValueError(f"m_min={m_min} > m_max={m_max}")

    files = _files_for_k(k, db_dir)

    rows: list[dict] = []
    seen_layers: set[int] = set()

    for f in files:
        chunk = pd.read_parquet(f, columns=["m", "n_m", "nm_hat", "T_m"])
        df    = chunk[(chunk["m"] >= m_min) & (chunk["m"] <= m_max)]
        if df.empty:
            continue

        for m_val, grp in df.groupby("m"):
            if m_val in seen_layers:
                continue
            seen_layers.add(m_val)
            rows.append({
                "m":      int(m_val),
                "n_m":    int(grp["n_m"].iloc[0]),
                "nm_hat": float(grp["nm_hat"].iloc[0]),
                "T_m":    float(grp["T_m"].iloc[0]),
            })

    if not rows:
        return 0.0, pd.DataFrame(columns=["m", "n_m", "nm_hat", "T_m", "T_m_sq"])

    detail = pd.DataFrame(rows).sort_values("m").reset_index(drop=True)
    detail["T_m_sq"] = detail["T_m"] ** 2

    n_layers = len(detail)
    c_val    = float(detail["T_m_sq"].sum()) / n_layers

    return c_val, detail


# ─────────────────────────────────────────────────────────────────────────────
# UTILITAIRES ADDITIONNELS
# ─────────────────────────────────────────────────────────────────────────────

def layer_stats(m: int, k: int, db_dir: Path = DEFAULT_DB_DIR) -> dict:
    """
    Retourne un dictionnaire de statistiques pour une couche (k, m) :
      n_m, nm_hat, T_m, sum_q, mean_q, std_q, min_q, max_q, Q_r_last
    """
    df = query_layer(m, k, db_dir)
    if df.empty:
        return {"m": m, "k": k, "n_m": 0}

    return {
        "m":       m,
        "k":       k,
        "H_m":     int(df["H_m"].iloc[0]),
        "n_m":     int(df["n_m"].iloc[0]),
        "nm_hat":  float(df["nm_hat"].iloc[0]),
        "T_m":     float(df["T_m"].iloc[0]),
        "sum_q":   float(df["q"].sum()),
        "mean_q":  float(df["q"].mean()),
        "std_q":   float(df["q"].std(ddof=1)) if len(df) > 1 else 0.0,
        "min_q":   int(df["q"].min()),
        "max_q":   int(df["q"].max()),
        "Q_r_last": float(df["Q_r"].iloc[-1]),
    }


def prime_in_family(p: int, k: int, db_dir: Path = DEFAULT_DB_DIR) -> Optional[pd.Series]:
    """
    Retourne la ligne de la base correspondant au premier p dans la k-famille,
    ou None si p n'est pas dans la base.
    """
    files = _files_for_k(k, db_dir)
    for f in files:
        tbl = pq.read_table(f, filters=[("p", "=", p)])
        if tbl.num_rows > 0:
            return tbl.to_pandas().iloc[0]
    return None


def c_scan(
    k: int,
    m_step: int = 1000,
    db_dir: Path = DEFAULT_DB_DIR,
) -> pd.DataFrame:
    """
    Balayage de c(m, m+m_step) sur toutes les couches disponibles.
    Utile pour visualiser l'évolution de la dispersion par tranches.

    Retourne un DataFrame [m_min, m_max, c, n_layers].
    """
    files = _files_for_k(k, db_dir)
    if not files:
        return pd.DataFrame()

    # Trouver m_min et m_max dans la base
    first_row = pd.read_parquet(files[0],  columns=["m"]).head(1)
    last_row  = pd.read_parquet(files[-1], columns=["m"]).tail(1)
    m_global_min = int(first_row["m"].iloc[0])
    m_global_max = int(last_row["m"].iloc[0])

    records = []
    m = m_global_min
    while m <= m_global_max:
        m_end = min(m + m_step - 1, m_global_max)
        c_val, det = compute_c(m, m_end, k, db_dir)
        records.append({
            "m_min":    m,
            "m_max":    m_end,
            "c":        c_val,
            "n_layers": len(det),
        })
        m += m_step

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────────────────────
# DÉMONSTRATION CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Requêtes sur la base paramétrique")
    parser.add_argument("--k",     type=int, default=3,  help="Paramètre k (défaut 3)")
    parser.add_argument("--m",     type=int, default=10, help="Couche m pour query_layer")
    parser.add_argument("--mmin",  type=int, default=1,  help="m_min pour compute_c")
    parser.add_argument("--mmax",  type=int, default=100, help="m_max pour compute_c")
    parser.add_argument("--db",    type=Path, default=DEFAULT_DB_DIR)
    args = parser.parse_args()

    print(f"\n--- query_layer(m={args.m}, k={args.k}) ---")
    df = query_layer(args.m, args.k, db_dir=args.db)
    if df.empty:
        print("  (couche vide)")
    else:
        print(df[["p", "m", "eps", "q", "r", "Q_r", "X_r", "D_r"]].to_string(index=False))

    print(f"\n--- compute_c(m_min={args.mmin}, m_max={args.mmax}, k={args.k}) ---")
    c_val, detail = compute_c(args.mmin, args.mmax, args.k, db_dir=args.db)
    print(f"  c = {c_val:.6f}  ({len(detail)} couches)")
    if not detail.empty:
        print(detail.head(10).to_string(index=False))
