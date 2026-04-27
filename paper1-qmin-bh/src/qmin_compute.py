# =============================================================================
# Auteur  : Hassane BAKKAOUI — Chercheur indépendant
# Contact : bakkahassa@hotmail.com
# Dépôt   : https://github.com/BAKKAOUI/parametric-primes-2026
# Papier  : "Une famille paramétrique de nombres premiers via les formes
#            quadratiques — Distribution heuristique du décalage minimal,
#            validation numérique à grande échelle et contrôle statistique
#            des corrélations spectrales spurieuses"
# Auteur  : H. Bakkaoui (2026)
# arXiv   : math.NT [à compléter après soumission]
# Licence : MIT License — voir LICENSE dans le dépôt
# Date    : Avril 2026
# Version : 1.0.0
# =============================================================================
"""
qmin_compute.py
===============
Calcul du décalage minimal q_min(21, m) sur m ∈ [0, 10^6].

Famille paramétrique :
    p_{k,m,ε,q} = k·m·(m+1) + ε + 2k·q,   k=21, ε ∈ {±1}, q ∈ Z

Décalage minimal :
    q_min(k,m) = argmin_{q ∈ Z} { |q| : p_{k,m,ε,q} premier }

Stratégie de recherche :
    - Bidirectionnelle : on teste q = 0, ±1, ±2, ... jusqu'à trouver un premier.
    - À q fixé, on teste les quatre candidats : (ε=+1, q), (ε=-1, q),
      (ε=+1, -q), (ε=-1, -q).
    - En cas d'égalité |q|, on choisit ε=+1 par convention.
    - Borne garantie : |q_min| ≤ m (Proposition 2.3 du papier).

Test de primalité :
    BPSW (SymPy isprime), déterministe pour n < 3.3 × 10^24.

Résultats reproduits dans le papier (Section 5) :
    - Plage de q      : [-46, +50]
    - Moyenne |q|     : 3.711
    - Médiane |q|     : 3.0
    - Écart-type |q|  : 3.623
    - Taux q=0        : 13.59 %
    - Pente log. (OLS): 0.320  [IC 95 % : 0.31, 0.33],  R² = 0.984
    - ρ_21            : 0.320 × C_21 = 0.320 × 8.206 = 2.626 ± 0.08

Usage :
    python qmin_compute.py              # calcul complet M=10^6 (long)
    python qmin_compute.py --quick      # démonstration M=10^4

Dépendances : sympy, numpy, scipy, matplotlib
"""

import argparse
import time
import numpy as np
from scipy import stats
from sympy import isprime

# ---------------------------------------------------------------------------
# Paramètres globaux
# ---------------------------------------------------------------------------
K      = 21                # valeur de k étudiée
C_K    = 8.206301          # constante de Bateman–Horn agrégée (Table 2)
Q_BOUND = 1000             # borne de recherche (jamais atteinte en pratique)


# ---------------------------------------------------------------------------
# Calcul de q_min pour un m donné
# ---------------------------------------------------------------------------

def qmin_single(m: int, k: int = K, q_bound: int = Q_BOUND) -> int:
    """
    Retourne q_min(k, m) : le plus petit |q| tel que
    p_{k,m,ε,q} = k·m·(m+1) + ε + 2k·q soit premier.

    En cas d'égalité de |q|, ε=+1 est prioritaire.
    Retourne q signé (positif si ε=+1 gagne, négatif sinon).
    """
    base = k * m * (m + 1)
    step = 2 * k

    # q = 0 : tester ε=+1 puis ε=-1
    if isprime(base + 1):
        return 0
    if isprime(base - 1):
        return 0

    for q in range(1, q_bound + 1):
        # tester +q avec ε=+1 (priorité convention)
        if isprime(base + 1 + step * q):
            return q
        # tester -q avec ε=+1
        if isprime(base + 1 - step * q):
            return -q
        # tester +q avec ε=-1
        if isprime(base - 1 + step * q):
            return q
        # tester -q avec ε=-1
        if isprime(base - 1 - step * q):
            return -q

    raise RuntimeError(f"q_bound={q_bound} atteint pour m={m} — augmenter Q_BOUND")


# ---------------------------------------------------------------------------
# Calcul vectorisé sur [0, M]
# ---------------------------------------------------------------------------

def compute_qmin_array(M: int, k: int = K,
                       verbose: bool = True) -> np.ndarray:
    """
    Calcule q_min(k, m) pour m = 0, 1, ..., M.
    Retourne un tableau numpy de longueur M+1.
    """
    qmin = np.zeros(M + 1, dtype=np.int32)
    t0   = time.time()

    for m in range(M + 1):
        qmin[m] = qmin_single(m, k)

        if verbose and m > 0 and m % max(1, M // 20) == 0:
            elapsed = time.time() - t0
            rate    = m / elapsed
            eta     = (M - m) / rate if rate > 0 else 0
            print(f"  m = {m:>8,} / {M:,}   |q|_moy = {np.mean(np.abs(qmin[:m+1])):.4f}"
                  f"   vitesse = {rate:.0f} m/s   ETA = {eta:.0f}s")

    return qmin


# ---------------------------------------------------------------------------
# Statistiques descriptives (Table 3 du papier)
# ---------------------------------------------------------------------------

def descriptive_stats(qmin: np.ndarray) -> dict:
    """Reproduit les statistiques de la Table 3 du papier."""
    aq = np.abs(qmin)
    N  = len(qmin)
    s  = {
        "N"             : N,
        "min_q"         : int(qmin.min()),
        "max_q"         : int(qmin.max()),
        "mean_abs_q"    : float(aq.mean()),
        "median_abs_q"  : float(np.median(aq)),
        "std_abs_q"     : float(aq.std(ddof=1)),
        "n_zero"        : int((qmin == 0).sum()),
        "n_pos"         : int((qmin > 0).sum()),
        "n_neg"         : int((qmin < 0).sum()),
        "pct_zero"      : 100 * (qmin == 0).sum() / N,
        "pct_pos"       : 100 * (qmin > 0).sum() / N,
        "pct_neg"       : 100 * (qmin < 0).sum() / N,
        "imbalance_pct" : 100 * abs((qmin > 0).sum() - (qmin < 0).sum()) / N,
    }
    return s


def print_stats(s: dict):
    print("\n--- Statistiques descriptives (Table 3) ---")
    print(f"  N                 : {s['N']:,}")
    print(f"  Plage de q        : [{s['min_q']}, +{s['max_q']}]")
    print(f"  Moyenne |q|       : {s['mean_abs_q']:.3f}")
    print(f"  Médiane |q|       : {s['median_abs_q']:.1f}")
    print(f"  Écart-type |q|    : {s['std_abs_q']:.3f}")
    print(f"  q = 0             : {s['n_zero']:,}  ({s['pct_zero']:.2f} %)")
    print(f"  q > 0             : {s['n_pos']:,}  ({s['pct_pos']:.2f} %)")
    print(f"  q < 0             : {s['n_neg']:,}  ({s['pct_neg']:.2f} %)")
    print(f"  Déséquilibre ±    : {s['imbalance_pct']:.2f} %")


# ---------------------------------------------------------------------------
# Loi logarithmique empirique (Observation 5.1 du papier)
# ---------------------------------------------------------------------------

def fit_log_law(qmin: np.ndarray, M: int,
                n_bins: int = 50,
                m_min: int = 1000) -> dict:
    """
    Ajuste E[|q_min|](m) ≈ a·ln(m) + b par OLS sur 50 tranches log.

    Reproduit l'Observation 5.1 : pente = 0.320, R² = 0.984.
    """
    # Grille logarithmique
    edges = np.logspace(np.log10(m_min), np.log10(M), n_bins + 1).astype(int)
    edges = np.unique(edges)

    m_centers = []
    mean_aq   = []

    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        idx = np.arange(lo, min(hi + 1, M + 1))
        if len(idx) < 10:
            continue
        m_centers.append(np.sqrt(lo * hi))           # centre géométrique
        mean_aq.append(np.mean(np.abs(qmin[idx])))

    m_centers = np.array(m_centers)
    mean_aq   = np.array(mean_aq)
    ln_m      = np.log(m_centers)

    slope, intercept, r, p_val, se = stats.linregress(ln_m, mean_aq)
    r2    = r ** 2
    rho21 = slope * C_K

    return {
        "slope"    : slope,
        "intercept": intercept,
        "r2"       : r2,
        "se_slope" : se,
        "rho_21"   : rho21,
        "m_centers": m_centers,
        "mean_aq"  : mean_aq,
    }


def print_log_law(fit: dict):
    slope = fit["slope"]
    b     = fit["intercept"]
    r2    = fit["r2"]
    se    = fit["se_slope"]
    rho   = fit["rho_21"]
    t975  = 1.96                        # approximation IC 95 %
    print("\n--- Loi logarithmique empirique (Observation 5.1) ---")
    print(f"  E[|q_min|](m) ≈ {slope:.4f}·ln(m) + ({b:.4f})")
    print(f"  R²             = {r2:.4f}")
    print(f"  IC 95 % pente  ≈ [{slope - t975*se:.3f}, {slope + t975*se:.3f}]")
    print(f"  ρ_21 = pente × C_21 = {slope:.4f} × {C_K} = {rho:.4f}")
    print(f"  Valeur du papier : pente = 0.320, R² = 0.984, ρ_21 = 2.626")


# ---------------------------------------------------------------------------
# Convergence de la moyenne cumulative
# ---------------------------------------------------------------------------

def cumulative_mean(qmin: np.ndarray) -> np.ndarray:
    """Moyenne cumulative de |q_min| — converge vers la vraie espérance."""
    aq = np.abs(qmin).astype(float)
    return np.cumsum(aq) / (np.arange(1, len(aq) + 1))


# ---------------------------------------------------------------------------
# Visualisation (optionnelle)
# ---------------------------------------------------------------------------

def plot_results(qmin: np.ndarray, fit: dict, M: int):
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.rcParams.update({
            "font.size": 11, "axes.labelsize": 12,
            "figure.dpi": 150, "savefig.dpi": 300,
        })
    except ImportError:
        print("matplotlib non disponible — visualisation ignorée.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # --- Panneau gauche : convergence de la moyenne cumulative ---
    ax = axes[0]
    cm = cumulative_mean(qmin)
    ax.plot(cm, lw=0.6, color="#2166ac", alpha=0.8)
    ax.axhline(cm[-1], color="red", ls="--", lw=1.2,
               label=f"Valeur finale = {cm[-1]:.3f}")
    ax.set_xlabel("m")
    ax.set_ylabel(r"$\bar{|q|}$ cumulatif")
    ax.set_title(r"Convergence de $\mathbb{E}[|q_{\min}|]$, $k=21$")
    ax.legend(fontsize=9)
    ax.set_xscale("log")

    # --- Panneau droit : loi logarithmique ---
    ax = axes[1]
    m_c = fit["m_centers"]
    maq = fit["mean_aq"]
    ln_m_grid = np.linspace(np.log(m_c[0]), np.log(m_c[-1]), 200)
    fit_line  = fit["slope"] * ln_m_grid + fit["intercept"]

    ax.scatter(np.log(m_c), maq, s=20, color="#4daf4a", label="Données (50 tranches)")
    ax.plot(ln_m_grid, fit_line, color="red", lw=1.5,
            label=rf"Ajustement : {fit['slope']:.3f}·ln(m){fit['intercept']:+.3f}"
                  f"\n$R^2$={fit['r2']:.4f}")
    ax.set_xlabel(r"$\ln m$")
    ax.set_ylabel(r"$\mathbb{E}[|q_{\min}|]$")
    ax.set_title(rf"Loi logarithmique — $k={K}$, $M={M:,}$")
    ax.legend(fontsize=9)

    plt.tight_layout()
    out = f"qmin_k{K}_M{M}.png"
    plt.savefig(out)
    print(f"\n  Figure sauvegardée : {out}")
    plt.close()


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Calcul de q_min(21,m) et validation de la loi logarithmique")
    parser.add_argument("--quick", action="store_true",
                        help="Mode rapide : M=10^4 pour tests")
    parser.add_argument("--M", type=int, default=None,
                        help="Valeur maximale de m (défaut : 10^6 ou 10^4 si --quick)")
    parser.add_argument("--no-plot", action="store_true",
                        help="Désactiver la génération de figures")
    args = parser.parse_args()

    M = args.M if args.M else (10_000 if args.quick else 1_000_000)

    print("=" * 60)
    print(f"q_min(k={K}, m),  m ∈ [0, {M:,}]")
    print(f"Test de primalité : BPSW (SymPy isprime)")
    print("=" * 60)

    # Calcul
    t0   = time.time()
    qmin = compute_qmin_array(M, K, verbose=True)
    dur  = time.time() - t0
    print(f"\nTemps total : {dur:.1f}s  ({(M+1)/dur:.0f} valeurs/s)")

    # Statistiques descriptives
    s = descriptive_stats(qmin)
    print_stats(s)

    # Loi logarithmique
    m_min = max(100, M // 1000)
    fit   = fit_log_law(qmin, M, n_bins=50, m_min=m_min)
    print_log_law(fit)

    # Sauvegarde des données brutes
    out_npy = f"qmin_k{K}_M{M}.npy"
    np.save(out_npy, qmin)
    print(f"\n  Données sauvegardées : {out_npy}")

    # Figures
    if not args.no_plot:
        plot_results(qmin, fit, M)

    print("\nFin du calcul.")


if __name__ == "__main__":
    main()
