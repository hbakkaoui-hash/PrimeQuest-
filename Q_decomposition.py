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
Q_decomposition.py
==================
Décomposition structurelle Q(p_N) ≈ α√p + β·D_N (Section 5 du papier).

Objets centraux :
    Q(p_N) = Σ_{n=1}^{N} c_n,    c_n = q°_n + ε°_n
    D_N    = Σ_{n=1}^{N} ({√(p'_n/3)} - 1/2),    p'_n = p_n - ε°_n

Le Théorème 5.1 (indépendance de convention) établit que le modèle
    Q = α·√p + β·D_N + γ
atteint R² ≥ 0.967 sous les trois conventions viables,
et que la corrélation r(Q_r, D_N) ≥ 0.84 pour toutes les trois.

Cinq conventions testées (Section 5.1) :
    (1) Min |q|   : m* = argmin |q_m|
    (2) Min |c_i| : minimiser |q + ε|
    (3) Canonique : m° = round((-1 + √(1+4p'/3))/2)    ← référence
    (4) Moyenne   : diverge, écarté
    (5) Médiane   : diverge, écarté

Résultats reproduits (Table 6 et Table 7 du papier) :
    Convention    R²(√p+D_N)   r(Q_r,D_N)   Annulation   Hurst H
    Min |q|       0.979        0.846         0.0101 %      0.398
    Min |c_i|     0.980        0.843         0.0111 %      0.395
    Canonique     0.967        0.869         0.0036 %      0.392

Décomposition de variance (Table 7) :
    α√p          →  51.9 % de Var(Q)
    β·D_N        →  44.8 % de Var(Q)
    Résidu pur   →   3.3 % de Var(Q)

Annulation arithmétique (Table 8) :
    Σ p_n ≈ Σ H_{m°}  avec  H_m = 3m(m+1)+1
    Précision croissante avec N : jusqu'à 99.99999997 % à N=10^8.

Usage :
    python Q_decomposition.py               # N=10000 (papier)
    python Q_decomposition.py --quick       # N=2000  (rapide)
    python Q_decomposition.py --N 5000

Dépendances : numpy, scipy, matplotlib
"""

import argparse
import time
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Paramètres
# ---------------------------------------------------------------------------
K     = 3        # famille k=3, p ≡ 1 (mod 6)
N_REF = 10_000   # N du papier


# ---------------------------------------------------------------------------
# Génération des premiers p ≡ 1 (mod 6)
# ---------------------------------------------------------------------------

def primes_mod6(N: int) -> np.ndarray:
    """Retourne les N premiers p > 3 avec p ≡ 1 (mod 6), crible vectorisé."""
    limit = max(int(N * (np.log(N + 10) + np.log(np.log(N + 10))) * 1.3), 200)
    while True:
        sieve = np.ones(limit, dtype=bool)
        sieve[0] = sieve[1] = False
        for i in range(2, int(limit ** 0.5) + 1):
            if sieve[i]:
                sieve[i*i::i] = False
        p_all = np.where(sieve)[0]
        p_sel = p_all[(p_all % 6 == 1) & (p_all > 3)]
        if len(p_sel) >= N:
            return p_sel[:N]
        limit = int(limit * 1.5)


# ---------------------------------------------------------------------------
# Coordonnées paramétriques — trois conventions
# ---------------------------------------------------------------------------

def coords_canonical(p_arr: np.ndarray) -> tuple:
    """
    Convention canonique (3) :
        m° = round((-1 + √(1+4p'/3)) / 2),   p' = p - ε° = p - 1
        q° = (p - 3m°(m°+1) - ε°) / 6
    """
    eps    = np.ones(len(p_arr), dtype=np.int64)      # p ≡ 1 mod 6 → ε° = +1
    p_pr   = p_arr - eps
    m      = np.round((-1 + np.sqrt(1 + 4 * p_pr / 3)) / 2).astype(np.int64)
    q      = ((p_arr - 3 * m * (m + 1) - eps) // 6).astype(np.int64)
    return m, q, eps


def coords_min_absq(p_arr: np.ndarray) -> tuple:
    """
    Convention (1) — Min |q| :
        Pour chaque p, choisir le m qui minimise |q|
        parmi les deux candidats m = floor(...) et ceil(...).
    """
    eps    = np.ones(len(p_arr), dtype=np.int64)
    p_pr   = p_arr - eps
    m_f    = np.floor((-1 + np.sqrt(1 + 4 * p_pr / 3)) / 2).astype(np.int64)
    m_c    = m_f + 1

    q_f    = ((p_arr - 3 * m_f * (m_f + 1) - eps) // 6).astype(np.int64)
    q_c    = ((p_arr - 3 * m_c * (m_c + 1) - eps) // 6).astype(np.int64)

    use_c  = np.abs(q_c) < np.abs(q_f)
    m      = np.where(use_c, m_c, m_f)
    q      = np.where(use_c, q_c, q_f)
    return m, q, eps


def coords_min_absc(p_arr: np.ndarray) -> tuple:
    """
    Convention (2) — Min |c_i| = Min |q + ε| :
        Choisir le m qui minimise |q + ε| = |q + 1|.
    """
    m, q, eps = coords_min_absq(p_arr)
    m_alt      = m + 1
    q_alt      = ((p_arr - 3 * m_alt * (m_alt + 1) - eps) // 6).astype(np.int64)

    ci         = np.abs(q + eps)
    ci_alt     = np.abs(q_alt + eps)
    use_alt    = ci_alt < ci
    m          = np.where(use_alt, m_alt, m)
    q          = np.where(use_alt, q_alt, q)
    return m, q, eps


CONVENTIONS = {
    "Canonique" : coords_canonical,
    "Min |q|"   : coords_min_absq,
    "Min |c_i|" : coords_min_absc,
}


# ---------------------------------------------------------------------------
# Fonctionnelle Q et terme d'équidistribution D_N
# ---------------------------------------------------------------------------

def compute_Q(q: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """Q(p_N) = Σ c_n,   c_n = q°_n + ε°_n."""
    c = q + eps
    return np.cumsum(c.astype(float))


def compute_DN(p_arr: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """
    D_N = Σ_{n=1}^{N} ({√(p'_n/3)} - 1/2),   p'_n = p_n - ε°_n
    {x} = partie fractionnaire de x.
    """
    p_pr  = (p_arr - eps).astype(float)
    frac  = np.modf(np.sqrt(p_pr / 3))[0]    # partie fractionnaire
    return np.cumsum(frac - 0.5)


# ---------------------------------------------------------------------------
# Régression OLS : Q ~ α√p + β·D_N + γ
# ---------------------------------------------------------------------------

def ols_decomposition(Q: np.ndarray, p_arr: np.ndarray,
                      D_N: np.ndarray) -> dict:
    """
    Décomposition en deux étapes (conforme au papier) :

    Étape 1 : Régression Q ~ α·√p + γ₀  → résidu Q_r = Q - α√p - γ₀
              → R²(√p), correlation r(Q_r, D_N)
    Étape 2 : Régression Q ~ α·√p + β·D_N + γ  (modèle complet)
              → R²(√p + D_N), résidu pur Q_rr = Q - ajustement complet
              → Décomposition de variance.

    r(Q_r, D_N) calculé sur le résidu de l'étape 1 (pas de l'étape 2),
    sinon il est nul par construction (D_N est régresseur de l'étape 2).
    """
    sqrt_p = np.sqrt(p_arr)

    # --- Étape 1 : tendance seule ---
    X1     = np.column_stack([sqrt_p, np.ones(len(p_arr))])
    b1, _, _, _ = np.linalg.lstsq(X1, Q, rcond=None)
    alpha1, gamma1 = b1
    Q_r    = Q - X1 @ b1              # résidu après suppression de la tendance

    # R²(√p)
    r2_trend = 1 - np.var(Q_r) / np.var(Q)

    # Corrélation résidu (étape 1) avec D_N
    r_QrDN = float(np.corrcoef(Q_r, D_N)[0, 1])

    # --- Étape 2 : modèle complet ---
    X2     = np.column_stack([sqrt_p, D_N, np.ones(len(p_arr))])
    b2, _, _, _ = np.linalg.lstsq(X2, Q, rcond=None)
    alpha, beta_d, gamma = b2
    Q_rr   = Q - X2 @ b2             # résidu pur

    # R²(√p + D_N)
    r2_full = 1 - np.var(Q_rr) / np.var(Q)

    # Décomposition de variance (modèle complet)
    trend_var = np.var(alpha * sqrt_p)
    equid_var = np.var(beta_d * D_N)
    resid_var = np.var(Q_rr)
    total_var = np.var(Q)

    var_decomp = {
        "trend_pct" : 100 * trend_var / total_var,
        "equid_pct" : 100 * equid_var / total_var,
        "resid_pct" : 100 * resid_var / total_var,
    }

    return {
        "alpha"    : alpha,
        "beta"     : beta_d,
        "gamma"    : gamma,
        "r2"       : r2_full,          # R² du modèle complet
        "r2_trend" : r2_trend,         # R² tendance seule
        "r_QrDN"   : r_QrDN,          # corrélation résidu-tendance vs D_N
        "Q_r"      : Q_rr,             # résidu pur (modèle complet)
        "var_decomp": var_decomp,
    }


# ---------------------------------------------------------------------------
# Exposant de Hurst (statistique R/S)
# ---------------------------------------------------------------------------

def hurst_rs(x: np.ndarray, block_size: int = None) -> float:
    """
    Estime l'exposant de Hurst H via la statistique R/S.
    block_size ≈ √N selon [Beran 1994].
    H < 0.5 : anti-persistance (attendu H ≈ 0.39).
    """
    N = len(x)
    if block_size is None:
        block_size = max(10, int(np.sqrt(N)))

    n_blocks = N // block_size
    rs_vals  = []

    for b in range(n_blocks):
        seg   = x[b * block_size : (b + 1) * block_size].astype(float)
        seg   = seg - seg.mean()                 # centrer
        cs    = np.cumsum(seg)                   # somme cumulée
        R     = cs.max() - cs.min()              # étendue
        S     = seg.std(ddof=1)                  # écart-type
        if S > 0:
            rs_vals.append(R / S)

    if len(rs_vals) < 2:
        return float("nan")

    # H = ln(R/S) / ln(block_size)
    H = np.log(np.mean(rs_vals)) / np.log(block_size)
    return float(H)


# ---------------------------------------------------------------------------
# Annulation arithmétique Σp_n ≈ Σ H_{m°}
# ---------------------------------------------------------------------------

def arithmetic_cancellation(p_arr: np.ndarray, m_arr: np.ndarray) -> dict:
    """
    Calcule l'annulation :   |Σ p_n - Σ H_{m°}| / Σ p_n
    avec H_m = 3m(m+1)+1.

    Reproduit la Table 8 du papier.
    """
    H_m   = 3 * m_arr * (m_arr + 1) + 1
    A     = float(p_arr.sum())
    B     = float(H_m.sum())
    diff  = abs(A - B)
    rel   = diff / A
    return {
        "sum_p"       : A,
        "sum_H"       : B,
        "abs_diff"    : diff,
        "rel_diff"    : rel,
        "cancellation": 1 - rel,
    }


# ---------------------------------------------------------------------------
# Affichage des résultats
# ---------------------------------------------------------------------------

def print_results(results_all: dict, N: int):
    sep = "-" * 64

    print("\n" + "=" * 64)
    print(f"DÉCOMPOSITION Q(p_N) ≈ α√p + β·D_N   (N = {N:,} premiers, k=3)")
    print("=" * 64)

    # Table 6 : R², r(Q_r, D_N), annulation, Hurst
    print(f"\n{'Convention':<14} {'R²':>7} {'r(Qr,DN)':>10} "
          f"{'Annulation':>12} {'Hurst H':>9}")
    print(sep)
    for name, res in results_all.items():
        r2      = res["r2"]
        r_qrdn  = res["r_QrDN"]
        cancel  = res["cancellation"]["cancellation"] * 100
        H       = res["hurst"]
        print(f"  {name:<14} {r2:.4f}   {r_qrdn:>+.4f}   {cancel:.4f} %   {H:.4f}")

    print(f"\n  Valeurs du papier (Table 6) :")
    ref = {
        "Canonique" : (0.967, 0.869, 99.9964, 0.392),
        "Min |q|"   : (0.979, 0.846, 99.9899, 0.398),
        "Min |c_i|" : (0.980, 0.843, 99.9889, 0.395),
    }
    for name, (r2, r, c, h) in ref.items():
        print(f"    {name:<14} R²={r2:.3f}  r={r:+.3f}  Ann.={c:.4f}%  H={h:.3f}")

    # Table 7 : décomposition de variance (convention canonique)
    print(f"\n{'='*64}")
    print("DÉCOMPOSITION DE VARIANCE (convention canonique, Table 7)")
    print("=" * 64)
    vd = results_all["Canonique"]["var_decomp"]
    print(f"  α√p          → {vd['trend_pct']:5.1f} % de Var(Q)  "
          f"[papier : 51.9 %]")
    print(f"  β·D_N        → {vd['equid_pct']:5.1f} % de Var(Q)  "
          f"[papier : 44.8 %]")
    print(f"  Résidu pur   → {vd['resid_pct']:5.1f} % de Var(Q)  "
          f"[papier :  3.3 %]")

    # Annulation (convention canonique)
    ca = results_all["Canonique"]["cancellation"]
    print(f"\n{'='*64}")
    print("ANNULATION ARITHMÉTIQUE Σp_n ≈ Σ H_m°  (Table 8)")
    print("=" * 64)
    print(f"  N = {N:,} premiers")
    print(f"  Σ p_n          = {ca['sum_p']:.2f}")
    print(f"  Σ H_{{m°}}      = {ca['sum_H']:.2f}")
    print(f"  |A-B|          = {ca['abs_diff']:.2f}")
    print(f"  |A-B|/A        = {ca['rel_diff']:.2e}")
    print(f"  Annulation     = {ca['cancellation']*100:.4f} %")
    print(f"  [papier N=10⁴] : 99.9964 %  (|A-B|/A = 3.6×10⁻⁵)")

    # Coefficients
    print(f"\n{'='*64}")
    print("COEFFICIENTS OLS (convention canonique)")
    print("=" * 64)
    res = results_all["Canonique"]
    print(f"  α (√p)  = {res['alpha']:.6f}")
    print(f"  β (D_N) = {res['beta']:.6f}")
    print(f"  γ (cte) = {res['gamma']:.4f}")
    print(f"  R²       = {res['r2']:.6f}  [papier : 0.967]")
    print(f"  r(Qr,DN) = {res['r_QrDN']:+.6f}  [papier : 0.869]")


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def plot_results(results_all: dict, p_arr: np.ndarray, N: int,
                 out: str = "fig_Q_decomp_reprod.png"):
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.rcParams.update({
            "font.size": 10, "axes.labelsize": 11,
            "figure.dpi": 150, "savefig.dpi": 300,
        })
    except ImportError:
        print("matplotlib non disponible — visualisation ignorée.")
        return

    res    = results_all["Canonique"]
    Q      = res["Q"]
    D_N    = res["D_N"]
    Q_r    = res["Q_r"]
    sqrt_p = np.sqrt(p_arr)
    alpha  = res["alpha"]
    beta_d = res["beta"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # (A) Q vs α√p
    ax = axes[0, 0]
    trend = alpha * sqrt_p
    ax.plot(sqrt_p, Q, lw=0.5, color="#2166ac", alpha=0.7, label="Q(p_N)")
    ax.plot(sqrt_p, trend, lw=1.2, color="red", label=rf"α√p  (R²={1-np.var(Q-trend)/np.var(Q):.4f})")
    ax.set_xlabel(r"$\sqrt{p}$")
    ax.set_ylabel("Q")
    ax.set_title(r"Tendance dominante $Q \approx \alpha\sqrt{p}$")
    ax.legend(fontsize=8)

    # (B) Résidu Q_r vs β·D_N
    ax = axes[0, 1]
    ax.scatter(beta_d * D_N, Q_r, s=0.3, alpha=0.3, color="#4daf4a")
    lim = max(abs(Q_r).max(), abs(beta_d * D_N).max()) * 1.05
    ax.plot([-lim, lim], [-lim, lim], "r--", lw=1.0, label=f"r={res['r_QrDN']:+.4f}")
    ax.set_xlabel(r"$\beta \cdot D_N$")
    ax.set_ylabel(r"$Q_r = Q - \alpha\sqrt{p}$")
    ax.set_title(r"Résidu vs équidistribution $D_N$")
    ax.legend(fontsize=8)

    # (C) Décomposition de variance — barres
    ax = axes[1, 0]
    vd     = res["var_decomp"]
    labels = ["α√p\n(tendance)", "β·D_N\n(équidist.)", "Résidu pur"]
    vals   = [vd["trend_pct"], vd["equid_pct"], vd["resid_pct"]]
    colors = ["#1f78b4", "#33a02c", "#e31a1c"]
    bars   = ax.bar(labels, vals, color=colors, width=0.5)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f"{val:.1f} %", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("% de Var(Q)")
    ax.set_title(f"Décomposition de variance (N={N:,})")
    ax.set_ylim(0, 65)

    # (D) D_N : trajectoire et enveloppe de marche aléatoire
    ax = axes[1, 1]
    n_idx  = np.arange(1, N + 1)
    rw_env = 2 * np.sqrt(n_idx)
    ax.plot(n_idx, D_N, lw=0.6, color="#984ea3", alpha=0.8, label=r"$D_N$")
    ax.fill_between(n_idx, -rw_env, rw_env, alpha=0.15, color="#4daf4a",
                    label=r"Enveloppe $\pm 2\sqrt{N}$")
    ax.set_xlabel("N")
    ax.set_ylabel(r"$D_N$")
    ax.set_title(r"Équidistribution de Weyl : $\{\sqrt{p/3}\}$ sur $[0,1]$")
    ax.legend(fontsize=8)

    plt.suptitle(f"Décomposition $Q(p_N) \\approx \\alpha\\sqrt{{p}} + \\beta D_N$"
                 f" — k=3, N={N:,} premiers", fontsize=12, y=1.01)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    print(f"  Figure sauvegardée : {out}")
    plt.close()


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=r"Décomposition Q(p_N) ≈ α√p + β·D_N")
    parser.add_argument("--N", type=int, default=None,
                        help=f"Nombre de premiers (défaut : {N_REF})")
    parser.add_argument("--quick", action="store_true",
                        help="Mode rapide : N=2000")
    parser.add_argument("--no-plot", action="store_true",
                        help="Désactiver les figures")
    args = parser.parse_args()

    N = args.N if args.N else (2_000 if args.quick else N_REF)

    print("=" * 64)
    print(f"Décomposition Q(p_N) ≈ α√p + β·D_N")
    print(f"k={K},  N={N:,} premiers p ≡ 1 (mod 6)")
    print("=" * 64)

    # Génération des premiers
    print(f"\nGénération de {N:,} premiers ...")
    t0    = time.time()
    p_arr = primes_mod6(N)
    print(f"  Fait en {time.time()-t0:.2f}s  "
          f"(p_1={p_arr[0]}, p_N={p_arr[-1]:,})")

    results_all = {}

    for name, coord_fn in CONVENTIONS.items():
        print(f"\nConvention : {name}")
        m, q, eps = coord_fn(p_arr)

        Q   = compute_Q(q, eps)
        D_N = compute_DN(p_arr, eps)

        # OLS
        res = ols_decomposition(Q, p_arr, D_N)
        print(f"  R² = {res['r2']:.4f}   r(Q_r, D_N) = {res['r_QrDN']:+.4f}")

        # Hurst
        c = (q + eps).astype(float)
        H = hurst_rs(c, block_size=int(np.sqrt(N)))
        print(f"  Hurst H = {H:.4f}  (papier : ~0.39, anti-persistance)")

        # Annulation
        ca = arithmetic_cancellation(p_arr, m)
        print(f"  Annulation : {ca['cancellation']*100:.4f} %  "
              f"(|A-B|/A = {ca['rel_diff']:.2e})")

        results_all[name] = {**res, "hurst": H, "cancellation": ca,
                             "Q": Q, "D_N": D_N}

    # Rapport complet
    print_results(results_all, N)

    # Figures
    if not args.no_plot:
        print("\nGénération des figures ...")
        plot_results(results_all, p_arr, N)

    print("\nFin de la décomposition.")


if __name__ == "__main__":
    main()
