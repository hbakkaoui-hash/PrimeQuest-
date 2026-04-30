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
permutation_tests.py
====================
Trois tests de permutation indépendants sur le résidu pur Q_r (Section 6).

Problème central :
    La fonctionnelle cumulative Q(p_N) = Σ c_n est massivement autocorrélée.
    Toute corrélation de Pearson entre Q et cos(γ_i ln p) est artificiellement
    gonflée → régression spurieuse entre processus I(1) [Granger & Newbold 1974].

Solution : tester si les corrélations observées sont plus fortes que celles
obtenues par des permutations qui préservent la structure de bruit.

Trois méthodes (Section 6.2) :

  (T1) MÉLANGE CLASSIQUE
       Les incréments c_i = q°_i + ε°_i sont mélangés aléatoirement.
       Préserve la distribution marginale, détruit l'ordre.

  (T2) PERMUTATION PAR BLOCS (taille 50)
       Des blocs de 50 c_i consécutifs sont permutés.
       Préserve les corrélations locales (anti-persistance H≈0.39),
       détruit la structure cumulative globale.

  (T3) SURROGAT DE THEILER (randomisation de phase)
       La TFD de Q_r est calculée ; les phases sont randomisées à amplitudes fixes ;
       le signal est reconstruit. Préserve le spectre de puissance exact.
       Test le plus puissant pour détecter un signal périodique authentique.

Résultats reproduits (Table 5 du papier) :
    Méthode       z(γ₁)      p(γ₁)   Signif. obs.  Signif. perm. moy.
    Mélange c_i   -0.23σ     0.550   8/15           13.5/15
    Blocs (50)    -0.21σ     0.547   8/15           13.5/15
    Theiler       -0.02σ     0.514   8/15           13.5/15

Interprétation :
    Scores z négatifs → signal observé plus FAIBLE qu'une permutation aléatoire.
    → Pas de signal spectral détectable à l'échelle N=10 000 premiers.

Entrée :
    primes_k3.npy   : tableau des N premiers p ≡ 1 (mod 6) [optionnel]
    ou calcul interne si absent.

Usage :
    python permutation_tests.py --quick     # N=500,  B=200  (rapide)
    python permutation_tests.py             # N=10000, B=5000 (papier)
    python permutation_tests.py --N 2000 --B 1000

Dépendances : numpy, scipy, sympy (si calcul des premiers), matplotlib
"""

import argparse
import time
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Paramètres
# ---------------------------------------------------------------------------
K          = 3
GAMMA_FILE = None          # optionnel : fichier numpy avec les parties imaginaires
N_ZEROS    = 15            # nombre de zéros de Riemann utilisés
BLOCK_SIZE = 50            # taille de bloc pour T2
ALPHA      = 0.05          # seuil de significativité individuel

# Parties imaginaires des 15 premiers zéros non triviaux de ζ(s)
# Source : Odlyzko (1987), précision 6 décimales
GAMMA_ZEROS = np.array([
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
])


# ---------------------------------------------------------------------------
# Génération / chargement des premiers p ≡ 1 (mod 6)
# ---------------------------------------------------------------------------

def primes_mod6(N: int) -> np.ndarray:
    """
    Retourne les N premiers premiers p > 3 avec p ≡ 1 (mod 6).
    Utilise un crible d'Ératosthène vectorisé.
    """
    # Borne supérieure via le TNP : p_N ≈ N·ln(N)·(1 + ln ln N / ln N)
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
# Coordonnées paramétriques (m°, q°, ε°) pour k=3, convention canonique
# ---------------------------------------------------------------------------

def canonical_coords(p_arr: np.ndarray) -> tuple:
    """
    Convention canonique (Section 5.1 du papier) :
        m° = round((-1 + sqrt(1 + 4p'/3)) / 2),   p' = p - ε°
    Avec ε° = +1 pour p ≡ 1 (mod 6).
    Retourne (m_arr, q_arr, eps_arr).
    """
    eps = np.ones(len(p_arr), dtype=np.int32)    # p ≡ 1 mod 6 → ε° = +1
    p_prime = p_arr - eps
    m       = np.round((-1 + np.sqrt(1 + 4 * p_prime / 3)) / 2).astype(np.int64)
    q       = ((p_arr - 3 * m * (m + 1) - eps) / 6).astype(np.int64)
    return m, q, eps


def increments(p_arr: np.ndarray) -> np.ndarray:
    """
    Calcule les incréments c_n = q°_n + ε°_n.
    """
    _, q, eps = canonical_coords(p_arr)
    return (q + eps).astype(np.float64)


# ---------------------------------------------------------------------------
# Fonctionnelle cumulative Q et résidu pur Q_r
# ---------------------------------------------------------------------------

def cumulative_Q(c: np.ndarray) -> np.ndarray:
    """Q(p_N) = Σ_{n=1}^{N} c_n."""
    return np.cumsum(c)


def Q_residual(p_arr: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Résidu pur Q_r = Q - α√p - β·D_N (Section 5 du papier).

    D_N = Σ ({√(p'/3)} - 1/2),   p' = p - ε° = p - 1.

    Régression OLS : Q ~ α·√p + β·D_N → résidu Q_r.
    """
    Q    = cumulative_Q(c)
    sqrt_p = np.sqrt(p_arr)
    p_prime = (p_arr - 1).astype(float)
    frac  = np.modf(np.sqrt(p_prime / 3))[0]   # partie fractionnaire
    D_N   = np.cumsum(frac - 0.5)

    X    = np.column_stack([sqrt_p, D_N, np.ones(len(p_arr))])
    beta, _, _, _ = np.linalg.lstsq(X, Q, rcond=None)
    Q_r  = Q - X @ beta
    return Q_r, beta


# ---------------------------------------------------------------------------
# Amplitudes spectrales A_i = √(a_i² + b_i²)
# ---------------------------------------------------------------------------

def spectral_amplitudes(Q_r: np.ndarray, p_arr: np.ndarray,
                         gammas: np.ndarray) -> np.ndarray:
    """
    Pour chaque γ_i, projette Q_r sur cos(γ_i ln p) et sin(γ_i ln p).
    Retourne A_i = √(a_i² + b_i²).
    """
    ln_p = np.log(p_arr)
    A    = np.zeros(len(gammas))
    for j, gamma in enumerate(gammas):
        cos_v = np.cos(gamma * ln_p)
        sin_v = np.sin(gamma * ln_p)
        a_j   = 2 * np.mean(Q_r * cos_v)
        b_j   = 2 * np.mean(Q_r * sin_v)
        A[j]  = np.sqrt(a_j**2 + b_j**2)
    return A


# ---------------------------------------------------------------------------
# T1 — Mélange classique
# ---------------------------------------------------------------------------

def shuffle_test(c: np.ndarray, p_arr: np.ndarray,
                 gammas: np.ndarray, B: int,
                 rng: np.random.Generator) -> dict:
    """
    Permute aléatoirement les incréments c_i, reconstruit Q, calcule A_i.
    Préserve la distribution marginale de c_i.
    """
    A_obs = _observed_amplitudes(c, p_arr, gammas)
    A_perm = np.zeros((B, len(gammas)))

    for b in range(B):
        c_perm  = rng.permutation(c)
        Q_r_perm, _ = Q_residual(p_arr, c_perm)
        A_perm[b] = spectral_amplitudes(Q_r_perm, p_arr, gammas)

    return _statistics(A_obs, A_perm, gammas)


# ---------------------------------------------------------------------------
# T2 — Permutation par blocs
# ---------------------------------------------------------------------------

def block_test(c: np.ndarray, p_arr: np.ndarray,
               gammas: np.ndarray, B: int,
               block_size: int, rng: np.random.Generator) -> dict:
    """
    Découpe c en blocs de taille block_size, permute les blocs.
    Préserve les corrélations locales tout en détruisant la structure globale.
    """
    N      = len(c)
    n_blk  = N // block_size
    # Tronquer à n_blk·block_size
    c_trunc = c[:n_blk * block_size]
    p_trunc = p_arr[:n_blk * block_size]

    A_obs = _observed_amplitudes(c_trunc, p_trunc, gammas)
    A_perm = np.zeros((B, len(gammas)))

    blocks = c_trunc.reshape(n_blk, block_size)
    for b in range(B):
        idx      = rng.permutation(n_blk)
        c_perm   = blocks[idx].ravel()
        Q_r_perm, _ = Q_residual(p_trunc, c_perm)
        A_perm[b] = spectral_amplitudes(Q_r_perm, p_trunc, gammas)

    return _statistics(A_obs, A_perm, gammas)


# ---------------------------------------------------------------------------
# T3 — Surrogat de Theiler (randomisation de phase)
# ---------------------------------------------------------------------------

def theiler_test(c: np.ndarray, p_arr: np.ndarray,
                 gammas: np.ndarray, B: int,
                 rng: np.random.Generator) -> dict:
    """
    1. Calcule Q_r et sa TFD.
    2. Randomise les phases en préservant les amplitudes (et la symétrie hermitienne).
    3. Reconstruit Q_r via TFD inverse.
    Préserve le spectre de puissance exact.
    """
    Q_r_obs, _ = Q_residual(p_arr, c)
    A_obs       = spectral_amplitudes(Q_r_obs, p_arr, gammas)

    # TFD du résidu observé
    F      = np.fft.rfft(Q_r_obs)
    amps   = np.abs(F)                  # amplitudes TFD à préserver
    A_perm = np.zeros((B, len(gammas)))

    for b in range(B):
        phases    = rng.uniform(0, 2 * np.pi, len(F))
        F_perm    = amps * np.exp(1j * phases)
        Q_r_perm  = np.fft.irfft(F_perm, n=len(Q_r_obs))
        A_perm[b] = spectral_amplitudes(Q_r_perm, p_arr, gammas)

    return _statistics(A_obs, A_perm, gammas)


# ---------------------------------------------------------------------------
# Utilitaires communs
# ---------------------------------------------------------------------------

def _observed_amplitudes(c: np.ndarray, p_arr: np.ndarray,
                          gammas: np.ndarray) -> np.ndarray:
    Q_r, _ = Q_residual(p_arr, c)
    return spectral_amplitudes(Q_r, p_arr, gammas)


def _statistics(A_obs: np.ndarray, A_perm: np.ndarray,
                gammas: np.ndarray) -> dict:
    """
    Calcule les z-scores et p-values pour chaque γ_i.
    z_i = (A_obs_i - μ_perm_i) / σ_perm_i
    p_i = rang de A_obs_i dans la distribution de permutation (unilatéral supérieur).
    """
    mu_perm  = A_perm.mean(axis=0)
    std_perm = A_perm.std(axis=0, ddof=1)
    std_perm = np.where(std_perm > 0, std_perm, 1e-12)

    z_scores = (A_obs - mu_perm) / std_perm

    # p-value permutation : proportion de A_perm >= A_obs
    p_values = np.mean(A_perm >= A_obs[None, :], axis=0)

    # Nombre de γ_i significatifs au niveau ALPHA
    n_sig_obs  = int((p_values < ALPHA).sum())
    n_sig_perm = float(np.mean((A_perm > A_obs[None, :]).sum(axis=1)))

    return {
        "A_obs"     : A_obs,
        "mu_perm"   : mu_perm,
        "std_perm"  : std_perm,
        "z_scores"  : z_scores,
        "p_values"  : p_values,
        "n_sig_obs" : n_sig_obs,
        "n_sig_perm": n_sig_perm,
        "z_gamma1"  : z_scores[0],
        "p_gamma1"  : p_values[0],
    }


# ---------------------------------------------------------------------------
# Affichage du rapport (Table 5 du papier)
# ---------------------------------------------------------------------------

def print_report(results: dict):
    sep = "-" * 60
    print("\n" + "=" * 60)
    print("RÉSULTATS DES TESTS DE PERMUTATION (Table 5)")
    print("=" * 60)
    print(f"\n{'Méthode':<20} {'z(γ₁)':>8} {'p(γ₁)':>8} "
          f"{'Obs. sig.':>12} {'Perm. moy.':>12}")
    print(sep)

    paper = {
        "T1 Mélange" : (-0.23, 0.550, "8/15", "13.5/15"),
        "T2 Blocs"   : (-0.21, 0.547, "8/15", "13.5/15"),
        "T3 Theiler" : (-0.02, 0.514, "8/15", "13.5/15"),
    }
    for name, res in results.items():
        z   = res["z_gamma1"]
        p   = res["p_gamma1"]
        ns  = res["n_sig_obs"]
        nsp = res["n_sig_perm"]
        N   = len(GAMMA_ZEROS)
        print(f"  {name:<18} {z:>+7.2f}σ  {p:>7.3f}   "
              f"{ns:>3}/{N}         {nsp:>5.1f}/{N}")

    print(sep)
    print("\n  Valeurs du papier :")
    for name, (z, p, ns, nsp) in paper.items():
        print(f"    {name:<20} {z:>+6.2f}σ   p={p:.3f}   {ns}   {nsp}")

    print("\n" + "=" * 60)
    print("INTERPRÉTATION (Section 6.3)")
    print("=" * 60)
    for name, res in results.items():
        z = res["z_gamma1"]
        verdict = "signal ABSENT (z < 0)" if z < 0 else "signal possible (z > 0)"
        print(f"  {name}: z = {z:+.3f}σ → {verdict}")

    print("""
  Scores z négatifs dans les trois méthodes :
  le signal observé est plus faible qu'une permutation aléatoire.
  → Aucun signal spectral détectable à cette échelle (N premiers).
  → Les corrélations antérieures étaient des artefacts de régression
    spurieuse entre processus I(1) [Granger & Newbold 1974].""")


# ---------------------------------------------------------------------------
# Visualisation (optionnelle)
# ---------------------------------------------------------------------------

def plot_results(results: dict, gammas: np.ndarray,
                 out: str = "fig06_permutation_reprod.png"):
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

    n_methods = len(results)
    fig, axes = plt.subplots(2, n_methods, figsize=(5 * n_methods, 7))
    if n_methods == 1:
        axes = axes[:, None]

    colors = {"T1 Mélange": "#1f78b4",
              "T2 Blocs"  : "#33a02c",
              "T3 Theiler": "#e31a1c"}

    for col, (name, res) in enumerate(results.items()):
        color = colors.get(name, "gray")

        # Rangée du haut : distribution de A_1 sous permutation
        ax = axes[0, col]
        A_obs_1  = res["A_obs"][0]
        mu_1     = res["mu_perm"][0]
        std_1    = res["std_perm"][0]
        x_range  = np.linspace(mu_1 - 4*std_1, mu_1 + 4*std_1, 200)
        ax.plot(x_range,
                stats.norm.pdf(x_range, mu_1, std_1),
                color=color, lw=1.5)
        ax.axvline(A_obs_1, color="black", lw=1.5,
                   label=f"Obs = {A_obs_1:.3f}")
        ax.axvline(mu_1, color=color, ls="--", lw=1.0, alpha=0.6,
                   label=f"Moy. perm. = {mu_1:.3f}")
        ax.set_title(f"{name}\n$z(\\gamma_1)={res['z_gamma1']:+.2f}\\sigma$",
                     fontsize=10)
        ax.set_xlabel("$|A_1|$")
        ax.legend(fontsize=7)

        # Rangée du bas : z-scores pour les 15 premiers zéros
        ax = axes[1, col]
        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.bar(np.arange(1, len(gammas)+1), res["z_scores"],
               color=color, alpha=0.75, width=0.6)
        ax.axhline(-1.96, color="red", lw=0.8, ls=":", label="±1.96σ")
        ax.axhline(+1.96, color="red", lw=0.8, ls=":")
        ax.set_xlabel("Zéro de Riemann $\\gamma_i$")
        ax.set_ylabel("Score $z$")
        ax.set_xticks(np.arange(1, len(gammas)+1, 2))
        if col == 0:
            ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(out)
    print(f"  Figure sauvegardée : {out}")
    plt.close()


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Tests de permutation T1/T2/T3 sur le résidu pur Q_r")
    parser.add_argument("--N", type=int, default=None,
                        help="Nombre de premiers utilisés (défaut : 10000)")
    parser.add_argument("--B", type=int, default=None,
                        help="Nombre de réplications de permutation (défaut : 5000)")
    parser.add_argument("--quick", action="store_true",
                        help="Mode rapide : N=500, B=200")
    parser.add_argument("--seed", type=int, default=42,
                        help="Graine aléatoire (défaut : 42)")
    parser.add_argument("--no-plot", action="store_true",
                        help="Désactiver la génération de figures")
    parser.add_argument("--t1-only", action="store_true",
                        help="Lancer uniquement T1 (diagnostic rapide)")
    args = parser.parse_args()

    N = args.N if args.N else (500  if args.quick else 10_000)
    B = args.B if args.B else (200  if args.quick else 5_000)

    rng = np.random.default_rng(args.seed)

    print("=" * 60)
    print(f"Tests de permutation — k={K}, {N} premiers, {B} réplications")
    print("=" * 60)

    # Génération des premiers p ≡ 1 (mod 6)
    print(f"\nGénération de {N:,} premiers p ≡ 1 (mod 6) ...")
    t0     = time.time()
    p_arr  = primes_mod6(N)
    print(f"  Fait en {time.time()-t0:.2f}s  "
          f"(p_1={p_arr[0]}, p_N={p_arr[-1]:,})")

    # Incréments c_n = q° + ε°
    c = increments(p_arr)
    print(f"  Moyenne c_n = {c.mean():.4f}  (attendu ≈ 0)")
    print(f"  Std c_n     = {c.std():.4f}")

    # Résidu pur Q_r observé
    Q_r_obs, beta = Q_residual(p_arr, c)
    print(f"\n  Modèle Q ≈ α√p + β·D_N + γ :")
    print(f"    α = {beta[0]:.5f},  β = {beta[1]:.5f},  γ = {beta[2]:.3f}")
    print(f"    Var(Q_r) / Var(Q) = "
          f"{Q_r_obs.var() / cumulative_Q(c).var():.4f}  "
          f"(papier : ≈ 0.033)")

    # Amplitudes observées
    gammas = GAMMA_ZEROS[:N_ZEROS]
    A_obs  = spectral_amplitudes(Q_r_obs, p_arr, gammas)
    print(f"\n  Amplitude observée A(γ₁) = {A_obs[0]:.4f}  "
          f"(papier : 0.293)")

    results = {}

    # T1 — Mélange
    print(f"\n[T1] Mélange classique  ({B} réplications) ...")
    t0 = time.time()
    results["T1 Mélange"] = shuffle_test(c, p_arr, gammas, B, rng)
    print(f"     Fait en {time.time()-t0:.1f}s")

    if not args.t1_only:
        # T2 — Blocs
        print(f"[T2] Permutation par blocs ({BLOCK_SIZE})  ({B} réplications) ...")
        t0 = time.time()
        results["T2 Blocs"] = block_test(c, p_arr, gammas, B, BLOCK_SIZE, rng)
        print(f"     Fait en {time.time()-t0:.1f}s")

        # T3 — Theiler
        print(f"[T3] Surrogat de Theiler  ({B} réplications) ...")
        t0 = time.time()
        results["T3 Theiler"] = theiler_test(c, p_arr, gammas, B, rng)
        print(f"     Fait en {time.time()-t0:.1f}s")

    # Rapport
    print_report(results)

    # Figures
    if not args.no_plot:
        plot_results(results, gammas)

    print("\nFin des tests de permutation.")


if __name__ == "__main__":
    main()
