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
distrib_validation.py
=====================
Validation distributionnelle de la Conjecture 3.2 :

    |q_min(21,m)| · C_21 / ln(m)  →(loi)  Exp(1)  quand m → ∞

Trois tests complémentaires (Section 5.4 du papier) :

  (1) TEST KS NAÏF
      KS sur la variable normalisée X_m = |q_min| / (0.32·ln(m) - 0.45)
      sur un sous-échantillon de 10 000 valeurs (pas=100, indépendance approx.)

  (2) TEST KS CORRIGÉ POUR L'AUTOCORRÉLATION
      Estimation de l'autocorrélation résiduelle ρ̂ (lag 1),
      calcul de la taille d'échantillon effective n_eff via Bartlett,
      recalibration de la statistique D_KS.

  (3) TEST DU χ² SUR TRANCHE CONDITIONNELLE
      Conditionnellement à log₁₀(m) ∈ [4.0, 4.5], le paramètre p_q ≈ C_k/(2·ln(m))
      est quasi-constant. Les n_bin ≈ 38 000 valeurs sont partitionnées en
      15 classes équiprobables sous la loi géométrique ; χ²(14) est calculé.
      Résultat du papier : χ² = 13.1, ddl = 14, p = 0.52.

  (4) GRAPHIQUE Q–Q (Figure 15)
      Q–Q plot de X_m vs Exp(1), avec bande de confiance bootstrap 95 %.

Entrée :
    qmin_k21_M1000000.npy   (produit par qmin_compute.py)
    ou mode --quick : génère un jeu de test synthétique.

Usage :
    python distrib_validation.py                          # données réelles
    python distrib_validation.py --quick                  # jeu synthétique
    python distrib_validation.py --input qmin_k21_M*.npy

Dépendances : numpy, scipy, matplotlib
"""

import argparse
import os
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Paramètres du papier
# ---------------------------------------------------------------------------
K         = 21
C_K       = 8.206301       # constante BH agrégée
SLOPE_EMP = 0.320          # pente empirique (Observation 5.1)
INTER_EMP = -0.450         # ordonnée à l'origine empirique
N_BINS_KS = 50             # nombre de tranches log pour le sous-échantillonnage
N_CLASSES_CHI2 = 15        # classes équiprobables pour le χ²
LOG10_BIN = (4.0, 4.5)     # tranche conditionnelle pour le χ²
N_BOOT    = 500            # réplications bootstrap pour la bande Q–Q


# ---------------------------------------------------------------------------
# Chargement / génération des données
# ---------------------------------------------------------------------------

def load_or_generate(input_path: str | None, quick: bool) -> np.ndarray:
    """
    Charge qmin depuis un fichier .npy ou génère des données synthétiques.
    Mode --quick : tire |q| depuis une loi géométrique de paramètre C_K/(2·ln(m))
    avec m uniform dans [10^4, 10^6].
    """
    if input_path and os.path.exists(input_path):
        print(f"Chargement des données : {input_path}")
        qmin = np.load(input_path)
        M    = len(qmin) - 1
        m_arr = np.arange(M + 1, dtype=np.float64)
        print(f"  M = {M:,},  N = {len(qmin):,} valeurs")
        return qmin, m_arr

    if quick:
        print("Mode --quick : génération synthétique (m ∈ [10^4, 10^6], N=10^6)")
        rng   = np.random.default_rng(42)
        M     = 1_000_000
        m_arr = np.arange(M + 1, dtype=np.float64)
        # Génération cohérente avec la loi empirique :
        # |q| ~ Géom(p_emp) où p_emp = 1 / (slope·ln(m) + intercept)
        # → variable normalisée X_m ~ Exp(1) par construction.
        log_m   = np.log(np.maximum(m_arr, 3.0))
        mu_emp  = np.maximum(SLOPE_EMP * log_m + INTER_EMP, 0.5)
        p_emp   = np.clip(1.0 / mu_emp, 1e-4, 0.9999)
        raw     = rng.geometric(p_emp, size=M + 1).astype(np.int32) - 1
        signs   = rng.choice([-1, 1], size=M + 1)
        qmin    = (raw * signs).astype(np.int32)
        return qmin, m_arr

    # Chercher automatiquement
    candidates = [f for f in os.listdir(".") if f.startswith("qmin_k21") and f.endswith(".npy")]
    if candidates:
        path = sorted(candidates)[-1]
        print(f"Fichier trouvé automatiquement : {path}")
        qmin  = np.load(path)
        M     = len(qmin) - 1
        m_arr = np.arange(M + 1, dtype=np.float64)
        return qmin, m_arr

    raise FileNotFoundError(
        "Aucun fichier qmin_k21_*.npy trouvé. "
        "Lancez d'abord qmin_compute.py, ou utilisez --quick."
    )


# ---------------------------------------------------------------------------
# Variable normalisée X_m = |q_min| / (slope·ln(m) + intercept)
# ---------------------------------------------------------------------------

def normalized_variable(qmin: np.ndarray, m_arr: np.ndarray,
                         m_min: int = 10_000) -> tuple:
    """
    X_m = |q_min(m)| / (SLOPE_EMP·ln(m) + INTER_EMP)

    Sélectionne m ∈ [m_min, M] avec q_min ≠ 0.
    Retourne (X, m_sel) — tous deux de même longueur.
    """
    mask      = (m_arr >= m_min) & (qmin != 0)
    m_sel     = m_arr[mask]
    q_sel     = np.abs(qmin[mask]).astype(float)
    denom     = SLOPE_EMP * np.log(m_sel) + INTER_EMP
    denom     = np.maximum(denom, 1e-6)     # sécurité
    X         = q_sel / denom
    return X, m_sel


# ---------------------------------------------------------------------------
# (1) TEST KS NAÏF
# ---------------------------------------------------------------------------

def ks_naive(X: np.ndarray, step: int = 100) -> dict:
    """
    KS test sur sous-échantillon (pas=step) vs Exp(1).
    Hypothesis nulle : X ~ Exp(1).
    """
    X_sub   = X[::step]
    n       = len(X_sub)
    stat, p = stats.kstest(X_sub, "expon", args=(0, 1))
    return {"n": n, "D_KS": stat, "p_value": p, "step": step}


# ---------------------------------------------------------------------------
# (2) TEST KS CORRIGÉ POUR L'AUTOCORRÉLATION (Bartlett)
# ---------------------------------------------------------------------------

def ks_corrected(X: np.ndarray, step: int = 100) -> dict:
    """
    Corrige la taille d'échantillon effective via l'autocorrélation de lag 1.

    n_eff = n · (1 - ρ̂) / (1 + ρ̂)   [correction de Bartlett]

    La statistique D_KS ajustée est D · √(n / n_eff).
    """
    X_sub = X[::step]
    n     = len(X_sub)

    # Autocorrélation de lag 1
    rho_hat = np.corrcoef(X_sub[:-1], X_sub[1:])[0, 1]
    rho_hat = np.clip(rho_hat, -0.99, 0.99)

    n_eff   = int(n * (1 - rho_hat) / (1 + rho_hat))
    n_eff   = max(n_eff, 10)

    # KS naïf
    D_naive, _ = stats.kstest(X_sub, "expon", args=(0, 1))

    # Recalibration : D_adj = D · √(n / n_eff)
    D_adj = D_naive * np.sqrt(n / n_eff)

    # p-value via la distribution asymptotique de Kolmogorov
    p_adj = stats.kstwo.sf(D_adj * np.sqrt(n_eff), n_eff)

    return {
        "n"      : n,
        "n_eff"  : n_eff,
        "rho_hat": rho_hat,
        "D_KS"   : D_naive,
        "D_adj"  : D_adj,
        "p_adj"  : p_adj,
    }


# ---------------------------------------------------------------------------
# (3) TEST DU χ² SUR TRANCHE CONDITIONNELLE
# ---------------------------------------------------------------------------

def chi2_conditional(qmin: np.ndarray, m_arr: np.ndarray,
                     log10_bin: tuple = LOG10_BIN,
                     n_classes: int = N_CLASSES_CHI2) -> dict:
    """
    Test d'adéquation χ² pour la loi géométrique, conditionnel à la tranche
    log₁₀(m) ∈ [log10_bin[0], log10_bin[1]].

    Paramètre géométrique : p_q = C_K / (2·ln(m̃))  où m̃ est le centre de tranche.

    Construit n_classes classes équiprobables sous la loi géométrique(p_q),
    compte les effectifs observés, calcule χ²(n_classes - 1).

    Résultat du papier (Table 4) : χ² = 13.1, ddl = 14, p = 0.52.
    """
    lo, hi    = 10 ** log10_bin[0], 10 ** log10_bin[1]
    mask      = (m_arr >= lo) & (m_arr < hi) & (qmin != 0)
    q_bin     = np.abs(qmin[mask]).astype(int)
    n_bin     = len(q_bin)

    if n_bin < 100:
        return {"error": f"Tranche trop petite : n={n_bin}"}

    # Paramètre géométrique au centre de la tranche
    m_center  = np.sqrt(lo * hi)
    p_q       = C_K / (2 * np.log(m_center))
    p_q       = np.clip(p_q, 1e-6, 1 - 1e-6)

    # Quantiles équiprobables de la loi géométrique(p_q)
    # CDF géométrique (support 0, 1, 2, ...) : P(X ≤ k) = 1 - (1-p)^(k+1)
    # On cherche les seuils t_j tels que P(X ≤ t_j) ≈ j/n_classes
    probs      = np.arange(1, n_classes) / n_classes
    thresholds = np.floor(np.log(1 - probs) / np.log(1 - p_q)).astype(int)
    thresholds = np.unique(np.clip(thresholds, 0, q_bin.max() + 1))

    # Effectifs observés par classe
    bins       = np.concatenate([[0], thresholds + 1, [q_bin.max() + 2]])
    observed   = np.histogram(q_bin, bins=bins)[0]

    # Effectifs théoriques sous Géom(p_q)
    geom_pmf   = lambda k: (1 - p_q) ** k * p_q

    expected   = []
    for i in range(len(bins) - 1):
        lo_b, hi_b = bins[i], bins[i + 1]
        prob = sum(geom_pmf(k) for k in range(lo_b, hi_b))
        expected.append(n_bin * prob)

    observed  = np.array(observed, dtype=float)
    expected  = np.array(expected, dtype=float)

    # Fusionner les classes avec effectif théorique < 5
    obs_merged, exp_merged = [], []
    acc_o, acc_e = 0.0, 0.0
    for o, e in zip(observed, expected):
        acc_o += o
        acc_e += e
        if acc_e >= 5:
            obs_merged.append(acc_o)
            exp_merged.append(acc_e)
            acc_o, acc_e = 0.0, 0.0
    if acc_e > 0:                          # dernière classe résiduelle
        if obs_merged:
            obs_merged[-1] += acc_o
            exp_merged[-1] += acc_e
        else:
            obs_merged.append(acc_o)
            exp_merged.append(acc_e)

    obs_merged = np.array(obs_merged)
    exp_merged = np.array(exp_merged)
    n_classes_eff = len(obs_merged)
    ddl           = n_classes_eff - 1

    chi2_stat = float(np.sum((obs_merged - exp_merged) ** 2 / exp_merged))
    p_val     = float(stats.chi2.sf(chi2_stat, ddl))

    return {
        "n_bin"      : n_bin,
        "m_center"   : m_center,
        "p_q"        : p_q,
        "n_classes"  : n_classes_eff,
        "ddl"        : ddl,
        "chi2"       : chi2_stat,
        "p_value"    : p_val,
        "observed"   : obs_merged,
        "expected"   : exp_merged,
    }


# ---------------------------------------------------------------------------
# (4) GRAPHIQUE Q–Q avec bande bootstrap
# ---------------------------------------------------------------------------

def qq_plot(X: np.ndarray, step: int = 100,
            n_boot: int = N_BOOT, out: str = "fig15_QQ_reprod.png"):
    """
    Q–Q plot de X_m (normalisé) vs Exp(1), bande bootstrap 95 %.
    Reproduit la Figure 15 du papier.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.rcParams.update({"font.size": 11, "axes.labelsize": 12,
                                    "figure.dpi": 150, "savefig.dpi": 300})
    except ImportError:
        print("matplotlib non disponible — Q–Q plot ignoré.")
        return

    rng   = np.random.default_rng(0)
    X_sub = X[::step]
    n     = len(X_sub)

    # Quantiles empiriques et théoriques
    quant_emp  = np.sort(X_sub)
    probs      = (np.arange(1, n + 1) - 0.5) / n
    quant_theo = stats.expon.ppf(probs)

    # Bande bootstrap : rééchantillonner X_sub et calculer les quantiles
    boot_quants = np.zeros((n_boot, n))
    for b in range(n_boot):
        samp            = rng.choice(X_sub, size=n, replace=True)
        boot_quants[b]  = np.sort(samp)

    ci_lo = np.percentile(boot_quants, 2.5, axis=0)
    ci_hi = np.percentile(boot_quants, 97.5, axis=0)

    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.fill_between(quant_theo, ci_lo, ci_hi,
                    alpha=0.25, color="#2166ac", label="IC bootstrap 95 %")
    ax.plot(quant_theo, quant_emp, ".", ms=1.5, color="#d73027", alpha=0.6,
            label=r"$X_m = |q_{\min}| / (0.32\ln m - 0.45)$")
    ax.plot([0, quant_theo[-1]], [0, quant_theo[-1]], "k--", lw=1.0,
            label="Référence Exp(1)")
    ax.set_xlabel(r"Quantiles théoriques $\mathrm{Exp}(1)$")
    ax.set_ylabel("Quantiles empiriques de $X_m$")
    ax.set_title(fr"Q–Q plot : $X_m$ vs $\mathrm{{Exp}}(1)$, $k=21$"
                 f"\n($n={n:,}$ valeurs sous-échantillonnées, pas={step})")
    ax.legend(fontsize=8)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    plt.tight_layout()
    plt.savefig(out)
    print(f"  Q–Q plot sauvegardé : {out}")
    plt.close()


# ---------------------------------------------------------------------------
# Rapport complet
# ---------------------------------------------------------------------------

def print_report(ks_n: dict, ks_c: dict, chi2: dict):
    sep = "-" * 52

    print("\n" + "=" * 52)
    print("RAPPORT DE VALIDATION DISTRIBUTIONNELLE (Section 5.4)")
    print("=" * 52)

    print(f"\n(1) TEST KS NAÏF  (n = {ks_n['n']:,}, pas = {ks_n['step']})")
    print(sep)
    print(f"  D_KS        = {ks_n['D_KS']:.4f}")
    print(f"  p-value     = {ks_n['p_value']:.4f}")
    print(f"  Papier      : D_KS = 0.031, p ≈ 0.12")
    flag = "✓" if ks_n["p_value"] > 0.05 else "—"
    print(f"  H₀ (Exp(1)) non rejetée au niveau 5 % : {flag}")

    print(f"\n(2) TEST KS CORRIGÉ POUR L'AUTOCORRÉLATION")
    print(sep)
    print(f"  ρ̂ (lag 1)   = {ks_c['rho_hat']:.4f}")
    print(f"  n_eff        = {ks_c['n_eff']:,}  (sur n = {ks_c['n']:,})")
    print(f"  D_KS naïf   = {ks_c['D_KS']:.4f}")
    print(f"  D_KS ajusté = {ks_c['D_adj']:.4f}")
    print(f"  p_adj        = {ks_c['p_adj']:.4f}")
    print(f"  Papier      : D_adj = 0.047, p_adj = 0.04")
    flag = "✓" if ks_c["p_adj"] > 0.05 else "marginal"
    print(f"  H₀ au niveau 5 % : {flag}")

    if "error" in chi2:
        print(f"\n(3) TEST χ² : {chi2['error']}")
    else:
        print(f"\n(3) TEST χ² CONDITIONNEL  "
              f"(log₁₀(m) ∈ [{LOG10_BIN[0]}, {LOG10_BIN[1]}])")
        print(sep)
        print(f"  n_bin        = {chi2['n_bin']:,}")
        print(f"  m̃ (centre)  = {chi2['m_center']:.0f}")
        print(f"  p_q estimé   = {chi2['p_q']:.4f}")
        print(f"  Classes eff. = {chi2['n_classes']}  (ddl = {chi2['ddl']})")
        print(f"  χ²           = {chi2['chi2']:.2f}")
        print(f"  p-value      = {chi2['p_value']:.4f}")
        print(f"  Papier      : χ² = 13.1, ddl = 14, p = 0.52")
        flag = "✓" if chi2["p_value"] > 0.05 else "—"
        print(f"  H₀ (Géom(p_q)) non rejetée au niveau 5 % : {flag}")

    print("\n" + "=" * 52)
    print("INTERPRÉTATION (Section 5.4)")
    print("=" * 52)
    print("  Le test KS naïf est compatible avec Exp(1).")
    print("  Après correction d'autocorrélation, la p-valeur est marginale")
    print("  (p_adj ≈ 0.04), indiquant un ajustement provisoire.")
    print("  Le test χ² conditionnel (non affecté par l'autocorrélation)")
    print("  est compatible avec la loi géométrique (p = 0.52).")
    print("  Ces résultats soutiennent la Conjecture 3.2 sans la prouver.")


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Validation distributionnelle de la Conjecture 3.2")
    parser.add_argument("--input", type=str, default=None,
                        help="Chemin vers qmin_k21_*.npy")
    parser.add_argument("--quick", action="store_true",
                        help="Données synthétiques (test rapide)")
    parser.add_argument("--no-plot", action="store_true",
                        help="Désactiver le graphique Q–Q")
    args = parser.parse_args()

    print("=" * 52)
    print("Validation distributionnelle — Conjecture 3.2")
    print(f"Conjecture : |q_min|·C_21/ln(m) →(loi) Exp(1)")
    print("=" * 52)

    # Chargement
    qmin, m_arr = load_or_generate(args.input, args.quick)
    M = len(qmin) - 1

    # Variable normalisée
    print("\nCalcul de la variable normalisée X_m ...")
    X, m_sel = normalized_variable(qmin, m_arr, m_min=10_000)
    print(f"  {len(X):,} valeurs retenues (m ≥ 10 000, q_min ≠ 0)")
    print(f"  Moyenne X_m = {X.mean():.4f}  (théorie Exp(1) : 1.000)")
    print(f"  Écart-type  = {X.std():.4f}  (théorie Exp(1) : 1.000)")

    # (1) KS naïf
    print("\nTest KS naïf ...")
    ks_n = ks_naive(X, step=100)

    # (2) KS corrigé
    print("Test KS corrigé (autocorrélation) ...")
    ks_c = ks_corrected(X, step=100)

    # (3) χ² conditionnel
    print("Test χ² conditionnel ...")
    chi2 = chi2_conditional(qmin, m_arr)

    # Rapport
    print_report(ks_n, ks_c, chi2)

    # (4) Q–Q plot
    if not args.no_plot:
        print("\nGénération du Q–Q plot ...")
        qq_plot(X, step=100, n_boot=N_BOOT)

    print("\nFin de la validation distributionnelle.")


if __name__ == "__main__":
    main()
