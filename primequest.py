#!/usr/bin/env python3
"""
PrimeQuest — Preuve de primalité pour n'importe quelle taille
Trois méthodes complémentaires, toutes fondées sur Pocklington N-1 :

  1. Famille paramétrique : p = 3m(m+1)+1, m = 2^a·3^b−1
  2. Proth               : p = k·2^n+1,  k impair < 2^n
  3. Pocklington récursif (Maurer) : universel, garantit toute taille

Usage : python3 primequest.py [DIGITS]
        Ou modifier DIGITS_CIBLE ci-dessous.
"""

import gmpy2
import time
import math
import random
import sys
from typing import Optional

# ═══════════════════════════════════════════════════════════════
DIGITS_CIBLE = 1000   # ← modifier ici ou passer en argument
MR_TOURS     = 25     # tours Miller-Rabin (sécurité ≥ 2^-50)
TEMOINS_MAX  = 300    # max témoins Pocklington
TOLERANCE    = 5      # tolérance ±digits pour méthodes 1 & 2
SIEVE_LIMIT  = 100_000
# ═══════════════════════════════════════════════════════════════

deux  = gmpy2.mpz(2)
trois = gmpy2.mpz(3)

# ── Crible d'Ératosthène ─────────────────────────────────────────────────────

def _petits_premiers(limite: int) -> list:
    t = bytearray([1]) * (limite + 1)
    t[0] = t[1] = 0
    for i in range(2, int(limite**0.5) + 1):
        if t[i]:
            t[i*i::i] = bytearray(len(t[i*i::i]))
    return [i for i in range(2, limite + 1) if t[i]]

PREMIERS_PETITS = _petits_premiers(SIEVE_LIMIT)

def crible_rapide(p: gmpy2.mpz) -> bool:
    """False = p est composé (éliminé). True = passe le crible."""
    for q in PREMIERS_PETITS:
        if q * q > p:
            return True
        if p % q == 0:
            return p == q
    return True

# ── Théorème de Pocklington (N-1) ────────────────────────────────────────────

def pocklington_test(p: gmpy2.mpz, facteurs_F: list) -> tuple:
    """
    Théorème de Pocklington : si F | p-1, F > sqrt(p), F complètement factorisé,
    et pour chaque facteur premier q de F on trouve w avec
       w^(p-1) ≡ 1 (mod p)  et  gcd(w^((p-1)/q) - 1, p) = 1,
    alors p est premier.

    facteurs_F : facteurs premiers distincts de F (int ou gmpy2.mpz).
    Retourne (True, {q: témoin}) ou (False, état_partiel).
    """
    p1     = p - 1
    valide = {int(q): None for q in facteurs_F}

    for w in range(2, 2 + TEMOINS_MAX):
        gw = gmpy2.mpz(w)
        if gmpy2.powmod(gw, p1, p) != 1:
            continue
        for q in list(valide):
            if valide[q] is not None:
                continue
            exp = p1 // gmpy2.mpz(q)
            val = gmpy2.powmod(gw, exp, p)
            if gmpy2.gcd(val - 1, p) == 1:
                valide[q] = w
        if all(v is not None for v in valide.values()):
            return True, valide

    return False, valide

# ═══════════════════════════════════════════════════════════════════════════════
# MÉTHODE 1 : Famille paramétrique  p = 3m(m+1)+1,  m = 2^a·3^b−1
# p−1 = 2^a · 3^(b+1) · m  →  F = 2^a · 3^(b+1) complètement factorisé
# ═══════════════════════════════════════════════════════════════════════════════

def _enumerer_paires(digits: int, tol: int) -> list:
    d_min, d_max = digits - tol, digits + tol
    # Borne supérieure sur a : p ≈ 3·(2^a·3^b)^2, donc log10(p) ≈ 2a·log10(2)
    # Pour nb ≤ d_max : a ≤ d_max / (2·log10(2)) ≈ d_max * 1.66
    a_max = int(d_max * 1.67) + 5
    paires = []
    for a in range(1, a_max):
        pow2a = deux ** a
        for b in range(1, a_max):
            m  = pow2a * trois**b - 1
            p  = 3 * m * (m + 1) + 1
            nb = int(gmpy2.num_digits(p))
            if nb < d_min:
                continue
            if nb > d_max:
                break
            F = pow2a * trois**(b + 1)
            paires.append((abs(nb - digits), a, b, m, p, F, nb))
    paires.sort(key=lambda x: x[0])
    return paires

def methode1_famille(digits: int, tol: int) -> Optional[dict]:
    print("\n  ── Méthode 1 : p = 3m(m+1)+1,  m = 2^a·3^b−1 ──")
    t0     = time.perf_counter()
    paires = _enumerer_paires(digits, tol)
    print(f"  {len(paires)} paires (a,b) énumérées  ({time.perf_counter()-t0:.2f}s)")

    if not paires:
        print("  Fenêtre vide — passage à la méthode 2.")
        return None

    for idx, (_, a, b, m, p, F, nb) in enumerate(paires, 1):
        if not crible_rapide(p):
            continue
        if not gmpy2.is_prime(p, MR_TOURS):
            continue
        elapsed = time.perf_counter() - t0
        print(f"  Candidat #{idx} (a={a}, b={b}) — {nb} chiffres "
              f"[{elapsed:.1f}s] — Pocklington…")
        ok, temoins = pocklington_test(p, [2, 3])
        if ok:
            return {"methode": "Famille p=3m(m+1)+1",
                    "p": p, "nb": nb, "a": a, "b": b, "m": m, "F": F,
                    "temoins": temoins}

    print("  Aucun premier prouvé dans cette famille.")
    return None

# ═══════════════════════════════════════════════════════════════════════════════
# MÉTHODE 2 : Proth  p = k·2^n+1,  k impair < 2^n
# Théorème de Proth : ∃a : a^((p-1)/2) ≡ −1 (mod p)  ⟹  p premier
# ═══════════════════════════════════════════════════════════════════════════════

def _proth_test(p: gmpy2.mpz) -> tuple:
    exp = (p - 1) // 2
    for a in range(2, 2 + TEMOINS_MAX):
        r = gmpy2.powmod(gmpy2.mpz(a), exp, p)
        if r == p - 1:
            return True, a
    return False, None

def methode2_proth(digits: int, tol: int) -> Optional[dict]:
    print("\n  ── Méthode 2 : Proth  p = k·2^n+1 ──")
    d_min, d_max = digits - tol, digits + tol
    # n ≈ (digits/2)·log2(10)  →  2^n a ≈ digits/2 chiffres
    n_cible   = int(digits / 2 * math.log2(10))
    MAX_ESSAIS = 200_000
    essais    = 0

    for delta_n in range(0, 40):
        for sign in [0, -1, 1]:
            n = n_cible + sign * delta_n
            if n < 1:
                continue
            pow2n = deux ** n
            k_min = max(gmpy2.mpz(1), (gmpy2.mpz(10)**(d_min - 1) - 1) // pow2n)
            k_max = (gmpy2.mpz(10)**d_max - 1) // pow2n
            if k_min > k_max or k_min >= pow2n:
                continue

            tentatives = min(3000, int(k_max - k_min) + 1)
            for _ in range(tentatives):
                k = gmpy2.mpz(random.randint(int(k_min), int(k_max)))
                if k % 2 == 0:
                    k += 1
                if k >= pow2n:
                    continue
                p  = k * pow2n + 1
                nb = int(gmpy2.num_digits(p))
                if not (d_min <= nb <= d_max):
                    continue
                essais += 1
                if essais > MAX_ESSAIS:
                    print(f"  Limite d'essais atteinte ({MAX_ESSAIS}).")
                    return None
                if not crible_rapide(p):
                    continue
                if not gmpy2.is_prime(p, MR_TOURS):
                    continue
                print(f"  Candidat Proth k={k}, n={n} — {nb} chiffres — test Proth…")
                ok, a_t = _proth_test(p)
                if ok:
                    return {"methode": "Proth k·2^n+1",
                            "p": p, "nb": nb, "k": k, "n": n, "a_temoin": a_t}

    print(f"  Aucun Proth prouvé ({essais} essais).")
    return None

# ═══════════════════════════════════════════════════════════════════════════════
# MÉTHODE 3 : Pocklington récursif (algorithme de Maurer)
# Garantit la certification pour n'importe quelle taille
#
# Principe : pour certifier p à bits_cible bits,
#   1. Certifier récursivement q à ≈ bits_cible/2 bits
#   2. Chercher p = 2·k·q + 1 (donc p−1 = 2·k·q, F = 2·q > √p)
#   3. Appliquer Pocklington avec facteurs {2, q}
# Profondeur ≈ log2(bits_cible / 64) niveaux
# ═══════════════════════════════════════════════════════════════════════════════

def _certifier_petit(bits: int) -> gmpy2.mpz:
    """Premier certifié de ~bits bits (cas de base, bits ≤ 64)."""
    if bits <= 2:
        return gmpy2.mpz(3)
    lo = gmpy2.mpz(1) << (bits - 1)
    hi = (gmpy2.mpz(1) << bits) - 1
    p  = lo | gmpy2.mpz(1)
    while p <= hi:
        if gmpy2.is_prime(p, 30):
            return p
        p += 2
    return gmpy2.mpz(hi | 1)

def pocklington_recursif(bits_cible: int, profondeur: int = 0) -> tuple:
    """
    Retourne (p, chaine) :
      p     : premier certifié de exactement bits_cible bits
      chaine: liste de dicts décrivant chaque niveau de preuve
    """
    ind = "    " + "  " * profondeur

    if bits_cible <= 64:
        p = _certifier_petit(bits_cible)
        print(f"{ind}Base [{profondeur}] : {bits_cible} bits → {p}")
        return p, [{"type": "base", "p": p}]

    # bits du sous-premier q : bits_cible//2 + 3 garantit F = 2q > √p
    bits_q = bits_cible // 2 + 3

    print(f"{ind}Niveau [{profondeur}] : {bits_cible} bits  "
          f"(sous-premier : {bits_q} bits)…")

    for tentative in range(10):
        q, chaine_q = pocklington_recursif(bits_q, profondeur + 1)

        lo    = gmpy2.mpz(1) << (bits_cible - 1)
        hi    = (gmpy2.mpz(1) << bits_cible) - 1
        deux_q = deux * q
        k_min  = max(gmpy2.mpz(1), (lo - 1) // deux_q)
        k_max  = hi // deux_q

        if k_min > k_max:
            print(f"{ind}  Plage vide, nouveau sous-premier…")
            continue

        found = False
        for _ in range(300 * bits_cible):
            k = gmpy2.mpz(random.randint(int(k_min), int(k_max)))
            p = deux_q * k + 1

            if not crible_rapide(p):
                continue
            if not gmpy2.is_prime(p, MR_TOURS):
                continue

            F = deux_q
            if F * F < p:
                continue

            ok, temoins = pocklington_test(p, [2, int(q)])
            if ok:
                nb = int(gmpy2.num_digits(p))
                print(f"{ind}  ✓ p ({nb} chiffres) prouvé  "
                      f"[témoins q=2:{temoins[2]}, q=large:{temoins[int(q)]}]")
                chaine = chaine_q + [{"type": "pocklington",
                                       "p": p, "q": q, "k": k,
                                       "F": F, "temoins": temoins}]
                return p, chaine
            found = True

        if not found:
            print(f"{ind}  Aucun candidat,  nouveau sous-premier…")

    raise RuntimeError(f"Impossible de certifier {bits_cible} bits après 10 tentatives")

def methode3_recursive(digits: int) -> Optional[dict]:
    print("\n  ── Méthode 3 : Pocklington récursif (Maurer) — universel ──")
    bits = int(digits * math.log2(10)) + 1
    profondeur_max = math.ceil(math.log2(max(bits / 64, 1))) + 1
    print(f"  Objectif : ~{bits} bits ≈ {digits} chiffres")
    print(f"  Profondeur de récursion estimée : {profondeur_max} niveaux")

    try:
        p, chaine = pocklington_recursif(bits)
    except RuntimeError as e:
        print(f"  Erreur : {e}")
        return None

    nb = int(gmpy2.num_digits(p))
    return {"methode": "Pocklington récursif (Maurer)",
            "p": p, "nb": nb, "chaine": chaine}

# ═══════════════════════════════════════════════════════════════════════════════
# Affichage & sauvegarde
# ═══════════════════════════════════════════════════════════════════════════════

def afficher_resultat(res: dict, t_total: float):
    p  = res["p"]
    nb = res["nb"]
    sp = str(p)
    extrait = f"{sp[:40]}…{sp[-20:]}" if len(sp) > 62 else sp

    print(f"\n{'═'*65}")
    print(f"  PREMIER DE {nb} CHIFFRES — PRIMALITÉ PROUVÉE")
    print(f"{'═'*65}")
    print(f"  Méthode : {res['methode']}")
    print(f"  Temps   : {t_total:.2f} sec\n")

    m = res["methode"]
    if m == "Famille p=3m(m+1)+1":
        a, b, fm, F = res["a"], res["b"], res["m"], res["F"]
        t = res["temoins"]
        print(f"  a = {a},  b = {b}")
        print(f"  m = 2^{a}·3^{b}−1   ({gmpy2.num_digits(fm)} chiffres)")
        print(f"  p = 3m(m+1)+1    ({nb} chiffres)")
        print(f"  F = 2^{a}·3^{{{b+1}}} ({gmpy2.num_digits(F)} chiffres) > √p  ✅")
        print(f"  Témoins Pocklington :")
        print(f"    q=2  → w={t[2]}  ✅")
        print(f"    q=3  → w={t[3]}  ✅")

    elif m == "Proth k·2^n+1":
        print(f"  k = {res['k']}")
        print(f"  n = {res['n']}")
        print(f"  p = k·2^{res['n']}+1  ({nb} chiffres)")
        print(f"  Théorème de Proth : a={res['a_temoin']}  ✅")
        print(f"    a^((p-1)/2) ≡ -1 (mod p)  ✅")

    elif "récursif" in m:
        chaine = res["chaine"]
        print(f"  Longueur de la chaîne de certification : {len(chaine)} niveaux")
        print(f"  Résumé (du plus petit au plus grand) :")
        for i, c in enumerate(chaine):
            if c["type"] == "base":
                print(f"    [{i}] Base      : p = {c['p']}")
            else:
                pn = int(gmpy2.num_digits(c["p"]))
                qn = int(gmpy2.num_digits(c["q"]))
                print(f"    [{i}] Pocklington: {pn} chiffres = 2·{c['k']}·q({qn} ch)+1")

    print(f"\n  p = {extrait}")

def sauvegarder(res: dict, t_total: float):
    p  = res["p"]
    nb = res["nb"]
    nom = f"premier_{nb}chiffres.txt"
    with open(nom, "w") as f:
        f.write(f"Premier de {nb} chiffres\n")
        f.write(f"Méthode  : {res['methode']}\n")
        f.write(f"Temps    : {t_total:.2f}s\n\n")
        if res["methode"] == "Famille p=3m(m+1)+1":
            a, b = res["a"], res["b"]
            t    = res["temoins"]
            f.write(f"a={a}, b={b}\n")
            f.write(f"Témoins Pocklington : q=2→w={t[2]}, q=3→w={t[3]}\n\n")
            f.write(f"m =\n{res['m']}\n\n")
        elif res["methode"] == "Proth k·2^n+1":
            f.write(f"k={res['k']}, n={res['n']}\n")
            f.write(f"Témoin Proth : a={res['a_temoin']}\n\n")
        elif "récursif" in res["methode"]:
            chaine = res["chaine"]
            f.write(f"Chaîne de certification ({len(chaine)} niveaux) :\n")
            for i, c in enumerate(chaine):
                if c["type"] == "base":
                    f.write(f"  [{i}] base : {c['p']}\n")
                else:
                    pn = int(gmpy2.num_digits(c["p"]))
                    qn = int(gmpy2.num_digits(c["q"]))
                    f.write(f"  [{i}] pocklington : {pn} chiffres, q={qn} chiffres, k={c['k']}\n")
            f.write("\n")
        f.write(f"p =\n{p}\n")
    print(f"\n  Sauvegardé → {nom}")

# ═══════════════════════════════════════════════════════════════════════════════
# Programme principal
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            DIGITS_CIBLE = int(sys.argv[1])
        except ValueError:
            print(f"Usage : {sys.argv[0]} [nombre_de_chiffres]")
            sys.exit(1)

    print("=" * 65)
    print(f"  PrimeQuest — Preuve de primalité  ≈{DIGITS_CIBLE} chiffres")
    print("=" * 65)
    print(f"  Crible     : jusqu'à {SIEVE_LIMIT} ({len(PREMIERS_PETITS)} premiers)")
    print(f"  Tolérance  : ±{TOLERANCE} chiffres pour méthodes 1 & 2")
    print(f"  Miller-Rabin : {MR_TOURS} tours")

    t_debut = time.perf_counter()
    res     = None

    # Méthode 1 : famille rapide avec structure connue
    res = methode1_famille(DIGITS_CIBLE, TOLERANCE)

    # Méthode 2 : Proth si la famille ne donne rien
    if res is None:
        res = methode2_proth(DIGITS_CIBLE, TOLERANCE)

    # Méthode 3 : Pocklington récursif — garantie universelle
    if res is None:
        res = methode3_recursive(DIGITS_CIBLE)

    t_total = time.perf_counter() - t_debut

    if res is None:
        print("\n  ❌ Échec de toutes les méthodes.")
        sys.exit(1)

    afficher_resultat(res, t_total)
    sauvegarder(res, t_total)
