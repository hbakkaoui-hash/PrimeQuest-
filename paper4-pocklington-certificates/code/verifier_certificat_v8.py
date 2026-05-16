#!/usr/bin/env python3
"""
Vérification indépendante du certificat Pocklington — PrimeQuest v8 (k=23)
───────────────────────────────────────────────────────────────────────────
Recompute p from scratch and checks all Pocklington conditions.

Usage : py verifier_certificat_v8.py
"""

import gmpy2
import math
import time

# ── Paramètres du certificat à vérifier ─────────────────────────────────────
A  = 24_894
B  =  5_525
W2 = 3       # témoin pour q=2
W23 = 2      # témoin pour q=23
K  = 23

print("=" * 65)
print("  Vérification certificat Pocklington — PrimeQuest v8 (k=23)")
print("=" * 65)
print(f"  a={A}, b={B}, k={K}")
print(f"  Témoins : w(q=2)={W2}, w(q=23)={W23}")
print()

t0 = time.perf_counter()

# ── Reconstruction de m et p ─────────────────────────────────────────────────
print("  1. Calcul de m = 2^a · 23^b − 1 …", end=" ", flush=True)
deux = gmpy2.mpz(2)
kk   = gmpy2.mpz(K)
m = deux**A * kk**B - 1
print(f"OK  ({gmpy2.num_digits(m)} chiffres)")

print("  2. Calcul de p = 23·m·(m+1)+1 …", end=" ", flush=True)
p = K * m * (m + 1) + 1
nb = int(gmpy2.num_digits(p))
print(f"OK  ({nb} chiffres)")

# ── Vérification nombre de chiffres ─────────────────────────────────────────
assert nb == 30_037, f"Nombre de chiffres inattendu : {nb}"
print(f"  3. Nombre de chiffres : {nb}  ✅")

# ── Vérification F > √p  [RIGOUREUSE, INCONDITIONNELLE] ─────────────────────
print("  4. Condition F > √p …", end=" ", flush=True)
F = deux**A * kk**(B + 1)
assert F * F > p, "ÉCHEC : F² ≤ p !"
ratio = float(gmpy2.log2(F)) / float(gmpy2.log2(p)) * 2
print(f"OK  F/√p ≥ √23 ≈ 4.796  ✅")

# ── Vérification p−1 = F · m ─────────────────────────────────────────────────
print("  5. Vérification p−1 = F·m …", end=" ", flush=True)
assert p - 1 == F * m, "ÉCHEC : p−1 ≠ F·m !"
print("OK  ✅")

# ── Témoin w=W2 pour q=2 ─────────────────────────────────────────────────────
print(f"  6. Témoin w={W2} pour q=2 …", end=" ", flush=True)
p1 = p - 1
gw = gmpy2.mpz(W2)
assert gmpy2.powmod(gw, p1, p) == 1, f"ÉCHEC : {W2}^(p−1) ≢ 1 mod p"
val2 = gmpy2.powmod(gw, p1 // 2, p)
assert gmpy2.gcd(val2 - 1, p) == 1, f"ÉCHEC : pgcd({W2}^((p−1)/2)−1, p) ≠ 1"
print(f"OK  ✅")

# ── Témoin w=W23 pour q=23 ───────────────────────────────────────────────────
print(f"  7. Témoin w={W23} pour q=23 …", end=" ", flush=True)
gw23 = gmpy2.mpz(W23)
assert gmpy2.powmod(gw23, p1, p) == 1, f"ÉCHEC : {W23}^(p−1) ≢ 1 mod p"
val23 = gmpy2.powmod(gw23, p1 // 23, p)
assert gmpy2.gcd(val23 - 1, p) == 1, f"ÉCHEC : pgcd({W23}^((p−1)/23)−1, p) ≠ 1"
print(f"OK  ✅")

# ── Résultat ─────────────────────────────────────────────────────────────────
dt = time.perf_counter() - t0
sp = str(p)
print()
print("=" * 65)
print(f"  p EST UN NOMBRE PREMIER  [PREUVE DÉTERMINISTE]  ✅")
print(f"  Vérification indépendante en {dt:.1f}s")
print(f"  p = {sp[:40]}…{sp[-20:]}")
print("=" * 65)
