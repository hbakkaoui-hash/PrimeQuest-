import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, floor, log10
from sympy import isprime
import pandas as pd

# Fonction pour générer les nombres premiers jusqu'à une limite donnée avec le crible d'Ératosthène
def crible_eratosthene(n):
    premiers = []
    nombres = [True] * (n - 1)
    for i in range(2, int(sqrt(n)) + 1):
        if nombres[i - 2]:
            premiers.append(i)
            for j in range(i * i, n + 1, i):
                nombres[j - 2] = False
    for i in range(int(sqrt(n)) + 1, n + 1):
        if nombres[i - 2]:
            premiers.append(i)
    return premiers

# Fonction pour calculer le nombre de chiffres d'un entier
def nombre_de_chiffres(n):
    return floor(log10(n)) + 1 if n > 0 else 1

# Paramètres
k = 21  # Facteur fixe
n_max = 10000  # n varie de 5001 à 10000
q_min_range = -50  # Minimum de q
q_max_range = 50  # Maximum de q

# Liste pour stocker les résultats
resultats = []
points = []

# Recherche des nombres premiers
for n in range(n_max + 1):  # n varie de 0 à n_max
    p_n = k * n * (n + 1)
    q_min = None
    p_min = None

    # Recherche de q
    for step in range(q_max_range + 1):
        for q in (step, -step):  # Explore q = 0, 1, -1, 2, -2, ...
            if step == 0 and q != 0:
                continue

            if q < q_min_range or q > q_max_range:
                continue

            p1 = p_n + 1 + 2 * k * q
            p2 = p_n - 1 + 2 * k * q

            if isprime(p1):
                q_min = q
                p_min = p1
                break
            if isprime(p2):
                q_min = q
                p_min = p2
                break
        if q_min is not None:
            break

    if q_min is not None:
        resultats.append((n, q_min, p_min, nombre_de_chiffres(p_min)))
        points.append((n, q_min))

# Affichage des résultats sous forme de grille
df = pd.DataFrame(resultats, columns=["n", "q", "p (nombre premier)", "Nombre de chiffres"])
print("\nRésultats sous forme de tableau :")
print(df)

# Sauvegarder les résultats dans un fichier CSV
df.to_csv("resultats_premiers.csv", index=False)
print("Les résultats ont été enregistrés dans le fichier : resultats_premiers.csv")

# Création de la figure pour afficher les points
fig, ax = plt.subplots(figsize=(12, 6))

# Ajout des points avec une couleur différente et taille fine
x_coords, y_coords = zip(*points) if points else ([], [])
ax.scatter(x_coords, y_coords, color="red", s=1, label="Points (n, q)")  # Change 's' to a smaller value for smaller points

# Configuration des axes
ax.set_xlabel("Valeurs de n (entiers)")
ax.set_ylabel("Valeurs de q (relatives)")
ax.set_title("Graphique : Correspondance entre n et q avec les plus petits q absolus")
ax.axhline(0, color="black", linewidth=0.8, linestyle="--")  # Ligne horizontale à y=0 pour référence

# Ajout de la légende
ax.legend()

# Afficher la grille
plt.grid(visible=True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

# Enregistrement de l'image
output_path = "graphique.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Graphique enregistré sous : {output_path}")

plt.show()
