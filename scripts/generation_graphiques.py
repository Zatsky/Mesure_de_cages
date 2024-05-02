import csv
import matplotlib.pyplot as plt
from collections import defaultdict
import sys

# Fonction pour lire les données du fichier CSV
def lire_donnees_csv(nom_fichier):
    mesures = []
    with open(nom_fichier, 'r') as fichier:
        lecteur_csv = csv.reader(fichier)
        for i,ligne in enumerate(lecteur_csv):
            try:
                mesure = float(ligne[1])
                mesures.append(mesure)
            except ValueError:
                pass
    return mesures

arg1 = sys.argv[1]
# Lire les données du fichier CSV
result_dir = "data/" + arg1 + "/results/"
mesures = lire_donnees_csv(result_dir + "liste_mesure_alpha.csv")

# Compter le nombre de molécules pour chaque mesure
nombre_molecules_par_mesure = defaultdict(int)
for mesure in mesures:
    nombre_molecules_par_mesure[mesure] += 1

# Trier les mesures et le nombre de molécules par mesure
mesures_triees = sorted(nombre_molecules_par_mesure.keys())
nombres_molecules_cumulatifs = [0] * len(mesures_triees)

nombres_molecules = [nombre_molecules_par_mesure[mesure] for mesure in mesures_triees]

nombre_molecules_cumulatif = 0
for i, mesure in enumerate(mesures_triees):
    nombre_molecules_cumulatif += nombre_molecules_par_mesure[mesure]
    nombres_molecules_cumulatifs[i] = nombre_molecules_cumulatif
# Créer le graphique en barres
plt.bar(mesures_triees, nombres_molecules)

# Ajouter des titres et des étiquettes
plt.title("Nombre de molécules pour chaque mesure")
plt.xlabel("Cagitude")
plt.ylabel("Nombre de molécules")

plt.savefig(result_dir+"graph1.png")
# Afficher le graphique

plt.close()
plt.plot(mesures_triees, nombres_molecules_cumulatifs)

# Ajouter des titres et des étiquettes
plt.title("Nombre de molécules inférieur ou égal à une mesure")

plt.xlabel("Cagitude")
plt.ylabel("Nombre de molécules")
plt.savefig(result_dir+"graph2.png")