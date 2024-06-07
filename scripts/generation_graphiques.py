import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
from collections import defaultdict
import numpy as np

# Fonction pour lire les données du fichier CSV
def lire_donnees_csv(nom_fichier):
    mesures = []
    with open(nom_fichier, 'r') as fichier:
        lecteur_csv = csv.reader(fichier)
        for i, ligne in enumerate(lecteur_csv):
            try:
                mesure = float(ligne[1])
                mesures.append(mesure)
            except ValueError:
                pass
    return mesures

# Lire les données du fichier CSV
def gen_graphique(arg1):
    if arg1 == "CHEBI":
        max_x1 = 1000
        max_x2 = 250
        max_y1 = 50
        max_y2 = 800
    elif arg1 == "LOTUS":
        max_x1 = 100
        max_x2 = 250
        max_y1 = 1000
        max_y2 = 13000
    elif arg1 == "CHIMISTE":
        max_x1 = 50
        max_x2 = 130
        max_y1 = 30
        max_y2 = 30
    else:
        max_x1 = 100
        max_x2 = 250
        max_y1 = 50
        max_y2 = 800

    result_dir = "data/" + arg1 + "/results/"
    mesures = lire_donnees_csv(result_dir + "liste_mesure_alpha_connexe.csv")

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

    plt.plot(mesures_triees, nombres_molecules_cumulatifs, drawstyle='steps-post')
    plt.grid(True)

    plt.xlim(0, max_x1)
    plt.ylim(0, max_y1)  # Fixer la limite supérieure de l'axe y

    plt.title("Nombre de molécules en fonction de la cagitude")
    plt.xlabel("Cagitude")
    plt.ylabel("Nombre de molécules")

    plt.savefig(result_dir + "graphique_ligne.png")
    plt.close()

    # Créer le graphique en barres
    bins = np.arange(0, 1200 + 5, 5)

    # Utilisation de np.histogram pour regrouper les données par intervalles de 5
    hist, bin_edges = np.histogram(mesures, bins=bins)

    # Création du graphique à barres
    plt.bar(bin_edges[:-1], hist, width=5, align='edge')
    plt.grid(True)

    plt.xlim(0, max_x1)
    plt.ylim(0, max_y1)  # Fixer la limite supérieure de l'axe y

    # Ajouter des flèches pour indiquer que certaines valeurs dépassent la limite supérieure
    for i, count in enumerate(hist):
        if count > max_y1:
            plt.annotate(
                '⬆',
                (bin_edges[i], max_y1),
                ha='center',
                va='bottom',
                fontsize=12,
                color='red'
            )

    plt.title("Nombre de molécules en fonction de la cagitude")
    plt.xlabel("Cagitude")
    plt.ylabel("Nombre de molécules")

    plt.savefig(result_dir + "histogramme_discret.png")
    plt.close()

    # Créer le graphique de distribution cumulative
    plt.plot(mesures_triees, nombres_molecules_cumulatifs)
    plt.grid(True)

    plt.xlim(-10, max_x2)
    plt.ylim(-10, max_y2)  # Fixer la limite supérieure de l'axe y

    plt.title("Nombre de molécules inférieur ou égal à une cagitude")
    plt.xlabel("Cagitude")
    plt.ylabel("Nombre de molécules")

    plt.savefig(result_dir + "graphique_distribution_cumulative.png")
    plt.close()
gen_graphique("CHEBI")