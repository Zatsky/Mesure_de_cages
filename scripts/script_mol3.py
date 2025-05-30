import os
from openbabel import pybel
import sys

def extract_mol_files_from_sdf(sdf_file, output_directory,id):
    # Création du répertoire de sortie s'il n'existe pas
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Charger le fichier SDF
    molecules = pybel.readfile("sdf", sdf_file)
    # Parcourir les molécules et les écrire dans des fichiers .mol
    for mol in molecules:
        try:
            # Extraire l'ID ChEBI des propriétés de la molécule
            chebi_id = mol.data[id]
            
            # Remplacer les caractères non autorisés dans les noms de fichiers
            chebi_id = chebi_id.replace(":", "_")
            output_file = f"{output_directory}{chebi_id}.mol"
            with open(output_file, 'w') as f:
                f.write(mol.write("mol"))
        except KeyError:
            # Gérer les cas où l'ID ChEBI n'est pas trouvé
            print("ChEBI ID not found for one of the molecules.")

arg1 = sys.argv[1]
if arg1 == "CHEBI":
    sdf_file = 'ChEBI_lite_3star.sdf'
    id = "ChEBI ID"
elif arg1 == "LOTUS":
    sdf_file = 'LOTUS_2021_03_simple.sdf'
    id = "lotus_id"

# Remplacez 'output_directory' par le répertoire où vous souhaitez enregistrer les fichiers .mol
output_directory = 'data/'+arg1+'/mol_files/'

# Appel de la fonction pour extraire les fichiers .mol
extract_mol_files_from_sdf(sdf_file, output_directory, id)
