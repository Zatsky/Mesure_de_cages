import os
import shutil

# Répertoire source et destination
repertoire_source = 'data/smi_files_reduit'
repertoire_destination = 'data/chebi_smi'

# Parcourir le répertoire source
for fichier in os.listdir(repertoire_source):
    # Vérifier si le fichier a l'extension .mol
    if fichier.endswith('.mol'):
        # Extraire le numéro de fichier
        numero = fichier.split('.')[0]
        
        # Créer le chemin complet du fichier .mol et .smi
        chemin_mol = os.path.join(repertoire_source, fichier)
        chemin_smi = os.path.join(repertoire_source, numero + '.smi')
        
        # Vérifier si le fichier .smi existe et le déplacer
        if os.path.exists(chemin_smi):
            # Déplacer les fichiers .mol et .smi vers le répertoire destination
            shutil.move(chemin_mol, repertoire_destination)
            shutil.move(chemin_smi, repertoire_destination)
            print(f"Fichiers {numero}.mol et {numero}.smi déplacés avec succès.")
        else:
            # Déplacer seulement le fichier .mol
            shutil.move(chemin_mol, repertoire_destination)
            print(f"Fichier {numero}.mol déplacé avec succès.")