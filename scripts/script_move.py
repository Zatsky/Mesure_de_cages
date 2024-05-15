import os
import shutil

# Répertoire source et destination
repertoire_source = 'data/DEFAULT/smi_files'
repertoire_destination = 'data/DEFAULT/mol_files'

# Parcourir le répertoire source
for fichier in os.listdir(repertoire_source):
    # Vérifier si le fichier a l'extension .mol
    if fichier.endswith('.mol'):
        chemin_mol = os.path.join(repertoire_source, fichier)
        #os.remove(chemin_mol)
        
        shutil.move(chemin_mol, repertoire_destination)