from openbabel import pybel

def extract_mol_files_from_sdf(sdf_file, output_directory):
    # Charger le fichier SDF
    molecules = pybel.readfile("sdf", sdf_file)
    
    # Parcourir les molécules et les écrire dans des fichiers .mol
    for i, mol in enumerate(molecules):
        output_file = f"{output_directory}{i}.mol"
        with open(output_file, 'w') as f:
            f.write(mol.write("mol"))

sdf_file = 'ChEBI_lite_3star.sdf'

# Replace 'output_directory' with the directory where you want to save the .mol files
output_directory = 'data/CHEBI/mol_files/'
# Appel de la fonction pour extraire les fichiers .mol
extract_mol_files_from_sdf(sdf_file, output_directory)