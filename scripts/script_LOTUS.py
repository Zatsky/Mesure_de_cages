from rdkit import Chem
from rdkit.Chem import PandasTools
import urllib.request
import requests
import zipfile
import shutil
import os
import csv

#--- Fonctions ---#

def sdf_to_csv(sdf_file, output_file, delimiter=','):
    """
    Convert an SDF file to a CSV or TSV file.
    
    Args:
        sdf_file (str): Path to the input SDF file.
        output_file (str): Path to the output CSV or TSV file.
        delimiter (str, optional): Delimiter to use in the output file. Default is ',' (CSV).
    """
    # Read the SDF file
    suppl = Chem.SDMolSupplier(sdf_file)

    # Convert to a pandas DataFrame
    df = PandasTools.LoadSDF(sdf_file, smilesName='SMILES', molColName='Molecule', includeFingerprints=False)

    # Drop the 'Molecule' column if present
    if 'Molecule' in df.columns:
        df = df.drop('Molecule', axis=1)

    # Save the DataFrame as CSV or TSV
    df.to_csv(output_file, sep=delimiter, index=False)
    
    print(f"Conversion complete! Output file: {output_file}")
#--- Script ---#

# Vérifier si la bibliothèque RDKit est disponible
try:
    from rdkit import Chem

except ImportError :
    print("Téléchargement de la bibliothèque RDKit...")
    os.system("pip install rdkit")
    from rdkit import Chem


# URL du fichier ZIP à télécharger
url = 'https://lotus.naturalproducts.net/download/sdf'
downloaded_file_name = 'LOTUS_DB_LATEST.zip'

# Télécharger le fichier ZIP
urllib.request.urlretrieve(url, downloaded_file_name)

# Chemin vers le répertoire où vous voulez extraire le contenu du fichier ZIP
extracted_dir_path = os.getcwd()  # Utiliser le répertoire courant

# Extraction du contenu du fichier ZIP
with zipfile.ZipFile(downloaded_file_name, 'r') as zip_ref:
    zip_ref.extractall(extracted_dir_path)

unzipped_file_name = 'LOTUS_2021_03_simple.sdf'
os.remove(downloaded_file_name)

sdf_file = unzipped_file_name  # Replace with the actual path to your SDF file
csv_file = 'data/output.csv'

sdf_to_csv(sdf_file, csv_file)

with open("scripts/script_smi_to_mol.sh", "w") as f_script :
	with open(csv_file, 'r') as f :
		csvreader = csv.reader(f)
		i = 0
		for row in csvreader :
			if row[1] != "lotus_id":
				petit_smile = row[4]
				name = "_".join(row[0].split(" "))
				if name+".mol" not in os.listdir("data/lotus_smi/") :
					with open("data/lotus_smi/%s.smi"%name, 'w') as f :
						f.write("%s %s\n"%(petit_smile, name))
					i += 1
				
					f_script.write("obgen data/lotus_smi/%s.smi > data/lotus_smi/%s.mol\n"%(name, name))