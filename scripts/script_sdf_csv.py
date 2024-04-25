from rdkit import Chem
from rdkit.Chem import PandasTools
import urllib.request
import gzip
import shutil
import os
import csv

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


url = 'https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_lite_3star.sdf.gz'
downloaded_file_name = 'ChEBI_lite_3star.sdf.gz'

urllib.request.urlretrieve(url, downloaded_file_name)
unzipped_file_name = 'ChEBI_lite_3star.sdf'
# Example usage

with gzip.open(downloaded_file_name, 'rb') as f_in, open(unzipped_file_name, 'wb') as f_out :
    shutil.copyfileobj(f_in, f_out)

os.remove(downloaded_file_name)
sdf_file = unzipped_file_name  # Replace with the actual path to your SDF file
csv_file = 'data/output.csv'

sdf_to_csv(sdf_file, csv_file)

with open("scripts/script_smi_to_mol.sh", "w") as f_script :
	with open(csv_file, 'r') as f :
		csvreader = csv.reader(f)
		i = 0
		for row in csvreader :
			if row[1] != "ChEBI Name":
				petit_smile = row[4]
				name = "_".join(row[0].split(" "))
				if name+".mol" not in os.listdir("data/smi_files_reduit/") :
					with open("data/smi_files_reduit/%s.smi"%name, 'w') as f :
						f.write("%s %s\n"%(petit_smile, name))
					i += 1
				
					f_script.write("obgen data/smi_files_reduit/%s.smi > data/smi_files_reduit/%s.mol\n"%(name, name))