from rdkit import Chem

def generate_mol_files(input_sdf, output_directory):
    suppl = Chem.SDMolSupplier(input_sdf)
    for mol in suppl:
        if mol is not None:
            
            # Retrieve ChEBI ID
            chebi_id = mol.GetProp('ChEBI ID')
            
            # Output MOL file
            output_file = f"{output_directory}/{chebi_id}.mol"
            
            # Write molecule to MOL file
            Chem.MolToMolFile(mol, output_file)

# Replace 'input.sdf' with the path to your input .sdf file
input_sdf = 'ChEBI_lite_3star.sdf'

# Replace 'output_directory' with the directory where you want to save the .mol files
output_directory = 'data/CHEBI/mol_files'

generate_mol_files(input_sdf, output_directory)