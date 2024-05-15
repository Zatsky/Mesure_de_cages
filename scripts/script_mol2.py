from rdkit import Chem
import os

def extract_mol_blocks_from_sdf(sdf_file):
    mol_blocks = []
    with open(sdf_file, 'r') as f:
        current_mol_block = ''
        for line in f:
            if line.strip() == '$$$$':
                mol_blocks.append(current_mol_block)
                current_mol_block = ''
            else:
                current_mol_block += line
    return mol_blocks

def write_mol_files(mol_blocks,output):
    for i, mol_block in enumerate(mol_blocks):
        mol = Chem.MolFromMolBlock(mol_block)
        if has_connectivity(mol):
            # Vous pouvez ajouter d'autres conditions de filtrage ici si nécessaire
            if mol.GetNumAtoms() > 0:
                with open(output + f'molecule_{i}.mol', 'w') as f:
                    f.write(mol_block)

def has_connectivity(mol):
    if mol:
        bonds = mol.GetBonds()
        if bonds:
            return True
    return False

def remove_files_with_no_connectivity(sdf_file,output):
    with open(sdf_file, 'r') as f:
        current_mol_block = ''
        i = 0
        for line in f:
            if line.strip() == '$$$$':
                i += 1
                if not has_connectivity(current_mol_block):
                    # Si le bloc de molécule actuel n'a pas de connectivité, supprimer le fichier
                    file_name = sdf_file.split('.')[0] + f'_molecule_{i}.sdf'  # Nom du fichier SDF à supprimer
                    os.remove(file_name)
                current_mol_block = ''
            else:
                current_mol_block += line

# Replace 'input.sdf' with the path to your input .sdf file
sdf_file = 'ChEBI_lite_3star.sdf'

# Replace 'output_directory' with the directory where you want to save the .mol files
output_directory = 'data/CHEBI/mol_files/'
mol_blocks = extract_mol_blocks_from_sdf(sdf_file)
write_mol_files(mol_blocks, output_directory)
#remove_files_with_no_connectivity(sdf_file,output_directory)