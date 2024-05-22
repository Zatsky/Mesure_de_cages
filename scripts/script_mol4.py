from openbabel import pybel

mol = next(pybel.readfile("sdf", "ChEBI_lite_3star.sdf"))
output_file = "data/op/molfile.mol"
with open(output_file, 'w') as f:
    f.write(mol.write("mol"))
