import argparse 
import pathlib 
import matplotlib.pyplot as plt 

from rdkit import Chem
from rdkit.Chem import AllChem

# ----------------------------------------------------------------------- #
#                 READING AND STORING INFORMATION FROM FILES              #
# ----------------------------------------------------------------------- #

# Parse command line
parser = argparse.ArgumentParser()

# Inputs
parser.add_argument("-f", "--file", type=str, dest = "input_file",\
    default = None, help="Name of the file to analyze.")

args = parser.parse_args()

input_file = args.input_file

# Creating a variable with the path of the file
path = str(pathlib.Path().absolute())
filepath = path + '/' + input_file

# Print ---
print(' ')
print(' ----------------------------------------------')
print('|                READING INPUT                 |')
print(' ----------------------------------------------')
print(' ')
print('- File to analyze:',filepath)

m = Chem.MolFromPDBFile(filepath) 

mh = AllChem.AddHs(m, addCoords=True)
mp = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant='MMFF94s')
ff = AllChem.MMFFGetMoleculeForceField(mh, mp)
e_bound = ff.CalcEnergy()

# Print
print('- Energy of the bound conformation:',e_bound)

for atom in mh.GetAtoms():
    if not atom.GetAtomicNum() == 1:

        idx = atom.GetIdx()
        ff.MMFFAddPositionConstraint(idx, maxDispl=0.5, forceConstant=100)

ff.Minimize(maxIts=10000)
e_min = ff.CalcEnergy()


print('- Energy of the minimized conformation:',e_min)
print('- Difference:', e_min-e_bound)
print(' ')
