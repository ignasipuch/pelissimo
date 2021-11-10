import argparse 
import pathlib 
import matplotlib.pyplot as plt 
import time

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
from sklearn.cluster import MeanShift, estimate_bandwidth
import numpy as np
from itertools import cycle


# ----------------------------------------------------------------------- #
#                 READING AND STORING INFORMATION FROM FILES              #
# ----------------------------------------------------------------------- #

start_time = time.time()

# Parse command line
parser = argparse.ArgumentParser()

# Inputs
parser.add_argument("-f", "--file", type=str, dest = "input_file",\
    default = None, help="Name of the file to analyze.")
parser.add_argument("-cn", "--conformations_number", type=int, dest = "conformations_number",\
    default = 10, help="Number of conformers to be generated for the inputed file.")

args = parser.parse_args()

input_file = args.input_file
conformations_number = args.conformations_number

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
print('- Creating ' + str(conformations_number) + ' conformations.')

# Original conformation
m = Chem.MolFromPDBFile(filepath) 
weight = ExactMolWt(m)
rotatable_bonds = CalcNumRotatableBonds(m)

mh = AllChem.AddHs(m, addCoords=True)
mp = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant='MMFF94s')
ff = AllChem.MMFFGetMoleculeForceField(mh, mp)
e_bound = ff.CalcEnergy()

# Print
print('- Molecular weight:', weight)
print('- Energy of the bound conformation:',e_bound)
print('- Number of rotatable bonds:',rotatable_bonds)

# Creating conformations
cids = AllChem.EmbedMultipleConfs(mh, numConfs=conformations_number, numThreads=0) 
mp = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant='MMFF94s')

# RMSD calculation
rmslist = []
AllChem.AlignMolConformers(mh, RMSlist=rmslist)                           
rms = AllChem.GetConformerRMS(mh, 1, conformations_number-1, prealigned=True) 

# Optimization
res = AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=20000, numThreads=0) 

# Collecting energies
energies_list = [e[1] for e in res]
energies_list = [e_bound - e for e in energies_list ]
del energies_list[0]

# Clustering results
X = np.column_stack((rmslist,energies_list))

bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=len(energies_list))
ms = MeanShift().fit(X)
labels = ms.labels_
cluster_centers = ms.cluster_centers_
n_clusters_ = len(np.unique(labels))
cluster_centers = ms.cluster_centers_

# Print ---
print(' ')
print(' ----------------------------------------------')
print('|                   RESULTS                    |')
print(' ----------------------------------------------')
print(' ')
print('- Biggest energy difference: ', max(energies_list))
print('- Plotting figure')
print('- Number of clusters:',n_clusters_)


# Plots ---------------
colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(X[my_members, 0], X[my_members, 1], col + '.')

plt.xlabel('RMSD ($\mathring{A}$)')
plt.ylabel('$\Delta$Energy kcal/mol ')
plt.ylim(0,max(energies_list) + max(energies_list)/10)
plt.title('Conformers ' + input_file.strip('_lig.pdb'))
plt.savefig('conformer_' + input_file.strip('_lig.pdb') + '.png', format='png')

# ---------------------
#plt.scatter(rmslist, energies_list)
#plt.title('Conformers')
#plt.ylim(0,max(energies_list) + max(energies_list)/10)
#plt.xlabel('RMSD ($\mathring{A}$)')
#plt.ylabel('$\Delta$Energy kcal/mol ')
#plt.savefig('conformer_' + input_file.strip('.pdb') + '.png', format='png')

print("- Duration of the execution: %s seconds" % (time.time() - start_time))
print(' ')

print(' ----------------------------------------------')