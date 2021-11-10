# -*- coding: utf-8 -*-

# Imports
import sys
import os
import pathlib 
import argparse
import shutil
import matplotlib.pyplot as plt 
import time

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
from sklearn.cluster import MeanShift, estimate_bandwidth
import numpy as np
from itertools import cycle


# *********************************************************************** #
#                        PREPARING READING OF FILES                       #
# *********************************************************************** #

# Time
start_time = time.time()

# Parse command line
parser = argparse.ArgumentParser()

# Inputs
parser.add_argument("-r", "--residue_name", type=str, dest = "residue_name",\
    default = 'LIG', help="Ligand's chain name.")
parser.add_argument("-f", "--file", type=str, dest = "input_file",\
    default = None, help="Name of the file to analyze.")
parser.add_argument("-cn", "--conformations_number", type=int, dest = "conformations_number",\
    default = 50, help="Number of conformers to be generated for the inputed file.")

args = parser.parse_args()

residue_name = args.residue_name
input_file = args.input_file
conformations_number = args.conformations_number

# Creating a variable with the path of the file
path = str(pathlib.Path().absolute())
filepath = path + '/' + input_file

path_energies = path + '/' + residue_name + '_peleConf'
path_energies_input = path_energies + '/' + 'input'
path_energies_simulation = path_energies + '/' + 'simulation'

if  os.path.exists(path_energies) == False:
    os.mkdir(path_energies)

if  os.path.exists(path_energies_input) == False:
    os.mkdir(path_energies_input)

if  os.path.exists(path_energies_simulation) == False:
    os.mkdir(path_energies_simulation)

# *********************************************************************** #
#                    COPYING FILES TO NEW DIRECTORY                       #
# *********************************************************************** #

# Printing header
print(' ')
print('*******************************************************************')
print('*                    peleLigandConformations                      *')
print('* --------------------------------------------------------------- *')
print('*      Ligand\'s internal energy from conformer generation         *')
print('*******************************************************************')
print(' ')

# Copying complexed crystal
shutil.copy(filepath, path_energies_input)

# Ligand's path
path_ligand = os.path.join(path_energies_simulation,'ligand.pdb')

# Removing the protein from the pdb.
with open(filepath) as filein:
    
    lines = (l for l in filein if residue_name in l)
        
    with open(path_ligand, 'w') as fileout:

        fileout.writelines(lines)

#
print(' -   Input:')
print('     -   File to analyze:',input_file + '.')
print('     -   Number of conformations:',str(conformations_number) +'.')
#

# Original conformation information
m = Chem.MolFromPDBFile(path_ligand) 
weight = ExactMolWt(m)
rotatable_bonds = CalcNumRotatableBonds(m)

m = AllChem.RemoveHs(m)
mh = AllChem.AddHs(m, addCoords=True)
mp = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant='MMFF94s')
ff = AllChem.MMFFGetMoleculeForceField(mh, mp)
e_bound = ff.CalcEnergy()

# 
print(' -   Information:')
print('     -   Molecular weight:', weight)
print('     -   Number of rotatable bonds:',rotatable_bonds)
print(' -   Generating conformations...')
# 

# Creating conformations
cids = AllChem.EmbedMultipleConfs(mh, numConfs=conformations_number, numThreads=0) 

# RMSD calculation
rmslist = []
AllChem.AlignMolConformers(mh, RMSlist=rmslist)                           
rms = AllChem.GetConformerRMS(mh, 1, conformations_number-1, prealigned=True) 

# Optimization
res = AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=20000, numThreads=0) 

# Collecting energies

energies_list = [e[1] for e in res]
convergence = [e[0] for e in res]
energies_list = [e - min(energies_list)  for e in energies_list ]
del energies_list[0]

# Clustering results
X = np.column_stack((energies_list,rmslist))

bandwidth = estimate_bandwidth(X, quantile=0.3, n_samples=len(energies_list))
ms = MeanShift(bandwidth=bandwidth, bin_seeding=True).fit(X)
labels = ms.labels_
cluster_centers = ms.cluster_centers_
n_clusters_ = len(np.unique(labels))
cluster_centers = ms.cluster_centers_

# 
print(' -   Plotting figure')
print(' -   Number of clusters:',n_clusters_)
#

# ------- #
#  Plot 1 #
# ------- #
colors = cycle('grcmykbgrcmykbgrcmykbgrcmykb')

plt.figure()
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(X[my_members, 1], X[my_members, 0], col + '.')

plt.xlabel('RMSD ($\mathring{A}$)')
plt.ylabel('$\Delta$Energy kcal/mol ')
plt.ylim(- max(energies_list)/10,max(energies_list) + max(energies_list)/10)
plt.title('Conformers ' + input_file.strip('_lig.pdb'))
plt.plot(cluster_centers[:,1],cluster_centers[:,0],'Dk',markersize=6)
plt.savefig(path_energies + '/conformer_' + input_file.strip('_lig.pdb') + '_clusters.png', format='png')
plt.close()

# ------- #
#  Plot 2 #
# ------- #
plt.figure()
plt.scatter(rmslist, energies_list)
plt.title('Conformers')
plt.ylim(0- max(energies_list)/10,max(energies_list) + max(energies_list)/10)
plt.xlabel('RMSD ($\mathring{A}$)')
plt.ylabel('$\Delta$Energy kcal/mol ')
plt.savefig(path_energies + '/conformer_' + input_file.strip('.pdb') + '.png', format='png')
plt.close()

#
print(' ')
print('                    --Duration of the execution--                   ')
print('                      %s seconds' % (time.time() - start_time)  )
print(' ')
print('*******************************************************************')
#


for cluster in rms_clusters:

    cont += 1

    for cid in cluster:

        for name in mh.GetPropNames():

            mh.ClearProp(name)

        conformer_properties = conformer_properties_dictionary[cid]
        mh.setIntProp("conformer_id", cid + 1)

        for key in conformer_properties.keys():

            mh.setProp(key, str(conformer_properties[key]))
            e = conformer_properties["energy_abs"]

            if e:

                mh.SetDoubleProp("energy_delta", e - min_energy)

        Chem.rdmolfiles.MoltoPDBFile(mh, path_energies_simulation + '/cluster_' + labels[cont] + '.pdb', confId = cid)