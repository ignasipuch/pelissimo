# -*- coding: utf-8 -*-

# Imports
import sys
import os
import pathlib 
import argparse
import shutil
import matplotlib.pyplot as plt 
import time
from itertools import cycle

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
from rdkit.ML.Cluster import Butina
import numpy as np



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

# Print
print(' -   Input:')
print('     -   File to analyze:',input_file + '.')
print('     -   Number of conformations:',str(conformations_number) +'.')
#

# Original conformation information
m = Chem.MolFromPDBFile(path_ligand) 
print(Chem.MolToSmiles(m))
weight = ExactMolWt(m)
rotatable_bonds = CalcNumRotatableBonds(m)

# Making sure all the Hydrogen's are right.
m = AllChem.RemoveHs(m)
mh = AllChem.AddHs(m, addCoords=True)

# Print 
print(' -   Information:')
print('     -   Molecular weight:', weight)
print('     -   Number of rotatable bonds:',rotatable_bonds)
print(' -   Generating conformations...')
# 

# Creating conformations
cids = AllChem.EmbedMultipleConfs(mh, numConfs=conformations_number, numThreads=0)

# Dictionary to store information
conformer_properties_dictionary = {}

# For all the conformations created calculate the energy of minimized conformation.
# Each entry of the dictionary will have two subentries with the convergence and the energy value.
for cid in cids:
    
    ff = AllChem.MMFFGetMoleculeForceField(mh, AllChem.MMFFGetMoleculeProperties(mh), confId=cid)
    ff.Initialize()
    ff.CalcEnergy()
    results = {}
    results["converged"] = ff.Minimize(maxIts=2000)
    results["energy_abs"] = ff.CalcEnergy()

    conformer_properties_dictionary[cid] = results

# Calculation of the RMSD of all the conformations to the original molecule.
dmat = AllChem.GetConformerRMSMatrix(mh, prealigned=False)

# RDKit Machine learning module to cluster the conformations and create clusters.
# The result is a tuple of tuples where each tuple is a cluster and the content are the 
# the conformation numbers.
rms_clusters = Butina.ClusterData(dmat, mh.GetNumConformers(), 2.0, isDistData=True, reordering=True)

print(' -   Number of clusters:',len(rms_clusters))

cluster_number = 0
min_energy = 10e10

# For all the entry in the tuple
for cluster in rms_clusters:
    
    cluster_number += 1
    rmslist = []

    # Aligning
    AllChem.AlignMolConformers(mh, confIds=cids, RMSlist=rmslist)
    
    for cid in cluster:

        # Check minimium energies
        e = results["energy_abs"]

        if e < min_energy:

            min_energy = e
            results = conformer_properties_dictionary[cid]
            results["cluster_number"] = cluster_number
            results["cluster_centroid"] = cluster[0] + 1
            index = cluster.index(cid)

            if index > 0:

            	results["rms_to_centroid"] = rmslist[index-1]

            else:

            	results["rms_to_centroid"] = 0.0

#w = Chem.SDWriter(path_energies_simulation + '/clusters.pdb')

cont = -1
labels = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']

for cluster in rms_clusters:

    cont += 1

    for cid in cluster:

        for name in mh.GetPropNames():

            mh.ClearProp(name)

        conformer_properties = conformer_properties_dictionary[cid]
        mh.SetIntProp("conformer_id", cid + 1)

        for key in conformer_properties.keys():

            mh.SetProp(key, str(conformer_properties[key]))
            e = conformer_properties["energy_abs"]

            if e:

                mh.SetDoubleProp("energy_delta", e - min_energy)

        Chem.rdmolfiles.MolToPDBFile(mh, path_energies_simulation + '/cluster_' + labels[cont] + '.pdb', confId = cid)

#
print(' ')
print('                    --Duration of the execution--                   ')
print('                      %s seconds' % (time.time() - start_time)  )
print(' ')
print('*******************************************************************')
#