# -*- coding: utf-8 -*-
"""
This module is designed to run create N conformers of a ligand and prepare
a PELE energy calculation of them.
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import os
import pathlib 
import argparse
import shutil
import time
from distutils.dir_util import copy_tree
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
from rdkit.ML.Cluster import Butina

def parse_args(args):
    """
    Function
    ----------
    It parses the command-line arguments.

    Parameters
    ----------
    - args : list[str]
        List of command-line arguments to parse

    Returns
    ----------
    - parsed_args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-r", "--residue_name", type=str, dest = "residue_name",\
        default = 'LIG', help="Ligand's chain name.")
    parser.add_argument("-f", "--file", type=str, dest = "input_file",\
        default = None, help="Name of the file to analyze.")
    parser.add_argument("-cn", "--conformations_number", type=int, dest = "conformations_number",\
        default = 50, help="Number of conformers to be generated for the inputed file.")
    parser.add_argument("-co", "--conf_file_name", type=str, dest = "conf_file_name",\
        default = 'pele', help="Name of the .conf file used for the simulation.")

    parsed_args = parser.parse_args(args)

    return parsed_args

def path_definer(input_file,residue_name):
    """
    Function
    ----------
    Defines all the paths that are going to be used
    Parameters

    Parameters
    ----------
    - input_folder : str
        The path to the directory created by the induced fit simulation.

    Returns
    ----------
    - path_output: str
        The path to the output directory generated in the simulation.
    - path_reports : str
        The path to the generated directory that will contain the copied reports.
        cluster.
    """        
    
    path = str(pathlib.Path().absolute())
    path_pdb = os.path.join(path,input_file)

    path_energies = os.path.join(path,residue_name + '_linen_con')
    path_energies_input = os.path.join(path_energies,'input')
    path_energies_simulation = os.path.join(path_energies,'simulation')
    path_energies_clusters = os.path.join(path_energies_simulation,'clusters')
    path_energies_DataLocal = os.path.join(path_energies,'DataLocal')

    if  os.path.exists(path_energies) == False:
        os.mkdir(path_energies)

    if  os.path.exists(path_energies_input) == False:
        os.mkdir(path_energies_input)

    if  os.path.exists(path_energies_simulation) == False:
        os.mkdir(path_energies_simulation)
    
    if  os.path.exists(path_energies_clusters) == False:
        os.mkdir(path_energies_clusters)

    if  os.path.exists(path_energies_DataLocal) == False:
        os.mkdir(path_energies_DataLocal)

    return path, path_pdb, path_energies_input, path_energies, path_energies_clusters, path_energies_DataLocal

def linen_conformer(conformations_number,
                    conf_file_name,
                    path,
                    path_pdb,
                    path_energies_input,
                    path_energies,
                    path_energies_clusters,
                    path_energies_DataLocal):
    
    print(' ')
    print('*******************************************************************')
    print('*                    peleLigandConformations                      *')
    print('* --------------------------------------------------------------- *')
    print('*      Ligand\'s internal energy from conformer generation         *')
    print('*******************************************************************')
    print(' ')

    start_time = time.time()
    
    shutil.copy(path_pdb, path_energies_input)

    if os.path.isdir(os.path.join(path,'DataLocal')):
        copy_tree(os.path.join(path,'DataLocal'), path_energies_DataLocal)
    
    else: 
        print(
        '                              WARNING:                               \n'
        '     No DataLocal folder was found in this directory. If linen_a.py  \n'
        '     is to be used, the DataLocal folder will have to be copied in   \n'
        '     the /LIG_linen_con folder.                                      \n'
        '\n'
        )
    
    if os.path.isfile(os.path.join(path,conf_file_name + '.conf')):
        shutil.copy(os.path.join(path,conf_file_name + '.conf'), path_energies)

    else: 
        print(
        '                              WARNING:                               \n'
        '     No pele.conf file was found in this directory. If linen_a.py    \n'
        '     is to be used, the script will assign default values to the     \n'
        '     forcefield and the solvent.                                     \n'
        '\n'
        )

    if 'CONECT' not in open(path_pdb).read():
        raise Exception('NoConnectivityError: The pdb of the ligand must have connectivity. Otherwise it cannot be assigned and the conformations will make no sense.')

    # Preparing molecule
    m = Chem.MolFromPDBFile(path_pdb) 
    weight = ExactMolWt(m)
    rotatable_bonds = CalcNumRotatableBonds(m)

    m = AllChem.RemoveHs(m)
    mh = AllChem.AddHs(m, addCoords=True)

    # Print
    print(' -   Input:')
    print('     -   File to analyze:',path_pdb.split('/')[-1] + '.')
    print('     -   Number of conformations:',str(conformations_number) + '.')
    print(' -   Information:')
    print('     -   Molecular weight:', weight)
    print('     -   Number of rotatable bonds:',rotatable_bonds)
    print(' -   Generating conformations...')
    # 

    # Creating conformations and calculating
    cids = AllChem.EmbedMultipleConfs(mh, numConfs=conformations_number, numThreads=0)

    conformer_properties_dictionary = {}

    for cid in cids:

        ff = AllChem.MMFFGetMoleculeForceField(mh, AllChem.MMFFGetMoleculeProperties(mh), confId=cid)
        ff.Initialize()
        ff.CalcEnergy()
        results = {}
        results["converged"] = ff.Minimize(maxIts=2000)
        results["energy_abs"] = ff.CalcEnergy()

        conformer_properties_dictionary[cid] = results

    dmat = AllChem.GetConformerRMSMatrix(mh, prealigned=False)

    rms_clusters = Butina.ClusterData(dmat, mh.GetNumConformers(), 2.0, isDistData=True, reordering=True)

    #
    print(' -   Number of clusters:', len(rms_clusters))
    #

    cluster_number = 0
    min_energy = 10e10

    for cluster in rms_clusters:

        cluster_number += 1
        rmslist = []

        AllChem.AlignMolConformers(mh, confIds=cids, RMSlist=rmslist)

        for cid in cluster:

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

    labels = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N',\
    'O','P','Q','R','S','T','U','V','W','X','Y','Z']

    cont = -1

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

            Chem.rdmolfiles.MolToPDBFile(mh, path_energies_clusters + '/cluster_' + labels[cont] + '.pdb', confId = cid)

    #
    print(' ')
    print('                    --Duration of the execution--                   ')
    print('                      %s seconds' % (time.time() - start_time)  )
    print(' ')
    print('*******************************************************************')
    #

def main(args):
    """
    Function
    ----------
    It runs statistics function. 

    Parameters
    ----------
    - args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    path, path_pdb, path_energies_input, path_energies, path_energies_clusters, path_energies_DataLocal =\
    path_definer(input_file=args.input_file,residue_name=args.residue_name)

    linen_conformer(conformations_number = args.conformations_number,
                    conf_file_name = args.conf_file_name,
                    path=path,
                    path_pdb = path_pdb,
                    path_energies_input = path_energies_input,
                    path_energies = path_energies,
                    path_energies_clusters = path_energies_clusters,
                    path_energies_DataLocal= path_energies_DataLocal
                    )

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)