# -*- coding: utf-8 -*-
"""
This module is designed to run PELE wth a single ligand to then cluster
positions and obtain energies.
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import os
import pathlib 
import argparse
import shutil
from distutils.dir_util import copy_tree
import numpy as np
from collections import Counter

def parse_args(args):
    """
    It parses the command-line arguments.
    Parameters
    ----------
    args : list[str]
        List of command-line arguments to parse
    Returns
    -------
    parsed_args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, dest = "input_folder",\
        default = 'LIG_Pele', help="Name of the directory where the simulation\
        is located.")
    parser.add_argument("-r", "--residue_name", type=str, dest = "residue_name",\
        default = 'LIG', help="Ligand's residue name.")
    parser.add_argument("-pdb", "--pdb_name", type=str, dest = "pdb_name",\
        default = None, help="Original pdb name.")

    parsed_args = parser.parse_args(args)

    return parsed_args

def linen_prepare(input_folder,
                  residue_name,
                  pdb_name):

    def path_definer(input_folder,
                     residue_name):

        """
        Defines all the paths that are going to be used

        Parameters
        ----------
        input_folder : str
            The path to the directory created by the induced fit simulation.
        clusters_folder : str
            Name of the directory where the directory clusters is located (results/analysis).

        Returns
        -------
        path_previous_simulation: str
            The path to the directory generated by the simulation we want to analyze the clusters
            from.
        path_clusters : str
            The path to the directory containing the pdbs of the representative poses for each
            cluster.
        path_energies_input : str 
            The path to the generated directory containing the input proportioned.
        path_energies_simulation : str
            The path to the generated directory containing all the necessary files to perform the
            PELE energy calculation.
        """        
        
        path = str(pathlib.Path().absolute())
        path_previous_simulation = path + '/' + input_folder

        if os.path.isdir(path_previous_simulation) == False:
            raise Exception('PathError: There is no folder with this name: '\
            + path_previous_simulation + '. Please check the path and the folder name.')

        path_energies = path + '/' + residue_name + '_linen_sim'
        path_energies_simulation = path_energies + '/simulation'

        if  os.path.exists(path_energies) == False:
            os.mkdir(path_energies)

        if  os.path.exists(path_energies_simulation) == False:
            os.mkdir(path_energies_simulation)

        return  path,path_previous_simulation, path_energies_simulation

    #---

    path, path_previous_simulation, path_energies_simulation = \
    path_definer(input_folder,residue_name)

    with open (os.path.join(path,pdb_name)) as filein:

        lines = (l for l in filein if residue_name in l)
        new_path = os.path.join(path_energies_simulation,'ligand.pdb')
        path_DataLocal = path_energies_simulation + '/DataLocal'

        if  os.path.exists(path_DataLocal) == False:
            os.mkdir(path_DataLocal)

        copy_tree(os.path.join(path_previous_simulation,'DataLocal'), path_DataLocal)
        
        with open(new_path, 'w') as fileout:
            
            fileout.writelines(lines)

def main(args):
    """
    It reads the command-line arguments and runs linen_results.
    Parameters
    ----------
    args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    Examples
    --------
    """

    linen_prepare(input_folder = args.input_folder,
                  residue_name = args.residue_name,
                  pdb_name = args.pdb_name)

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
       