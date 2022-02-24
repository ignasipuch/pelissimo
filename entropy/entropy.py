# -*- coding: utf-8 -*-
"""
This script ensambles other scripts and calculates entropy
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import dihedral_clustering as dc
import lice 

import sys
import os
import pathlib
import argparse

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

    parser.add_argument("-f", "--file", type=str, dest="input_file",
                        default=None, required=True, help="Name of the file corresponding to the isolated ligand with connectivity.")
    parser.add_argument("-d", "--directory", type=str, dest="input_folder",
                        default='LIG_Pele', help="Name of the directory where the simulation\
        is located.")
    parser.add_argument("-o", "--output_directory", type=str, dest="output_folder",
                        default='output', help="Name of the output directory where the simulation\
        is located.")
    parser.add_argument("-r", "--residue_name", type=str, dest="residue_name",
                        default='LIG', help="Ligand's residue name.")
    parser.add_argument("-cm", "--clustering_method", type=str, dest="clustering_method",
                        default='bin', help="Method to cluster data: bin or kmeans.")
    parser.add_argument("-nc", "--n_clusters", type=int, dest="n_clusters",
                        default=0, help="Number of clusters to cluster the data.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def ensambler(input_folder,
              residue_name,
              input_file,
              output_folder,
              clustering_method,
              n_clusters):

    path = str(pathlib.Path().absolute())
    path_pl_simulation = os.path.join(path,input_folder)
    path_l_simulation = os.path.join(path,residue_name + '_linen_cry')

    os.chdir(path_pl_simulation)
    os.system('python /home/bsc72/bsc72825/projects/code/dihedral_clustering.py -f ' + input_file + ' -d ' + output_folder  + ' -r ' + residue_name + ' -cm ' + clustering_method + ' -nc ' + str(n_clusters))
    os.chdir(path_l_simulation)
    os.system('python /home/bsc72/bsc72825/projects/code/dihedral_clustering.py -f ' + input_file + ' -d ' + output_folder  + ' -r ' + residue_name + ' -cm ' + clustering_method + ' -nc ' + str(n_clusters))
    os.chdir(path)
    os.system('python /home/bsc72/bsc72825/projects/code/lice.py -d ' + input_folder  + ' -r ' + residue_name)
    

def main(args):

    ensambler(input_folder=args.input_folder,
              residue_name=args.residue_name,
              input_file=args.input_file,
              output_folder=args.output_folder,
              clustering_method=args.clustering_method,
              n_clusters=args.n_clusters)    

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)