# -*- coding: utf-8 -*-
"""
This script ensambles other scripts and calculates entropy
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

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
                        default=None, help="Name of the file corresponding to the isolated ligand with connectivity.")
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
    parser.add_argument("-d2", "--second_directory", type=str, dest="second_directory",
                        default='linen', help="Name of the second folder you want the entropic information from.")
    parser.add_argument("--cpus", type=int, dest="cpus",
                        help="Flag to choose number of cpus.")

    parser.add_argument("--evolution", dest="evolution_bool",
                        default=False, action='store_true', help="Flag to choose if dihedral evolution is wanted.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def ensambler(input_folder,
              second_directory,
              residue_name,
              input_file,
              output_folder,
              clustering_method,
              n_clusters,
              evolution_bool,
              cpus):

    path = str(pathlib.Path().absolute())
    path_pl_simulation = os.path.join(path, input_folder)
    path_l_simulation = os.path.join(path, residue_name + f'_{second_directory}')

    if second_directory == 'linen':

        os.chdir(path_pl_simulation)
        os.system('python /home/bsc72/bsc72825/projects/code/dihedral_clustering.py -f ' + input_file + ' -d ' +
                  output_folder + ' -r ' + residue_name + ' -cm ' + clustering_method + ' -nc ' + str(n_clusters))
        os.chdir(path_l_simulation)

        if evolution_bool:
            os.system('python /home/bsc72/bsc72825/projects/code/dihedral_clustering.py -f ' + input_file + ' -d ' +
                  output_folder + ' -r ' + residue_name + ' -cm ' + clustering_method + ' -nc ' + str(n_clusters) + ' --evolution')
        else: 
            os.system('python /home/bsc72/bsc72825/projects/code/dihedral_clustering.py -f ' + input_file + ' -d ' +
                  output_folder + ' -r ' + residue_name + ' -cm ' + clustering_method + ' -nc ' + str(n_clusters))

        os.chdir(path)
        os.system('python /home/bsc72/bsc72825/projects/code/lice.py -d ' +
                  input_folder + ' -r ' + residue_name + f'-d2 LIG_{second_directory}')


    if second_directory == 'prot':

        os.chdir(path_pl_simulation)
        os.system('python /home/bsc72/bsc72825/projects/code/propy.py -d ' +
                  output_folder + ' -r ' + residue_name + f' --cpus {cpus}')

        os.chdir(path_l_simulation)
        os.system('python /home/bsc72/bsc72825/projects/code/propy.py -d ' +
                  output_folder + ' -r ' + residue_name + f' --cpus {cpus}')

        os.chdir(path)
        os.system('python /home/bsc72/bsc72825/projects/code/lice.py -d ' +
                  input_folder + ' -r ' + residue_name + f' -d2 LIG_{second_directory}')


def main(args):

    ensambler(input_folder=args.input_folder,
              second_directory=args.second_directory,
              residue_name=args.residue_name,
              input_file=args.input_file,
              output_folder=args.output_folder,
              clustering_method=args.clustering_method,
              n_clusters=args.n_clusters,
              evolution_bool=args.evolution_bool,
              cpus=args.cpus)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
