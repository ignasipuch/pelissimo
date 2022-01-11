# -*- coding: utf-8 -*-
"""
This module is designed to run an analysis on PELE pltform resutls of a simulation.
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import os
import pathlib
import argparse
import numpy as np

# ----------------------------------------------------------------------- #
# Constants:
T = 298.
R = 1.985e-3
# ----------------------------------------------------------------------- #


def parse_args(args):
    """
    Function
    ----------
    It parses the command-line arguments.
    Parameters

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

    parser.add_argument("-d", "--directory", type=str, dest="input_folder",
                        default='LIG_Pele', help="Name of the directory where the simulation\
        is located.")
    parser.add_argument("-cl", "--clusters_folder", type=str, dest="clusters_folder",
                        default='results', help="Name of the directory containing the folder: clusters.")
    parser.add_argument("-r", "--residue_name", type=str, dest="residue_name",
                        default='LIG', help="Ligand's residue name.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def lice_results(input_folder,
                 clusters_folder,
                 residue_name):
    """
    Function
    ----------
    It calculates the change entropy given a set of clusters with their populations.

    Parameters
    ----------
    - input_folder : str
        The path to the directory created by the induced fit simulation.
    - clusters_folder : str
        Name of the directory where the directory clusters is located (results/analysis).

    """

    def path_definer(input_folder,
                     clusters_folder,
                     residue_name):
        """
        Function
        ----------
        Defines all the paths that are going to be used

        Parameters
        ----------
        - input_folder : str
            The path to the directory created by the induced fit simulation.
        - clusters_folder : str
            Name of the directory where the directory clusters is located (results/analysis).
        - residue_name : str
            Residue name of the ligand in the pdb of each cluster.

        Returns
        ----------
        - path_pl_clusters : str
            The path to the clustering folder generated by the platform of the 
            protein-ligand simulation.
        - path_l_clusters : str
            The path to the clustering folder generated by the platform of the 
            ligand simulation.

        """

        path = str(pathlib.Path().absolute())

        path_pl_simulation = os.path.join(path, input_folder)
        path_pl_clusters = os.path.join(
            path_pl_simulation, clusters_folder, 'clusters')

        path_l_simulation = os.path.join(path, residue_name + '_linen_cry')
        path_l_clusters = os.path.join(
            path_l_simulation, clusters_folder, 'clusters')

        if os.path.isdir(path_pl_simulation) == False:
            raise Exception('PathError: There is no folder with this name: ' +
                            path_pl_simulation + '. Please check the path and the folder name.')

        return path_pl_simulation, path_pl_clusters, path_l_clusters

    def entropy_calculator(path):
        """
        Function
        ----------
        Calculates entropy with a given path to a info.csv file generated by the platform.

        Parameters
        ----------
        - path : str
            Path to the info.csv location.

        Returns
        ----------
        - S : float
            Conformational entropy of the simulation.

        """

        cont = 0
        S_sum = 0.
        population = []
        cluster = []

        with open(os.path.join(path, 'info.csv')) as filein:

            for line in filein:

                if cont != 0:

                    line = line.split(',')
                    cluster.append(int(line[0]))
                    population.append(float(line[1]))

                cont += 1

            num_clusters = max(cluster)

            for i in range(num_clusters):
                S_sum += population[i]*np.log(population[i])

            S = -R*S_sum

        return S

    #
    print(' ')
    print('*******************************************************************')
    print('*                           peleLiCE                              *')
    print('* --------------------------------------------------------------- *')
    print('*        Ligand\'s conformational entropy from simulations        *')
    print('*******************************************************************')
    print(' ')
    #

    path_pl_simulation,\
        path_pl_clusters,\
        path_l_clusters = path_definer(
            input_folder, clusters_folder, residue_name)

    if os.path.isdir(path_pl_clusters) == False:
        raise Exception('ClusterFolderNotFound: ' +
                        path_pl_clusters + ' does not exist.')

    if os.path.isdir(path_l_clusters) == False:
        raise Exception('ClusterFolderNotFound: ' +
                        path_l_clusters + 'does not exist.')

    S_out = entropy_calculator(path_l_clusters)
    S_in = entropy_calculator(path_pl_clusters)
    DS = S_in - S_out
    DG_s = -T*DS

    #
    print(' -   Entropy unbound:', S_out, ' kcal/(mol K)')
    print(' -   Entropy bound:', S_in, ' kcal/(mol K)')
    print(' -   Entropy change:', DS, ' kcal/(mol K)')
    print(' ')
    print(' -   Temperature: ', T)
    print(' -   Entropic contribution to the free energy: ',
          DG_s, ' kcal/mol')
    print(' ')
    print(' -   Writing entropy.csv inside ' + input_folder + '.')
    print(' ')
    #

    with open(os.path.join(path_pl_simulation, 'entropy.csv'), 'w') as fileout:
        fileout.writelines(
            'S_out,S_in,DS,-TDS\n'
            '' + str(S_out) + ',' + str(S_in) + ',' +
            str(DS) + ',' + str(DG_s) + '\n'
        )


def main(args):
    """
    Function
    ----------
    It reads the command-line arguments and runs lice_results.

    Parameters
    ----------
    - args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    lice_results(input_folder=args.input_folder,
                 clusters_folder=args.clusters_folder,
                 residue_name=args.residue_name)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)