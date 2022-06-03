# -*- coding: utf-8 -*-
"""
This script is designed to obtain the conformational entropy contribution
to the free energy change by using the results obtained from the dihedral 
clustering.
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import os
import pathlib
import argparse
from distutils.dir_util import copy_tree

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
    parser.add_argument("-r", "--residue_name", type=str, dest="residue_name",
                        default='LIG', help="Ligand's residue name.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def lice_results(input_folder,
                 residue_name):
    """
    Function
    ----------
    It calculates the change entropy given a set of clusters with their populations.

    Parameters
    ----------
    - input_folder : str
        The path to the directory created by the induced fit simulation.
    - residue_name : str
        Residue name of the ligand in the pdb of each cluster.

    """

    def path_definer(input_folder,
                     residue_name):
        """
        Function
        ----------
        Defines all the paths that are going to be used

        Parameters
        ----------
        - input_folder : str
            The path to the directory created by the induced fit simulation.
        - residue_name : str
            Residue name of the ligand in the pdb of each cluster.

        Returns
        ----------
        - path_pl_simulation : str
            The path to the induced fit simulation.
        - path_pl_dihedrals : str
            The path to the clustering folder generated by the platform of the 
            protein-ligand simulation.
        - path_l_dihedrals : str
            The path to the clustering folder generated by the platform of the 
            ligand simulation.

        """

        path = str(pathlib.Path().absolute())

        path_pl_simulation = os.path.join(path, input_folder)
        path_pl_dihedrals = os.path.join(
            path_pl_simulation, 'dihedrals')

        path_l_simulation = os.path.join(path, residue_name + '_linen')
        path_l_dihedrals = os.path.join(
            path_l_simulation, 'dihedrals')
        path_entropy = os.path.join(path_pl_simulation,'entropy')

        if os.path.isdir(path_pl_simulation) == False:
            raise Exception('PathError: There is no folder with this name: ' +
                            path_pl_simulation + '. Please check the path and the folder name.')

        if os.path.exists(path_entropy) == False:
            os.mkdir(path_entropy)

        return path_pl_simulation, path_pl_dihedrals, path_l_dihedrals, path_entropy

    def entropy_retriever(path):
        """
        Function
        ----------
        Retireves entropy of the simulation from dihedral clustering results.

        Parameters
        ----------
        - path : str
            Path to the entropy.csv location.

        Returns
        ----------
        - S : float
            Conformational entropy of the simulation.
        """

        cont = 0

        with open(os.path.join(path, 'entropy.csv')) as filein:

            for line in filein:

                if cont != 0:

                    line = line.split(',')
                    S = float(line[0])

                cont += 1

        return S

    #
    print(' ')
    print('* LiCE *')
    print(' ')
    #

    path_pl_simulation,\
        path_pl_dihedrals,\
        path_l_dihedrals,\
        path_entropy = path_definer(
            input_folder, residue_name)

    if os.path.isdir(path_pl_dihedrals) == False:
        raise Exception('dihedralsFolderNotFound: ' +
                        path_pl_dihedrals + ' does not exist.')

    if os.path.isdir(path_l_dihedrals) == False:
        raise Exception('dihedralsFolderNotFound: ' +
                        path_l_dihedrals + 'does not exist.')

    S_out = entropy_retriever(path_l_dihedrals)
    S_in = entropy_retriever(path_pl_dihedrals)

    copy_tree(path_l_dihedrals, os.path.join(path_entropy,'dihedrals_lig'))
    copy_tree(path_pl_dihedrals, os.path.join(path_entropy,'dihedrals_if'))

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
    print(' -   Writing entropy.csv inside ' + input_folder + '/entropy.')
    print(' -   Folders containing intermediate results are also in ' + input_folder + '/entropy.\n'
          '         /dihedrals_lig -> ligand simulation\n'
          '         /dihedrals_if  -> induced fit simulation')
    print(' ')
    #

    with open(os.path.join(path_entropy, 'entropy.csv'), 'w') as fileout:
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
                 residue_name=args.residue_name)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
