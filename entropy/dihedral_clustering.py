# -*- coding: utf-8 -*-
"""
This cluster trajectories with dihedrals
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
                        default='output', help="Name of the output directory where the simulation\
        is located.")
    parser.add_argument("-r", "--residue_name", type=str, dest="residue_name",
                        default='LIG', help="Ligand's residue name.")
    parser.add_argument("-c", "--conf_name", type=str, dest="conf_name",
                        default='pele', help="Name of the conf file used for the simulation.")

    parsed_args = parser.parse_args(args)

    return parsed_args

def main_function(conf_name,
                  residue_name):

    def forcefield_retriever(conf_name):

        path = str(pathlib.Path().absolute())

        with open(os.path.join(path,conf_name + '.conf')) as filein:

            for line in filein:

                if "ForceField" in line:

                    line = line.split(':')
                    line = line[1].split('     ')
                    line = line[0].split('"')
                    forcefield = line[1]

        return forcefield


    def path_definer(forcefield):

        path = str(pathlib.Path().absolute())

        if forcefield == 'OPLS2005':

            path_template = os.path.join(path,'DataLocal','Templates',forcefield,'HeteroAtoms')

        elif forcefield == 'OpenFF-OPLS2005':

            path_template = os.path.join(path,'DataLocal','Templates','OpenFF','Parsley')

        else: 

            raise Exception('ForceFieldError: The force field read in the .conf is not valid: OPLS2005 or OpenFF-OPLS2005.')

        return path_template
        

    def template_info_retriever(path_template,
                                residue_name):

        with open(os.path.join(path_template,residue_name.lower() + 'z'))

           





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

    main_function(conf_name=args.conf_name)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)