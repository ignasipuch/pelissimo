# -*- coding: utf-8 -*-
"""
This script is designed to fetch snapshots.
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

    parser.add_argument("-o", "--output_folder", type=str, dest="output_folder",
                        default='output', help="Name of the output directory where the simulation\
        is located.")
    parser.add_argument("-s", "--save", type=str, dest="save_name",
                        default=None, help="Name of the pdb that is saved.")

    parser.add_argument("-e", "--epoch", required=True, type=str, dest="epoch",
                        help="Epoch where the snapshot is located.")
    parser.add_argument("-t", "--trajectory", required=True, type=str, dest="trajectory",
                        help="Trajectory where the snapshot is located.")
    parser.add_argument("-m", "--model", required=True, type=str,  dest="model",
                        help="Model of the snapshot.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def fetcher(output_folder,
            epoch,
            trajectory,
            model,
            save_name):
    """
    Function
    ----------
    Retrieves the desired snapshot of a simulation.

    Parameters
    ----------
    - output_folder : str
        Path to the output folder where the interesting snapshot is located.
    - epoch : str
        Epoch where the snapshot is located.
    - trajectory : str
        Trajectory where the snapshot is located.
    - model : str
        Model of the snapshot.
    - save_name : str
        Name of the final pdb file.
    """

    def path_definer(output_folder):
        """
        Function
        ----------
        Defines the important paths that are going to be used throughout the simulation.

        Parameters
        ----------
        - output_folder : str
            Name of the folder where the output of the simulation is located.

        Returns
        ----------
        - path_output : str
            Path to the output folder of the simulation.
        - path_results : str
            Path to the results folder of the bootstrap analysis.
        """

        path = str(pathlib.Path().absolute())
        path_output = os.path.join(path, output_folder)
        path_results = os.path.join(path, 'fetched')

        if os.path.isdir(path_results) is False:
            os.mkdir(path_results)

        return path_output, path_results

    def retriever(path_output,
                  path_results,
                  epoch,
                  trajectory,
                  model,
                  save_name):
        """
        Function
        ----------
        Retrieves the wanted snapshot and saves it.

        Parameters
        ----------
        - path_output : str
            Path to the output folder where the interesting snapshot is located.
        - path_results : str
            Path to the result folder where the snapshot will be located.
        - epoch : str
            Epoch where the snapshot is located.
        - trajectory : str
            Trajectory where the snapshot is located.
        - model : str
            Model of the snapshot.
        - save_name : str
            Name of the final pdb file.
        """

        with open(os.path.join(path_output, epoch, 'trajectory_' + trajectory + '.pdb'), 'r') as filein:
            
            text = filein.readlines()

            begin_of_pdb = text.index('MODEL     '  + model + '\n') + 1 
            len_pdb = text.index('ENDMDL    \n') - 1
            end_of_pdb = begin_of_pdb + len_pdb

            output_file = open(os.path.join(path_results, save_name + '.pdb'), 'w')
            output_file.writelines(text[begin_of_pdb:end_of_pdb])


    print('\n* pele fetcher\n')

    path_output, path_results = path_definer(output_folder)

    if epoch.isnumeric() and trajectory.isnumeric() and model.isnumeric():

        if save_name is not None:
            retriever(path_output,
                      path_results,
                      epoch,
                      trajectory,
                      model,
                      save_name)
        else:
            save_name = epoch + '_' + trajectory + '_' + model
            retriever(path_output,
                      path_results,
                      epoch,
                      trajectory,
                      model,
                      save_name)

    else:

        epoch = epoch.strip('][').split(',')
        trajectory = trajectory.strip('][').split(',')
        model = model.strip('][').split(',')

        if save_name is not None:
            print(' -   Avoiding save name introduced due to overwriting problems')

        if len(epoch) == len(trajectory) and len(trajectory) == len(model):

            for i in range(len(epoch)):

                save_name = str(epoch[i]) + '_' + str(trajectory[i]) + '_' + str(model[i])
                retriever(path_output,
                          path_results,
                          epoch[i],
                          trajectory[i],
                          model[i],
                          save_name)
        
        else:
            raise Exception('LengthMismatch: Length mismatch between epoch, trajectory and model lists.')

    print(' -   Snapshot(s) stored in /fetched\n')




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

    fetcher(output_folder=args.output_folder,
            epoch=args.epoch,
            trajectory=args.trajectory,
            model=args.model,
            save_name=args.save_name)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
