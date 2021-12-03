# -*- coding: utf-8 -*-
"""
This script is made to analyse the trajectories of an independent
AdaptivePELE simulation.
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

    parser.add_argument("-d", "--directory", type=str, dest = "input_folder",\
        default = 'LIG_Pele', help="Name of the directory where the simulation\
        is located.")
    parser.add_argument("-rn", "--report_name", type=str, dest = "report_name",\
        default = 'report', help="Name of the reports to be analyzed.")

    parsed_args = parser.parse_args(args)

    return parsed_args

def path_definer(input_folder):
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
    path_previous_simulation = os.path.join(path,input_folder)
    path_output = os.path.join(path_previous_simulation,'output')

    if os.path.isdir(path_previous_simulation) == False:
        raise Exception('PathError: There is no folder with this name: ' + path_previous_simulation + '. Please check the path and the folder name.')
    elif os.path.isdir(path_output) == False:
        raise Exception('NoOutputError: There is no output folder inside ' + path_previous_simulation + '. Please check the path and the folder name.')

    path_reports = os.path.join(path,'/HYB_analysis')

    if  os.path.exists(path_reports) == False:
        os.mkdir(path_reports)

    return path_output, path_reports

def backtracker(path_output,
                path_reports,
                report_name):
    """
    Function
    ----------
    For each trajectory in the last epoch of the simulation, this function 
    backtracks all the correponding reports in each epoch and copies them
    in a folder: one for each trajectory. 

    Parameters
    ----------
    - path_output : str
        The path to the directory where the raw data is located.
    - path_reports : str 
        The path where the copied reports are goig to be located.
    - report_name : str 
        Name the reports you want to copy have.
    """   

    def data_sorter(path):
        """
        Function
        ----------
        Reads the data from processorMap.txt files and sorts it out.

        Parameters
        ----------
        - path : str
            The path to the directory where the raw data is located.

        Returns
        ----------
        - matrix_epoch : numpy.array
            Matrix of dimensions: 4 x Number of CPUS, where the four dimensions correspond to the 
            values: epoch, trajectory, snapshot and cpu.
        """    

        epochs = []
        trajectories = []
        snapshots = []
        cpus = []

        procMapping = open(os.path.join(path,"processorMapping.txt")).read().rstrip().split(':')
        procMapping = [eval(x) for x in procMapping]

        for tuple_number in range(len(procMapping)):

            epoch, trajectory, snapshot = procMapping[tuple_number][:]

            trajectories.append(trajectory)
            epochs.append(epoch)
            snapshots.append(snapshot)
            cpus.append(tuple_number + 1)

        matrix_epoch = np.vstack((np.array(epochs), np.array(trajectories), np.array(snapshots), np.array(cpus)))

        return matrix_epoch
    
    def tracker(matrix,epoch,trajectory):
        """
        Function
        ----------
        Retrieves the information about the trajectory in the previous epoch
        from which the current trajectory in the current epoch has begun.

        Parameters
        ----------
        - matrix : numpy.array
            Matrix of dimensions: # of epochs x 4 x Number of CPUS, where the four dimensions correspond to the 
            values: epoch, trajectory, snapshot and cpu.
        - epoch : int
            Epoch from which the backtrack begins.
        - trajectory : int
            Trajectory from which the backtrack begins.

        Returns
        ----------
        - previous_epoch : int
            Epoch from which the initial conformation of the current trajectory has begun. 
        - previous_trajectory : int
            Trajectory from which the initial conformation of the current trajectory has begun. 
        """

        previous_epoch, previous_trajectory, _, _ = matrix[epoch-1,:,trajectory-1] 
        # -1 is added due to python's indexing

        return previous_epoch, previous_trajectory       

    all_directories_in_output = os.listdir(path_output)
    unsorted_numeric_directories = [int(i) for i in all_directories_in_output if i.isnumeric()]
    sorted_directories = sorted(unsorted_numeric_directories)
    directories = [str(i) for i in sorted_directories]

    matrix_simulation_list = []

    # Reading data ---

    for directory in directories:

        if directory != '0':
            
            matrix_epoch = data_sorter(os.path.join(path_output,directory)) 
            matrix_simulation_list.append(matrix_epoch) 

    matrix_simulation = np.array(matrix_simulation_list)

    # ---

    # Finding trajectories ---

    number_of_trajectories = matrix_simulation.shape[2]

    for trajectory in range(1,number_of_trajectories + 1):

        path_trajectory = os.path.join(path_reports,str(trajectory))

        if  os.path.exists(path_trajectory) == False:
            os.mkdir(path_trajectory)

        epoch = matrix_simulation.shape[0]

        # Copying initial report
        shutil.copy(os.path.join(path_output,str(epoch),'report_' + str(trajectory)), path_trajectory)

        while epoch != 0:

            epoch, trajectory = tracker(matrix_simulation,epoch,trajectory)  

            # Copying backtracked reports
            shutil.copy(os.path.join(path_output,str(epoch),'report_' + str(trajectory)), path_trajectory)

        print(' ')

def main(args):
    """
    Function
    ----------
    It runs path definer to then run the backtracker function. 

    Parameters
    ----------
    - args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    # Defining paths
    path_output, path_reports = path_definer(input_folder=args.input_folder)

    # Tracking reports
    backtracker(path_output,
                path_reports,
                report_name = args.report_name)

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)





    
    