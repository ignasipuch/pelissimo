import sys
import os
import pathlib 
import argparse
import shutil
import time
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
    parser.add_argument("-cl", "--clusters_folder", type=str, dest = "clusters_folder",\
        default = 'results', help="Name of the directory containing the folder: clusters.")

    parsed_args = parser.parse_args(args)

    return parsed_args

labels = ['A','B','C','D','E','F','G','H',]

def linen_correction(input_folder, residue_name, clusters_folder):
    """
    It prepares everything to perform a PELE energy calculation 
    of all the clustered positions of a previous induced fit 
    PELE simulation.

    Parameters
    ----------
    input_folder : str
        The path to the directory created by the induced fit simulation.
    residue_name : str
        Residue name of the ligand in the pdb of each cluster.
    clusters_folder : str
        Name of the directory where the directory clusters is located (results/analysis).
    conf_file_name: str
        Name of the .conf file used to run pele (not adaptive.conf).
    """

    def path_definer(input_folder,
                     residue_name,clusters_folder):

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
        path_results = path_previous_simulation + '/' + clusters_folder
        path_output = path_previous_simulation + '/output' 

        if os.path.isdir(path_previous_simulation) == False:
            raise Exception('PathError: There is no folder with this name: ' + path_previous_simulation + '. Please check the path and the folder name.')

        path_energies = path + '/' + residue_name + '_linen'
        path_energies_input = path_energies + '/' + 'input'
        path_energies_simulation = path_energies + '/' + 'simulation'

        if  os.path.exists(path_energies) == False:
            os.mkdir(path_energies)

        if  os.path.exists(path_energies_input) == False:
            os.mkdir(path_energies_input)

        if  os.path.exists(path_energies_simulation) == False:
            os.mkdir(path_energies_simulation)

        return  path_previous_simulation, path_output, path_results, path_energies_input, path_energies_simulation
        
    _, path_output, path_results, _, _ = path_definer(input_folder,residue_name,clusters_folder)

    cont = 0

    step = []
    path = []
    cluster = []
    report_paths = []

    for i in range(len(labels)):
        labels[i] = str(i)

    with open(os.path.join(path_results,'data.csv')) as filein:

        for line in filein:

            if cont != 0:
      
                line = line.split(',')      
                line_1 = line[-1].split()[0]

                if any(label == line_1 for label in labels):

                    path_string = line[-2].split('trajectory_')[0]

                    report_num = str(line[-2])
                    report_num = report_num.split('/')[-1]  
                    report_num = report_num.split('.pdb')[0].split('_')[1]

                    step.append(int(line[0]))   
                    cluster.append(int(line[-1].split('\n')[0]))
                    report_paths.append(os.path.join(path_string,'report_' + str(report_num)))                    

            cont += 1

    report_paths_dictionary = Counter(report_paths)
    report_paths = sorted(report_paths)

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

    linen_correction(input_folder = args.input_folder,
                     residue_name = args.residue_name,
                     clusters_folder = args.clusters_folder)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)