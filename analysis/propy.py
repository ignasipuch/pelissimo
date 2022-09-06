# -*- coding: utf-8 -*-
"""
This script clusters trajectories depending on dihedral information
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import os
import pathlib
import argparse
from collections import defaultdict
import concurrent.futures
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    parser.add_argument("-cm", "--clustering_method", type=str, dest="clustering_method",
                        default='bin', help="Method to cluster data: bin or kmeans.")
    parser.add_argument("-nc", "--n_clusters", type=int, dest="n_clusters",
                        default=0, help="Number of clusters to cluster the data.")

    parser.add_argument("--evolution", dest="evolution_bool",
                        default=False, action='store_true', help="Flag to choose if dihedral evolution is wanted.")
    parser.add_argument("--cpus", type=int, dest="cpus", help="Flag to choose number of cpus.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def trajectory_processer(path,file,atom_list):

    def dihedral_angle_calculator(coordinates_atom_a,
                                  coordinates_atom_b,
                                  coordinates_atom_c,
                                  coordinates_atom_d):
        """
        Function
        ----------
        Calculate the dihedral angles of a conformation.
        Parameters
        ----------
        - dihedral_bond_dict : dict
            Dictionary with the atoms involved in the calculation of the dihedrals.
        - atoms_positions_dict : dict
            Dictionary with all the atoms of interest and their x, y, z positions.
        Returns
        ----------
        - dihedral_angle_dict : dict
            Dictionary with dihedral angles information of that conformation.
        """

        b1 = coordinates_atom_b - coordinates_atom_a
        b2 = coordinates_atom_c - coordinates_atom_b
        b2_u = b2/np.linalg.norm(b2)
        b3 = coordinates_atom_d - coordinates_atom_c

        n1 = np.cross(b1, b2)/np.linalg.norm(np.cross(b1, b2))
        n2 = np.cross(b2, b3)/np.linalg.norm(np.cross(b2, b3))

        m1 = np.cross(n1, b2_u)

        y = np.dot(m1, n2)
        x = np.dot(n1, n2)

        dihedral_angle = -np.arctan2(y, x)*(180./np.pi)
        
        return dihedral_angle

    dihedral_angles_trajectory = {}

    trajectory_number = file.split('_')[-1].split('.pdb')[0]

    filein = open(os.path.join(path, file), "r")
    content = filein.read()
    models_list = content.split('MODEL')
    models_list = models_list[1:len(models_list)]

    model_cont = 1

    for model in models_list:

        amino_dict = {}
        angles_dict = {}
        atoms_positions_dict = {}

        model_lines = model.splitlines(True)
        ligand_lines = [
            line for line in model_lines if 'ATOM' in line]

        cont_atoms = 0
        cont_amino = 0
        cont_lines = 0

        # Retrieving location of all the important atoms
        # and number of aminoacids

        for line in ligand_lines:
            if any(ext in line for ext in atom_list):
                line = line.split()

                if len(line[4]) == 5:
                    atoms_positions_dict[str(cont_atoms%3)] = [
                        float(line[5]), float(line[6]), float(line[7])]
                    cont_atoms +=1 

                else:
                    atoms_positions_dict[str(cont_atoms%3)] = [
                        float(line[6]), float(line[7]), float(line[8])]
                    cont_atoms +=1 

                if cont_atoms%3 == 0:
                    cont_amino += 1
                    
                    if cont_lines != 0:
                        key = str(int(cont_atoms/3))
                        amino_dict[key] = atoms_positions_dict.copy()
            
            cont_lines += 1 
        
        # From aminoacids and locations calculating all psi and phi 
        # angles for all conformations.

        for i in range(1,len(amino_dict)-1): 

            ni = np.array(amino_dict[str(i)][str(0)])
            cai = np.array(amino_dict[str(i)][str(1)])
            ci = np.array(amino_dict[str(i)][str(2)])
            
            ni1 = np.array(amino_dict[str(i + 1)][str(0)])
            cai1 = np.array(amino_dict[str(i + 1)][str(1)])
            ci1 = np.array(amino_dict[str(i + 1)][str(2)])

            psi = dihedral_angle_calculator(ni,cai,ci,ni1)
            phi = dihedral_angle_calculator(ci,ni1,cai1,ci1)

            angles_dict[i-1] = [psi,phi]

        dihedral_angles_trajectory[model_cont] = angles_dict

        model_cont += 1

    return trajectory_number, dihedral_angles_trajectory


def dihedral_angles_retriever_main(input_folder,
                                   cpus):
    """
    Function
    ----------
    Calculates the important dihedrals of all the conformations in a PELE 
    simulation. 

    Parameters
    ----------
    - input_folder : str
        Name of the folder where the output of the simulation is located.
    - residue_name : str
        Residue name of the ligand.
    - input_file : str
        Name of the file corresponding to the isolated ligand with connectivity.

    Returns
    ----------
    - dihedral_bond_df : pd.DataFrame
        Data frame with rotatable bonds, atoms conforming it and the index assigned.
    - simulation_df : pd.DataFrame
        Data frame with all the rotatable bonds' dihedral angle values of all the simulation
        with corresponding model, trajectory and epoch.
    - path_results : str 
        Path to the directory where the results will be stored.
    """

    def path_definer(input_folder):
        """
        Function
        ----------
        Defines the important paths that are going to be used throughout the simulation.

        Parameters
        ----------
        - input_folder : str
            Name of the folder where the output of the simulation is located.

        Returns
        ----------
        - path_template : str
            Path to the rotatable bonds file.
        - path_output : str
            Path to the output folder of the simulation.
        - path_results : str 
            Path to the directory where the results will be stored.
        """

        path = str(pathlib.Path().absolute())
        path_template = os.path.join(path, 'DataLocal', 'LigandRotamerLibs')
        path_output = os.path.join(path, input_folder)
        path_results = os.path.join(path, 'dihedrals')
        path_images = os.path.join(path_results,'images')

        if os.path.exists(path_results) == False:
            os.mkdir(path_results)

        if os.path.exists(path_images) == False:
            os.mkdir(path_images)

        return path, path_template, path_output, path_results, path_images

    def residency_function(path_output,
                           path):
        """
        Function
        ----------
        Retrieve the information about residency of each conformation. 

        Parameters
        ----------
        - path_output : str
            Path to the output folder of the simulation.

        Returns
        ----------
        - residency_df : pd.DataFrame
            Data Frame with the residency information of each model in the simulation.
        """

        def pelesteps_retriever(path):
            """
            Function
            ----------
            Reads the pelesteps from adaptive.conf or pele.conf.

            Returns
            ----------
            - pele_steps : int
                Number of pele steps of the simulation(s).
            """

            adaptive_path = os.path.join(path, 'adaptive.conf')

            if os.path.isfile(adaptive_path) == True:

                #
                print('     -   Retrieving information from adaptive.conf.')
                #

                with open(adaptive_path) as filein:

                    for line in filein:

                        if "peleSteps" in line:

                            peleSteps_string = line.split()[2]
                            pele_steps = int(peleSteps_string.split(',')[0])

            else:

                peleconf_path = os.path.join(path, 'pele.conf')

                #
                print('     -   Retrieving information from pele.conf.')
                #

                if os.path.isfile(peleconf_path) == True:

                    with open(peleconf_path) as filein:

                        for line in filein:

                            if "numberOfPeleSteps" in line:

                                pele_steps = int(line.split()[-1])

                else:

                    pele_steps = None

                    #
                    print('     -   No .conf was found.')
                    print('     -   The step weighted scoring function will not \n'
                          '         be taken into account.')
                    #

            return pele_steps

        def residency_retriever(path,
                                file_list,
                                pele_steps):
            """
            Function
            ----------
            Calculate the dihedral bond angles for an epoch.

            Parameters
            ----------
            - path : str
                Path to an epoch of the simulation.
            - file_list : list 
                List of all the files in directory.
            - atom_list : list
                List with all the atoms without repetition that have to be tracked.
            - dihedral_bond_dict : dict
                Dictionary with the atoms involved in the calculation of the dihedrals.

            Returns
            ----------
            - dihedral_angles_epoch : dict
                Dictionary with dihedral angles information of that epoch.
            """

            residency_epoch = {}
            file_list = [k for k in file_list if k.startswith('report')]

            for file in file_list:

                initial_step_list = []
                final_step = []
                model = []

                if file.startswith('report'):

                    trajectory_number = int(file.split('_')[-1])
     
                    with open(os.path.join(path,file)) as filein:

                        cont = 0 

                        for line in filein:

                            if cont != 0:

                                line = line.split()
                                model.append(cont)
                                initial_step_list.append(int(line[1]))

                            cont += 1

                    final_step = np.array(initial_step_list[1:])
                    final_step = np.append(final_step,pele_steps)
                    initial_step = np.array(initial_step_list)
                    residency = final_step - initial_step
                    residency = residency.astype(int)

                    if residency[-1] == 0: residency[-1] = 1

                    zip_iterator = zip(model, residency)
                    residency_epoch[trajectory_number] = dict(zip_iterator)  

            return residency_epoch

        pele_steps = pelesteps_retriever(path)
        files = os.listdir(path_output)
        numeric_files = [s for s in files if s.isnumeric()]

        #
        print('     -   Calculating residency.')
        #        

        residency_dict = {}

        if len(numeric_files) != 0:

            for epoch in numeric_files:

                new_directory = os.path.join(path_output, epoch)

                if os.path.isdir(new_directory) and epoch.isnumeric():

                    files = os.listdir(new_directory)
                    residency_dict[int(epoch)] = residency_retriever(new_directory,
                                                                     files,
                                                                     pele_steps)

        else:

            residency_dict[0] = residency_retriever(path_output,
                                                    files,
                                                    pele_steps)

        residency_df = pd.DataFrame([(epoch, trajectory, model, residency)
                                      for epoch, traj_mod_residency in residency_dict.items()
                                      for trajectory, mod_residency in traj_mod_residency.items()
                                      for model, residency in mod_residency.items()])

        residency_df.columns = ['epoch', 'trajectory',
                                 'model', 'residency']

        return residency_df        

    def trajectory_positions(path_output,
                             atom_list,
                             cpus):
        """
        Function
        ----------
        Calculate the dihedral bond angles for all the conformation in all epochs, 
        trajectories and models.

        Parameters
        ----------
        - path_output : str
            Path to the output folder of the simulation.
        - atom_list : list
            List with all the atoms without repetition that have to be tracked.
        - dihedral_bond_dict : dict
            Dictionary with the atoms involved in the calculation of the dihedrals.

        Returns
        ----------
        - simulation_df : pd.DataFrame
            Data frame with all the rotatable bonds' dihedral angle values of the 
            all the simulation with corresponding model, trajectory and epoch.
        """

        def trajectory_positions_retriever(path,
                                           file_list,
                                           atom_list):
            """
            Function
            ----------
            Calculate the dihedral bond angles for an epoch.

            Parameters
            ----------
            - path : str
                Path to an epoch of the simulation.
            - file_list : list 
                List of all the files in directory.
            - atom_list : list
                List with all the atoms without repetition that have to be tracked.
            - dihedral_bond_dict : dict
                Dictionary with the atoms involved in the calculation of the dihedrals.

            Returns
            ----------
            - dihedral_angles_epoch : dict
                Dictionary with dihedral angles information of that epoch.
            """

            dihedral_angles_epoch = {}
            file_list = [k for k in file_list if 'trajectory' in k]

            start_time = time.perf_counter()

            with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
                results = [executor.submit(trajectory_processer,path,file,atom_list) for file in file_list]

                for f in concurrent.futures.as_completed(results):
                    trajectory_number, dihedral_angles_trajectory = f.result()
                    dihedral_angles_epoch[int(trajectory_number)] = dihedral_angles_trajectory

            final_time = time.perf_counter()
            print('     -   Time to retrieve data: ' + str(final_time - start_time))

            return dihedral_angles_epoch

        files = os.listdir(path_output)
        numeric_files = [s for s in files if s.isnumeric()]

        simulation_dict = {}

        #
        print('     -   Retrieving trajectory data...')
        #  

        if len(numeric_files) != 0:

            for epoch in numeric_files:

                new_directory = os.path.join(path_output, epoch)

                if os.path.isdir(new_directory) and epoch.isnumeric():

                    files = os.listdir(new_directory)

                    simulation_dict[int(epoch)] = trajectory_positions_retriever(new_directory,
                                                                                 files,
                                                                                 atom_list)

        else:

            simulation_dict[0] = trajectory_positions_retriever(path_output,
                                                                files,
                                                                atom_list)

        simulation_df = pd.DataFrame([(epoch, trajectory, model, rot_bond, value)
                                      for epoch, traj_mod_rot_val in simulation_dict.items()
                                      for trajectory, mod_rot_val in traj_mod_rot_val.items()
                                      for model, rot_val in mod_rot_val.items()
                                      for rot_bond, value in rot_val.items()])

        simulation_df.columns = ['epoch', 'trajectory',
                                 'model', 'rotatable bond', 'value']

        return simulation_df

    def residency_to_simulation(residency_df,
                                simulation_df):
        """
        Function
        ----------
        Add the residency information to the simulation data frame.

        Parameters
        ----------
        - residency_df : pd.DataFrame
            Data Frame with the residency information of each model in the simulation.
        - simulation_df : pd.DataFrame
            Data frame with all the rotatable bonds' dihedral angle values of the 
            all the simulation with corresponding model, trajectory and epoch.

        Returns
        ----------
        - simulation_df : pd.DataFrame
            Data frame with all the rotatable bonds' dihedral angle values with residency
            information.
        """

        print(' -   Adding residency information to trajectories.')

        for _, row in residency_df.iterrows():
            epoch = row['epoch']
            trajectory = row['trajectory']
            model = row['model']
            residency = row['residency']

            epoch_df = simulation_df.loc[(simulation_df['epoch'] == epoch)]
            trajectory_df = epoch_df.loc[(epoch_df['trajectory'] == trajectory)]
            row_df = trajectory_df.loc[(trajectory_df['model'] == model)]
            
            if residency - 1 >= 1:  
                for _ in range(residency-1):
                    simulation_df = simulation_df.append(row_df,ignore_index=True)

        return simulation_df

    path, _, path_output, path_results, path_images = path_definer(input_folder)

    residency_df = residency_function(path_output,
                                      path)
                                    
    atom_list = ['  N   ','  CA  ','  C   ']

    simulation_df = trajectory_positions(path_output,
                                         atom_list,
                                         cpus)
    
    start_time = time.perf_counter()
    full_df = residency_to_simulation(residency_df,
                                      simulation_df)
    final_time = time.perf_counter()
    print(' -   Residency to simulation time: ', final_time - start_time)

    full_df = simulation_df

    return full_df, path_results, path_images


def clustering(simulation_df,
               path_results,
               path_images):
    """
    Function
    ----------
    Cluster the results obtained and stored in a data frame.

    Parameters
    ----------
    - n_cluster : int
        Number of clusters to cluster data.
    - clustering_method : str
        Method to cluster data: bin or kmeans.
    - simulation_df : pd.DataFrame
        Data frame with all the rotatable bonds' dihedral angle values of all the simulation
        with corresponding model, trajectory and epoch.
    - dihedral_bond_df : pd.DataFrame
        Data frame with rotatable bonds, atoms conforming it and the index assigned.
    - path_results : str 
        Path to the directory where the results will be stored. 

    """

    def binning(simulation_df,
                path_images,
                residency_df):
        """
        Function
        ----------
        Cluster the results obtained from simulation by binning the results.
        Entropic contributions are calculated from the binned data.

        Parameters
        ----------
        - simulation_df : pd.DataFrame
            Data frame with all the rotatable bonds' dihedral angle values of all the simulation
            with corresponding model, trajectory and epoch.
        - path_images : str 
            Path to the directory where the results will be stored. 
        """

        entropy_contribution_psi = []
        entropy_contribution_phi = []

        rotatable_bonds = simulation_df['rotatable bond'].to_numpy()
        values_psi = simulation_df['value'].to_numpy()[0]
        values_phi = simulation_df['value'].to_numpy()[1]

        results_psi = defaultdict(list)
        results_phi = defaultdict(list)

        for rot_bond, value in zip(rotatable_bonds, values_psi):
            results_psi[rot_bond].append(value)

        for rot_bond, value in zip(rotatable_bonds, values_phi):
            results_phi[rot_bond].append(value)

        rot_bond_values_psi = list(results_psi.items())
        rot_bond_values_phi = list(results_phi.items())
     
        for rot_bond, values in rot_bond_values_psi:

            bin_edges = np.histogram_bin_edges(values, bins=10)
            density, _ = np.histogram(
                values, bins=bin_edges, density=True)
            dense_bins = density[density != 0]

            entropy_contribution_psi.append(
                np.sum(np.array([p*np.log(p) for p in dense_bins])))

            # Plot
            plt.title('Dihedral ' + str(rot_bond) + ' distribution')
            plt.hist(values,
                     bins=bin_edges, density=True)
            plt.xlabel('Dihedral angle (º)')
            plt.ylabel('Density')
            plt.xlim(-180, 180)
            plt.xticks(list(np.arange(-180, 190, 30)))
            plt.savefig(os.path.join(path_images, 'dihedral_' +
                                     str(rot_bond) + '_distribution_psi.png'), format='png', transparent=True)
            plt.close()

        for rot_bond, values in rot_bond_values_phi:

            bin_edges = np.histogram_bin_edges(values, bins=10)
            density, _ = np.histogram(
                values, bins=bin_edges, density=True)
            dense_bins = density[density != 0]

            entropy_contribution_phi.append(
                np.sum(np.array([p*np.log(p) for p in dense_bins])))

            # Plot
            plt.title('Dihedral ' + str(rot_bond) + ' distribution')
            plt.hist(values,
                     bins=bin_edges, density=True)
            plt.xlabel('Dihedral angle (º)')
            plt.ylabel('Density')
            plt.xlim(-180, 180)
            plt.xticks(list(np.arange(-180, 190, 30)))
            plt.savefig(os.path.join(path_images, 'dihedral_' +
                                     str(rot_bond) + '_distribution_phi.png'), format='png', transparent=True)
            plt.close()

        entropy_contributions = np.array(entropy_contribution)
        S = -(R)*np.sum(entropy_contributions)

        entropy_percentages = 100. * \
            np.array(entropy_contribution)/np.sum(entropy_contributions)
        entropy_df = pd.DataFrame(entropy_percentages)
        entropy_df = entropy_df.round(decimals=2).T
        entropy_df.columns = [
            'Dihedral_' + str(value) + '_%' for value in list(np.arange(1, len(rot_bond_values) + 1))]
        entropy_df.insert(loc=0, column='S (kcal/mol K)',
                          value="{:3.6f}".format(S))
        entropy_df.to_csv(os.path.join(
            path_results, 'entropy.csv'), index=False)

        print(' -   Entropic information written in /dihedrals/entropy.csv.')

    print(' -   Clustering data obteined...')

    binning(simulation_df,
            path_images)

    simulation_df.to_csv(os.path.join(path_results, 'data.csv'), index=False)


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
    #
    print(' ')
    print('* Dihedral Clustering *')
    print(' ')
    print(' -   Gathering information')
    #

    full_df, \
    path_results,     \
    path_images = dihedral_angles_retriever_main(input_folder=args.input_folder,
                                                 cpus = args.cpus)

    #
    print(' ')
    print(' -   Beginning of clustering')
    #

    start_time = time.perf_counter()
    clustering(full_df,
               path_results,
               path_images)
    final_time = time.perf_counter()
    print(' -   Clustering time: ', final_time - start_time)

    #
    print(' ')
    print(' -   Results written in /dihedrals.\n')


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    start_time = time.perf_counter()
    main(args)
    final_time = time.perf_counter()
    print(' -   Total time: ', final_time - start_time)