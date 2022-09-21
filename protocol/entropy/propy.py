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
    parser.add_argument("--cpus", type=int, dest="cpus",
                        help="Flag to choose number of cpus.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def trajectory_processer(path, file, atom_list):
    """
    Function
    ----------
    Main loop to retrieve the positions of the important atoms to calculate
    the dihedral angles of psi and phi.

    Parameters
    ----------
    - path : str
        Dictionary with the atoms involved in the calculation of the dihedrals.
    - file : str
        PDB file of the trajectory from which to extract the data.
    - atom_list : list
        List of atoms that that are necessary to calculate both dihedral angles.

    Returns
    ----------
    - trajectory_number : int
        Number of trajectory that is analyzed.
    - dihedral_angles_trajectory : dict
        Dictionary with dihedral angles information of that trajectory.
    """

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
                    atoms_positions_dict[str(cont_atoms % 3)] = [
                        float(line[5]), float(line[6]), float(line[7])]
                    cont_atoms += 1

                else:
                    atoms_positions_dict[str(cont_atoms % 3)] = [
                        float(line[6]), float(line[7]), float(line[8])]
                    cont_atoms += 1

                if cont_atoms % 3 == 0:
                    cont_amino += 1

                    if cont_lines != 0:
                        key = str(int(cont_atoms/3))
                        amino_dict[key] = atoms_positions_dict.copy()

            cont_lines += 1

        # From aminoacids and locations calculating all psi and phi
        # angles for all conformations.

        for i in range(1, len(amino_dict)-1):

            ni = np.array(amino_dict[str(i)][str(0)])
            cai = np.array(amino_dict[str(i)][str(1)])
            ci = np.array(amino_dict[str(i)][str(2)])

            ni1 = np.array(amino_dict[str(i + 1)][str(0)])
            cai1 = np.array(amino_dict[str(i + 1)][str(1)])
            ci1 = np.array(amino_dict[str(i + 1)][str(2)])

            psi = dihedral_angle_calculator(ni, cai, ci, ni1)
            phi = dihedral_angle_calculator(ci, ni1, cai1, ci1)

            angles_dict[i-1] = [psi, phi]

        dihedral_angles_trajectory[model_cont] = angles_dict

        model_cont += 1

    return trajectory_number, dihedral_angles_trajectory


def dihedral_angles_retriever_main(input_folder,
                                   cpus):
    """
    Function
    ----------
    Calculates the important dihedrals of all the conformations in a PELE 
    simulation and stores the important information in three vectors: phi angles,
    psi angles and their occurrences in the simulation.

    Parameters
    ----------
    - input_folder : str
        Name of the folder where the output of the simulation is located.
    - cpus : int
        Number of cpus we want to retrieve the data.

    Returns
    ----------
    - weight_vector : np.array
        Vector with same length as psi_vector and phi_vector that informs about 
        the angles' weight, which is de residency.
    - psi_vector : np.array
        Vector with the psi values reached during the simulation.
    - phi_vector : np.array
        Vector with the psi values reached during the simulation.  
    - path_results : str 
        Path to the directory where the results will be stored.
    - path_images : str 
        Path to the directory where the images will be stored.
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
        - path: str
            Path to the working directory.
        - path_output : str
            Path to the output folder of the simulation.
        - path_results : str 
            Path to the directory where the results will be stored.
        - path images : str
            Path to the directory where the images are going to be stored.
        """

        path = str(pathlib.Path().absolute())
        path_output = os.path.join(path, input_folder)
        path_results = os.path.join(path, 'dihedrals')
        path_images = os.path.join(path_results, 'images')

        if os.path.exists(path_results) == False:
            os.mkdir(path_results)

        if os.path.exists(path_images) == False:
            os.mkdir(path_images)

        return path, path_output, path_results, path_images

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
        - path : str
            Path to the working directory.

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

            Parameters
            ----------
            - path : str
                Path to the working directory.

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
            - pele_steps : int
                Number of pele steps in the pele.conf.

            Returns
            ----------
            - residency_epoch : dict
                Dictionary with the information of residency of all the models.
            """

            residency_epoch = {}
            file_list = [k for k in file_list if k.startswith('report')]

            for file in file_list:

                initial_step_list = []
                final_step = []
                model = []

                if file.startswith('report'):

                    trajectory_number = int(file.split('_')[-1])

                    with open(os.path.join(path, file)) as filein:

                        cont = 0

                        for line in filein:

                            if cont != 0:

                                line = line.split()
                                model.append(cont)
                                initial_step_list.append(int(line[1]))

                            cont += 1

                    final_step = np.array(initial_step_list[1:])
                    final_step = np.append(final_step, pele_steps)
                    initial_step = np.array(initial_step_list)
                    residency = final_step - initial_step
                    residency = residency.astype(int)

                    if residency[-1] == 0:
                        residency[-1] = 1

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

        residency_df = residency_df.astype(
            {"epoch": int, "trajectory": int, "model": int, "residency": int})
        residency_df = residency_df.sort_values(['epoch', 'trajectory'])

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
        - cpus : int
            Number of cpus to retrieve the information from trajectories.

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

            Returns
            ----------
            - dihedral_angles_epoch : dict
                Dictionary with dihedral angles information of that epoch.
            """

            dihedral_angles_epoch = {}
            file_list = [k for k in file_list if 'trajectory' in k]

            with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
                results = [executor.submit(
                    trajectory_processer, path, file, atom_list) for file in file_list]

                for f in concurrent.futures.as_completed(results):
                    trajectory_number, dihedral_angles_trajectory = f.result()
                    dihedral_angles_epoch[int(
                        trajectory_number)] = dihedral_angles_trajectory

            return dihedral_angles_epoch

        files = os.listdir(path_output)
        numeric_files = [s for s in files if s.isnumeric()]

        simulation_dict = {}

        #
        print('     -   Retrieving trajectory data...')
        #

        start_time = time.perf_counter()

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

        simulation_df = simulation_df.astype(
            {"epoch": int, "trajectory": int, "model": int, "rotatable bond": int})
        simulation_df = simulation_df.sort_values(['epoch', 'trajectory'])

        final_time = time.perf_counter()
        print('     -   Time to retrieve data: ' + str(final_time - start_time))

        return simulation_df

    def dataframes_to_vectors(simulation_df,
                              residency_df):
        """
        Function
        ----------
        From the data frames extract the vector that we are going to use for 
        the histogram generation.

        Parameters
        ----------
        - simulation_df : pd.DataFrame
            Dataframe with psi and phi information for all rotatable bonds in the 
            protein's back bone in all the models.
        - residency_df : pd.DataFrame
            Dataframe with residency information for all models in all trajectories.

        Returns
        ----------
        - weight_vector : np.array
            Vector with same length as psi_vector and phi_vector that informs about 
            the angles' weight, which is de residency.
        - psi_vector : np.array
            Vector with the psi values reached during the simulation.
        - phi_vector : np.array
            Vector with the psi values reached during the simulation.       
        """

        len_protein = simulation_df.loc[simulation_df['rotatable bond'].idxmax(
        )]['rotatable bond'] + 1  # Length of the protein's backbone

        residency = residency_df['residency'].to_numpy()
        weight_vector = np.repeat(residency, len_protein)  # Repeat residency

        # Matrix of [len_prot*number_of_models,2]
        matrix = list(simulation_df['value'].to_numpy())

        psi_vector = np.array(matrix)[:, 0]
        phi_vector = np.array(matrix)[:, 1]

        return weight_vector, psi_vector, phi_vector

    path, path_output, path_results, path_images = path_definer(input_folder)       # Define paths

    residency_df = residency_function(path_output,                                  # Store residency of models
                                      path)

    atom_list = ['  N   ', '  CA  ', '  C   ']

    simulation_df = trajectory_positions(path_output,                               # Store dihedral angles calculated
                                         atom_list,
                                         cpus)

    weight_vector, psi_vector, phi_vector = dataframes_to_vectors(simulation_df,    # Tranforming data for the clustering
                                                                  residency_df)

    simulation_df.to_csv(os.path.join(path_results, 'data.csv'), index=False)       # Writing all the information into a csv

    return weight_vector, psi_vector, phi_vector, path_results, path_images


def clustering(weight_vector,
               psi_vector,
               phi_vector,
               path_results,
               path_images):
    """
    Function
    ----------
    Cluster the results obtained and stored in three vectors.

    Parameters
    ----------
    - weight_vector : np.array
        Vector with same length as psi_vector and phi_vector that informs about 
        the angles' weight, which is de residency.
    - psi_vector : np.array
        Vector with the psi values reached during the simulation.
    - phi_vector : np.array
        Vector with the psi values reached during the simulation.  
    - path_results : str 
        Path to the directory where the results will be stored.
    - path_images : str 
        Path to the directory where the images will be stored.

    """

    def binning(weight_vector,
                psi_vector,
                phi_vector,
                path_images):
        """
        Function
        ----------
        Cluster the results obtained from simulation with binning.
        Entropic contributions are calculated from the binned data.

        Parameters
        ----------
        - weight_vector : np.array
            Vector with same length as psi_vector and phi_vector that informs about 
            the angles' weight, which is de residency.
        - psi_vector : np.array
            Vector with the psi values reached during the simulation.
        - phi_vector : np.array
            Vector with the psi values reached during the simulation.  
        - path_images : str 
            Path to the directory where the results will be stored. 
        """

        bin_edges_psi = np.histogram_bin_edges(
            psi_vector, bins=10)                                # Entropy calculation psi
        density_psi, _ = np.histogram(
            psi_vector, bins=bin_edges_psi, density=True, weights=weight_vector)
        dense_bins_psi = density_psi[density_psi != 0]

        entropy_contribution_psi = np.sum(
            np.array([p*np.log(p) for p in dense_bins_psi]))    # Entropy psi

        bin_edges_phi = np.histogram_bin_edges(
            phi_vector, bins=10)                                # Entropy calculation phi
        density_phi, _ = np.histogram(
            phi_vector, bins=bin_edges_phi, density=True)
        dense_bins_phi = density_phi[density_phi != 0]

        entropy_contribution_phi = np.sum(
            np.array([p*np.log(p) for p in dense_bins_phi]))    # Entropy phi

        # Plot
        fig, axs = plt.subplots(1, 2, figsize=(15, 5))
        fig.suptitle('$\Psi$ and $\Phi$ dihedral distribution')

        # Sublot
        axs[0].set_title('$\Psi$')
        axs[0].hist(psi_vector,
                    bins=bin_edges_psi, density=True, weights=weight_vector)
        axs[0].set_xlabel('Dihedral angle (º)')
        axs[0].set_ylabel('Density')
        axs[0].set_xlim(-180, 180)
        axs[0].set_xticks(list(np.arange(-180, 190, 30)))

        # Sublot
        axs[1].set_title('$\Phi$')
        axs[1].hist(phi_vector,
                    bins=bin_edges_phi, density=True, weights=weight_vector)
        axs[1].set_xlabel('Dihedral angle (º)')
        axs[1].set_ylabel('Density')
        axs[1].set_xlim(-180, 180)
        axs[1].set_xticks(list(np.arange(-180, 190, 30)))

        plt.savefig(os.path.join(
            path_images, 'dihedral_distribution.png'), format='png', transparent=True)
        plt.close()

        entropy_contribution = entropy_contribution_psi + entropy_contribution_phi
        S = -(R)*entropy_contribution

        entropy_percentage_psi = 100.*entropy_contribution_psi/entropy_contribution
        entropy_percentage_phi = 100.*entropy_contribution_phi/entropy_contribution

        entropy_df = pd.DataFrame({'entropy (kcal/mol)': [S],
                                   'percentage phi': [entropy_percentage_phi],
                                   'percentage psi': [entropy_percentage_psi]}
                                   )

        entropy_df.to_csv(os.path.join(
            path_results, 'entropy.csv'), index=False)

        print(' -   Entropic information written in /dihedrals/entropy.csv.')

    print(' -   Clustering data obteined...')

    binning(weight_vector,  # Clustering results to calculate entropy and generate image.
            psi_vector,
            phi_vector,
            path_images)


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

    weight_vector, \
        psi_vector, \
        phi_vector, \
        path_results, \
        path_images = dihedral_angles_retriever_main(input_folder=args.input_folder,    # Retrieve data from the simulation
                                                     cpus=args.cpus)                  # and process it

    #
    print(' ')
    print(' -   Beginning of clustering')
    #

    start_time = time.perf_counter()

    clustering(weight_vector,                                                       # Clustering data to calculate entropic
               psi_vector,                                                          # contribution and generating image
               phi_vector,                                                          # of the simulation.
               path_results,
               path_images)

    final_time = time.perf_counter()
    print(' -   Clustering time: ', final_time - start_time)

    #
    print(' ')
    print(' -   Results written in /dihedrals.')


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])

    start_time = time.perf_counter()
    main(args)
    final_time = time.perf_counter()

    print(' -   Total time: ', final_time - start_time, '\n')
