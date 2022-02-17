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

import numpy as np
import pandas as pd
from rdkit import Chem
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

    parser.add_argument("-f", "--file", type=str, dest="input_file",
                        default=None, help="Name of the file corresponding to the isolated ligand with connectivity.")
    parser.add_argument("-d", "--directory", type=str, dest="input_folder",
                        default='output', help="Name of the output directory where the simulation\
        is located.")
    parser.add_argument("-r", "--residue_name", type=str, dest="residue_name",
                        default='LIG', help="Ligand's residue name.")
    parser.add_argument("-cm", "--clustering_method", type=str, dest="clustering_method",
                        default='bin', help="Method to cluster data: bin or kmeans.")
    parser.add_argument("-nc", "--n_clusters", type=int, dest="n_clusters",
                        default=None, help="Number of clusters to cluster the data.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def dihedral_angles_retriever_main(input_folder,
                                   residue_name,
                                   input_file):
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

        if os.path.exists(path_results) == False:
            os.mkdir(path_results)

        return path_template, path_output, path_results

    def template_info_retriever(path_template,
                                residue_name):
        """
        Function
        ----------
        Retrieve the information about the rotatable bonds of the ligand. 

        Parameters
        ----------
        - path_template : str
            Path to the rotatable bonds file.
        - residue_name : str
            Residue name of the ligand.

        Returns
        ----------
        - rotatable_bonds_dict : dict
            Dictionary with the rotatable bonds of the ligand.
        """

        rotatable_bonds_dict = {}

        with open(os.path.join(path_template, residue_name + '.rot.assign')) as filein:

            cont = 1

            for line in filein:

                line = line.split()

                if line[0].startswith('sidelib'):

                    rotatable_bonds_dict[cont] = tuple(
                        [line[2].split('_')[1], line[3].split('_')[1]])
                    cont += 1

        return rotatable_bonds_dict

    def atoms_to_track(residue_name,
                       rotatable_bonds_dict,
                       input_file):
        """
        Function
        ----------
        Selects the atoms that are going to be checked to calculate dihedrals involving
        the rotatable bonds. 

        Parameters
        ----------
        - residue_name : str
            Residue name of the ligand.
        - rotatable_bonds_dict : dict
            Dictionary with the rotatable bonds of the ligand.

        Returns
        ----------
        - dihedral_bond_dict : dict
            Dictionary with the atoms involved in the calculation of the dihedrals.
        - atom_list : list
            List with all the atoms without repetition that have to be tracked.
        """

        path = str(pathlib.Path().absolute())

        if input_file is None:
            raise Exception(
                'InputLigandStructure: Input ligand is missing and it is necessary.')

        path_ligand = os.path.join(path, input_file)

        atoms = {}
        H_neighbours = {}
        atoms_list = []
        atom_cont = 1

        with open(path_ligand) as filein:

            for line in filein:

                if residue_name in line:

                    atom = line.split()[2]
                    atoms[atom_cont] = atom
                    atoms_list.append(atom)
                    atom_cont += 1

                if 'CONECT' in line:

                    line = line.split()[1:]
                    neighbours_connect = [atoms[int(k)] for k in line]

                    for item in neighbours_connect:
                        if 'H' in item and len(neighbours_connect) == 2:
                            H_neighbours[neighbours_connect[abs(
                                neighbours_connect.index(item)-1)]] = item

                if 'ANISOU' in line:

                    raise Exception(
                        'LigandFileError: ANISOU lines detected in the pdb. This lines must be erased.')

        m = Chem.rdmolfiles.MolFromPDBFile(path_ligand)

        heavy_atoms = [k for k in atoms_list if 'H' not in k]

        neighbours_dict = {}
        dihedral_bond_dict = {}
        atom_list = []

        for atom in range(len(heavy_atoms)):

            neighbours_dict[heavy_atoms[atom]] = [
                atoms_list[x.GetIdx()] for x in m.GetAtomWithIdx(atom).GetNeighbors()]

        for key in rotatable_bonds_dict:

            all_neighbours_atom_1 = neighbours_dict[rotatable_bonds_dict[key][0]]
            all_neighbours_atom_2 = neighbours_dict[rotatable_bonds_dict[key][1]]

            neighbours_atom_1 = [
                x for x in all_neighbours_atom_1 if rotatable_bonds_dict[key][1] not in x]
            neighbours_atom_2 = [
                x for x in all_neighbours_atom_2 if rotatable_bonds_dict[key][0] not in x]

            if len(neighbours_atom_1) == 0:

                if not H_neighbours:
                    raise Exception(
                        'DihedralBondError: Information about the connectivity of the PDB will be required due to lack of heavy atoms to calculate certain dihedral bond angles.')

                neighbours_atom_1 = [
                    H_neighbours[rotatable_bonds_dict[key][0]]]

            if len(neighbours_atom_2) == 0:

                if not H_neighbours:
                    raise Exception(
                        'DihedralBondError: Information about the connectivity of the PDB will be required due to lack of heavy atoms to calculate certain dihedral bond angles.')

                neighbours_atom_2 = [
                    H_neighbours[rotatable_bonds_dict[key][1]]]

            dihedral_bond_dict[key] = neighbours_atom_1[0],\
                rotatable_bonds_dict[key][0],\
                rotatable_bonds_dict[key][1],\
                neighbours_atom_2[0]

            atom_list.append(neighbours_atom_1[0])
            atom_list.append(rotatable_bonds_dict[key][0])
            atom_list.append(rotatable_bonds_dict[key][1])
            atom_list.append(neighbours_atom_2[0])

        atom_list = sorted(set(atom_list))

        dihedral_bond_df = pd.DataFrame.from_dict(dihedral_bond_dict)
        dihedral_bond_df = dihedral_bond_df.T
        dihedral_bond_df.columns = ['A', 'B', 'C', 'D']
        index = dihedral_bond_df.index
        index.name = "Rot. bond index"

        return dihedral_bond_dict, atom_list, dihedral_bond_df

    def trajectory_positions(path_output,
                             atom_list,
                             dihedral_bond_dict):
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
                                           atom_list,
                                           dihedral_bond_dict,
                                           epoch):
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

            def dihedral_angle_calculator(dihedral_bond_dict,
                                          atoms_positions_dict):
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

                dihedral_angle_dict = {}

                for key in dihedral_bond_dict:

                    atom_name_a = dihedral_bond_dict[key][0]
                    atom_name_b = dihedral_bond_dict[key][1]
                    atom_name_c = dihedral_bond_dict[key][2]
                    atom_name_d = dihedral_bond_dict[key][3]

                    coordinates_atom_a = np.array(
                        atoms_positions_dict[atom_name_a])
                    coordinates_atom_b = np.array(
                        atoms_positions_dict[atom_name_b])
                    coordinates_atom_c = np.array(
                        atoms_positions_dict[atom_name_c])
                    coordinates_atom_d = np.array(
                        atoms_positions_dict[atom_name_d])

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
                    dihedral_angle_dict[key] = dihedral_angle

                return dihedral_angle_dict

            dihedral_angles_epoch = {}
            file_list = [k for k in file_list if 'trajectory' in k]

            for file in file_list:

                if 'trajectory' in file:

                    dihedral_angles_trajectory = {}
                    trajectory_number = file.split('_')[-1].split('.pdb')[0]

                    filein = open(os.path.join(path, file), "r")
                    content = filein.read()
                    models_list = content.split('MODEL')
                    models_list = models_list[1:len(models_list)]

                    model_cont = 1

                    for model in models_list:

                        atoms_positions_dict = {}

                        model_lines = model.splitlines(True)
                        ligand_lines = [
                            line for line in model_lines if 'LIG' in line]

                        for line in ligand_lines:

                            if any(ext in line for ext in atom_list):

                                line = line.split()

                                if len(line[4]) == 5:
                                    atoms_positions_dict[line[2]] = [
                                        float(line[5]), float(line[6]), float(line[7])]

                                else:
                                    atoms_positions_dict[line[2]] = [
                                        float(line[6]), float(line[7]), float(line[8])]

                        dihedral_angles_trajectory[model_cont] =\
                            dihedral_angle_calculator(
                                dihedral_bond_dict, atoms_positions_dict)

                        model_cont += 1

                dihedral_angles_epoch[int(
                    trajectory_number)] = dihedral_angles_trajectory

            return dihedral_angles_epoch

        files = os.listdir(path_output)
        numeric_files = [s for s in files if s.isnumeric()]

        simulation_dict = {}

        if len(numeric_files) != 0:

            for epoch in numeric_files:

                new_directory = os.path.join(path_output, epoch)

                if os.path.isdir(new_directory) and epoch.isnumeric():

                    files = os.listdir(new_directory)

                    simulation_dict[int(epoch)] = trajectory_positions_retriever(new_directory,
                                                                                 files,
                                                                                 atom_list,
                                                                                 dihedral_bond_dict,
                                                                                 epoch)

        else:

            epoch = 0
            simulation_dict[0] = trajectory_positions_retriever(path_output,
                                                                files,
                                                                atom_list,
                                                                dihedral_bond_dict,
                                                                epoch)

        simulation_df = pd.DataFrame([(epoch, trajectory, model, rot_bond, value)
                                      for epoch, traj_mod_rot_val in simulation_dict.items()
                                      for trajectory, mod_rot_val in traj_mod_rot_val.items()
                                      for model, rot_val in mod_rot_val.items()
                                      for rot_bond, value in rot_val.items()])

        simulation_df.columns = ['epoch', 'trajectory',
                                 'model', 'rotable bond', 'value']

        data_cluster = simulation_df.pivot(index=['epoch', 'trajectory', 'model'],
                                           columns='rotable bond',
                                           values='value')

        return data_cluster

    path_template, path_output, path_results = path_definer(input_folder)

    #
    print('     -   Retrieving information about rotatable bonds.')
    #

    rotatable_bonds_dict = template_info_retriever(path_template,
                                                   residue_name)

    dihedral_bond_dict, atom_list, dihedral_bond_df = atoms_to_track(residue_name,
                                                                     rotatable_bonds_dict,
                                                                     input_file)

    #
    print('     -   Calculating dihedral angles of all the conformations...')
    #

    simulation_df = trajectory_positions(path_output,
                                         atom_list,
                                         dihedral_bond_dict)

    return dihedral_bond_df, simulation_df, path_results


def clustering(n_cluster,
               clustering_method,
               simulation_df,
               dihedral_bond_df,
               path_results):
    """
    Function
    ----------
    Cluster the results obtained and stored in a data frame.

    Parameters
    ----------
    - n_cluster : int
        Number of clusters t ocluster data.
    - dihedral_bond_df : pd.DataFrame
        Data frame with rotatable bonds, atoms conforming it and the index assigned.
    - simulation_df : pd.DataFrame
        Data frame with all the rotatable bonds' dihedral angle values of all the simulation
        with corresponding model, trajectory and epoch.
    - path_results : str 
        Path to the directory where the results will be stored. 

    """

    def kmeans(simulation_df,
               dihedral_bond_df,
               path_results,
               n_cluster):

        def scaler(simulation_df,
                   dihedral_bond_df):
            """
            Function
            ----------
            Scale the data from the data frame.

            Parameters
            ----------
            - simulation_df : pd.DataFrame
                Data frame with all the rotatable bonds' dihedral angle values of all the simulation
                with corresponding model, trajectory and epoch.
            - dihedral_bond_df : pd.DataFrame
                Data frame with rotatable bonds, atoms conforming it and the index assigned.

            Returns
            ----------
            - simulation_df : pd.DataFrame
                Data frame with all the rotatable bonds' dihedral angle values of all the simulation
                scaled from 0 to 1.
            """

            from sklearn.preprocessing import MinMaxScaler

            column_list = list(np.arange(1, len(dihedral_bond_df) + 1))
            scaler = MinMaxScaler()
            scaler.fit_transform(simulation_df[column_list])
            simulation_df[column_list] = scaler.transform(
                simulation_df[column_list])

            return simulation_df

        def elbow_method(path_results,
                         simulation_df,
                         dihedral_bond_df):
            """
            Function
            ----------
            Perform the elbow method to determine the optimal number of clusters to cluster the 
            data.

            Parameters
            ----------
            - path_results : str 
                Path to the directory where the results will be stored.
            - simulation_df : pd.DataFrame
                Data frame with all the rotatable bonds' dihedral angle values of all the simulation
                with corresponding model, trajectory and epoch.

            Returns
            ----------
            - n_cluster : int
                Optimal number of clusters to cluster the data.
            """

            def plotter(path_results,
                        wcss,
                        derivative,
                        second_derivative,
                        n_clusters,
                        n_cluster):
                """
                Function
                ----------
                Plost the data obtained from the elbow method.

                Parameters
                ----------
                - path_results : str 
                    Path to the directory where the results will be stored.
                - wcss : list
                    List of Within-Cluster Sum of Square values for all the different cluster numbers
                    tried. 
                - derivative : list
                    Numerical derivative of the WCSS values.
                - second_derivative : list
                    Numerical second derivative of the WCSS values.
                - n_clusters : list
                    List with all the number of clusters tried for the elbow method.
                - _n_cluster : int
                    Optimal number of clusters to cluster the data.
                """

                fig1, ax1 = plt.subplots()
                fig2, ax2 = plt.subplots()
                fig3, ax3 = plt.subplots()

                ax1.set_title('WCSS')
                ax1.plot(n_clusters, wcss)
                ax1.scatter(n_cluster, wcss[n_cluster - 2],
                            color='red', marker='x', label='elbow')
                ax1.legend(loc='best')
                ax1.set_xlabel('Number of clusters')
                ax1.set_ylabel('WCSS')
                fig1.savefig(os.path.join(path_results, 'wcss.png'))

                ax2.set_title('WCSS derivative')
                ax2.plot(n_clusters[1:-1], derivative)
                ax2.set_xlabel('Number of clusters')
                ax2.set_ylabel('WCSS derivative')
                fig2.savefig(os.path.join(path_results, 'wcss_derivative.png'))

                ax3.set_title('WCSS 2nd derivative')
                ax3.plot(n_clusters[2:-2], second_derivative)
                ax3.set_xlabel('Number of clusters')
                ax3.set_ylabel('WCSS 2nd derivative')
                fig3.savefig(os.path.join(
                    path_results, 'wcss_2derivative.png'))

            from scipy.signal import argrelextrema

            wcss = []
            n_clusters = []
            column_list = list(np.arange(1, len(dihedral_bond_df) + 1))

            for n_cluster in range(2, 17):

                km = KMeans(init='k-means++',
                            n_clusters=n_cluster, max_iter=10000)
                km.fit(simulation_df[column_list])

                wcss.append(km.inertia_)
                n_clusters.append(n_cluster)

            derivative = list(np.gradient(wcss, 1))[1:-1]
            second_derivative = np.gradient(derivative, 1)[1:-1]

            maxs_location = list(argrelextrema(
                second_derivative, np.greater)[0])

            try:

                second_derivative_max = max(
                    [second_derivative[max_id] for max_id in maxs_location])
                max_location_array = np.where(
                    second_derivative == second_derivative_max)[0] + 2
                max_location = max_location_array[0]

            except ValueError:

                max_location_array = np.where(
                    second_derivative == max(second_derivative))[0] + 2
                max_location = max_location_array[0]

            n_cluster = n_clusters[max_location]

            plotter(path_results,
                    wcss,
                    derivative,
                    second_derivative,
                    n_clusters,
                    n_cluster)

            return n_cluster

        from sklearn.cluster import KMeans

        simulation_df_copy = simulation_df.copy()

        simulation_df = scaler(simulation_df,
                               dihedral_bond_df)

        if n_cluster == None:

            #
            print('     -   No information about the number of clusters was given.')
            print('     -   Determining the optimal number of clusters.')
            #

            n_cluster = elbow_method(path_results,
                                     simulation_df,
                                     dihedral_bond_df)

            #
            print('     -   Optimal number of clusters =', n_cluster)
            print('     -   Plots have been generated succesfully.')

        else:

            #
            print('     -   Number of clusters =', n_cluster)

        #
        print('     -   Clustering...')
        #

        km = KMeans(init='k-means++', n_clusters=n_cluster, max_iter=10000)
        y_predicted = km.fit_predict(simulation_df)
        simulation_df_copy['cluster'] = y_predicted
        simulation_df_copy.to_csv(os.path.join(path_results, 'data.csv'))
        dihedral_bond_df.to_csv(os.path.join(path_results, 'dihedrals.csv'))

    def binning(simulation_df,
                path_results):

        entropy_contribution = []

        for key in simulation_df:

            bin_edges = np.histogram_bin_edges(
                simulation_df[key].to_numpy(), bins='fd')
            density, _ = np.histogram(
                simulation_df[key].to_numpy(), bins=bin_edges, density=True)
            dense_bins = density[density != 0]

            entropy_contribution.append(
                np.sum(np.array([p*np.log(p) for p in dense_bins])))

            # Plot
            plt.title('Dihedral ' + str(key) + ' distribution')
            plt.hist(simulation_df[key].to_numpy(),
                     bins=bin_edges, density=True)
            plt.xlabel('Dihedral angle (º)')
            plt.ylabel('Density')
            plt.xlim(-180, 180)
            plt.xticks(list(np.arange(-180, 190, 30)))
            plt.savefig(os.path.join(path_results, 'dihedral_' +
                                     str(key) + '_strain.png'), format='png')
            plt.close()

        entropy_contributions = np.array(entropy_contribution)
        S = -(R/simulation_df.shape[1])*np.sum(entropy_contributions)

        entropy_percentages = 100. * \
            np.array(entropy_contribution)/np.sum(entropy_contributions)
        entropy_df = pd.DataFrame(entropy_percentages)
        entropy_df = entropy_df.round(decimals=2).T
        entropy_df.columns = [
            'Dihedral_' + str(value) + '_%' for value in list(np.arange(1, simulation_df.shape[1] + 1))]
        entropy_df.insert(loc=0, column='S (kcal/mol K)',
                          value="{:3.6f}".format(S))
        entropy_df.to_csv(os.path.join(
            path_results, 'entropy.csv'), index=False)

        print(' -   Entropic information written in /dihedrals/entropy.csv ')

    if clustering_method == 'kmeans':

        kmeans(simulation_df,
               dihedral_bond_df,
               path_results,
               n_cluster)

    elif clustering_method == 'bin':

        binning(simulation_df,
                path_results)

        simulation_df.to_csv(os.path.join(path_results, 'data.csv'))
        dihedral_bond_df.to_csv(os.path.join(path_results, 'dihedrals.csv'))


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
    print('*******************************************************************')
    print('*                     peleDihedralClustering                      *')
    print('*******************************************************************')
    print(' ')
    print(' -   Gathering information')
    #

    dihedral_bond_df, simulation_df, path_results = dihedral_angles_retriever_main(input_folder=args.input_folder,
                                                                                   residue_name=args.residue_name,
                                                                                   input_file=args.input_file)

    #
    print(' ')
    print(' -   Beginning of clustering')
    #

    clustering(args.n_clusters,
               args.clustering_method,
               simulation_df,
               dihedral_bond_df,
               path_results)

    #
    print(' ')
    print(' -   Results written in /dihedrals')
    print(' ')


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
