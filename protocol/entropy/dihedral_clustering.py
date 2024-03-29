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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
                        default=None, required=True, help="Name of the file corresponding to the isolated ligand with connectivity.")
    parser.add_argument("-d", "--directory", type=str, dest="input_folder",
                        default='output', help="Name of the output directory where the simulation\
        is located.")
    parser.add_argument("-r", "--residue_name", type=str, dest="residue_name",
                        default='LIG', help="Ligand's residue name.")
    parser.add_argument("--evolution", dest="evolution_bool",
                        default=False, action='store_true', help="Flag to choose if dihedral evolution is wanted.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def dihedral_angles_retriever_main(input_folder,
                                   residue_name,
                                   input_file,
                                   evolution_bool):
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

        residency_df.pivot(index=['epoch', 'trajectory', 'model'],
                            columns='residency')

        return residency_df        

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

        if len(rotatable_bonds_dict) == 0:
            raise Exception('NoRotatableBonds: No rotatable bonds were found in ' + 
            str(os.path.join(path_template, residue_name + '.rot.assign')))

        return rotatable_bonds_dict

    def atoms_to_track(rotatable_bonds_dict,
                       input_file):
        """
        Function
        ----------
        Selects the atoms that are going to be tracked to calculate dihedrals involving
        the rotatable bonds. 

        Parameters
        ----------
        - rotatable_bonds_dict : dict
            Dictionary with the rotatable bonds of the ligand.
        - input_file : str
            Name of the file corresponding to the isolated ligand with connectivity.

        Returns
        ----------
        - dihedral_bond_dict : dict
            Dictionary with the atoms involved in the calculation of the dihedrals.
        - atom_list : list
            List with all the atoms without repetition that have to be tracked.
        - dihedral_bond_df : pd.DataFrame
            Data frame with the rotatable bonds' information.
        """

        path = str(pathlib.Path().absolute())

        if input_file is None:
            raise Exception(
                'InputLigandStructure: Input ligand is missing and it is necessary.')

        path_ligand = os.path.join(path, input_file)

        atoms = {}
        neighbours = {}
        atoms_list = []
        atom_cont = 1

        with open(path_ligand) as filein:

            for line in filein:

                if 'HETATM' in line:

                    atom = line.split()[2]
                    atoms[atom_cont] = atom
                    atoms_list.append(atom)
                    atom_cont += 1

                if 'CONECT' in line:

                    line = line.split()[1:]
                    neighbours_connect = [atoms[int(k)] for k in line]

                    try:
                        
                        if neighbours[neighbours_connect[0]] is not None:  #Double bonds
                            pass

                    except KeyError:

                        neighbours[neighbours_connect[0]] = neighbours_connect[1:]
                    
                if 'ANISOU' in line:

                    raise Exception(
                        'LigandFileError: ANISOU lines detected in the pdb. This lines must be erased.')

        dihedral_bond_dict = {}
        atom_list = []

        for key in rotatable_bonds_dict:

            all_neighbours_atom_1 = neighbours[rotatable_bonds_dict[key][0]]
            all_neighbours_atom_2 = neighbours[rotatable_bonds_dict[key][1]]

            neighbours_atom_1 = [
                x for x in all_neighbours_atom_1 if rotatable_bonds_dict[key][1] not in x]
            neighbours_atom_2 = [
                x for x in all_neighbours_atom_2 if rotatable_bonds_dict[key][0] not in x]

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
                                           dihedral_bond_dict):
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
                            line for line in model_lines if 'HETATM' in line]

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
                                                                                 dihedral_bond_dict)

        else:

            simulation_dict[0] = trajectory_positions_retriever(path_output,
                                                                files,
                                                                atom_list,
                                                                dihedral_bond_dict)

        simulation_df = pd.DataFrame([(epoch, trajectory, model, rot_bond, value)
                                      for epoch, traj_mod_rot_val in simulation_dict.items()
                                      for trajectory, mod_rot_val in traj_mod_rot_val.items()
                                      for model, rot_val in mod_rot_val.items()
                                      for rot_bond, value in rot_val.items()])

        simulation_df.columns = ['epoch', 'trajectory',
                                 'model', 'rotatable bond', 'value']

        simulation_df.pivot(index=['epoch', 'trajectory', 'model'],
                            columns='rotatable bond',
                            values='value')

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

        for _, row in residency_df.iterrows():
            epoch = row['epoch']
            trajectory = row['trajectory']
            model = row['model']
            residency = row['residency']

            epoch_df = simulation_df.loc[(simulation_df['epoch'] == epoch)]
            trajectory_df = epoch_df.loc[(epoch_df['trajectory'] == trajectory)]
            row_df = trajectory_df.loc[(trajectory_df['model'] == model)]
            
            if residency - 1 >= 1:  
                for i in range(residency-1):
                    simulation_df = simulation_df.append(row_df,ignore_index=True)

        return simulation_df

    def dihedrals_evolution(simulation_df,
                            path_images):
        """
        Function
        ----------
        Plot the evolution of the dihedral angles throughout the simulation.

        Parameters
        ----------
        - simulation_df : pd.DataFrame
            Data frame with all the rotatable bonds' dihedral angle values of the 
            all the simulation with corresponding model, trajectory and epoch.
        - path_images : str
            Path to the directory where images are stored.

        """

        number_of_trajectories = simulation_df['trajectory'].max()
        number_of_dihedrals = simulation_df['rotatable bond'].max()

        for dihedral in range(1,number_of_dihedrals + 1):

            dihedral_df = simulation_df[simulation_df['rotatable bond'] == dihedral]
            #fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
            angles = []
            distance = []

            for trajectory in range(1, number_of_trajectories + 1):


                # Maybe it is not necessary to extract values to an array to then put them again in
                # a dataframe (?)
                dihedral_trajectory_df = dihedral_df[dihedral_df['trajectory'] == trajectory].reset_index(drop=True)

                r = dihedral_trajectory_df.index.to_numpy()
                distance.extend(list(r))  

                theta = dihedral_trajectory_df['value'].to_numpy() 
                angles.extend(list(theta))       

            df = pd.DataFrame({
                'radii': r,
                'angles': theta
                })

            g = sns.JointGrid(data=df, x='radii',y='angles')
            x, y = df.radii, df.angles

            sns.scatterplot(x=x, y=y, ax=g.ax_joint)
            sns.kdeplot(y=y, ax = g.ax_marg_y)
            g.fig.suptitle('Dihedral ' + str(dihedral) + ' evolution')
            g.set_axis_labels(xlabel='Pele Steps', ylabel='Angles (º)')
            g.ax_marg_x.set_xlim(0,1000)
            g.ax_marg_y.set_ylim(-180,180)
            plt.yticks(list(np.arange(-180, 190, 30)))
            g.savefig(os.path.join(path_images, 'dihedral_' +
                                    str(dihedral) + '_evolution.png'), format='png', transparent=True)
            plt.close()


    path, path_template, path_output, path_results, path_images = path_definer(input_folder)

    #
    print('     -   Retrieving information about rotatable bonds.')
    #

    rotatable_bonds_dict = template_info_retriever(path_template,
                                                   residue_name)

    dihedral_bond_dict, atom_list, dihedral_bond_df = atoms_to_track(rotatable_bonds_dict,
                                                                     input_file)

    #
    print('     -   Calculating dihedral angles of all the conformations...')
    #

    residency_df = residency_function(path_output,
                                      path)

    simulation_df = trajectory_positions(path_output,
                                         atom_list,
                                         dihedral_bond_dict)
    
    full_df = residency_to_simulation(residency_df,
                                      simulation_df)

    if evolution_bool: 

        #
        print('     -   Calculating dihedral\'s evolution...')
        #

        dihedrals_evolution(full_df,
                            path_images)


    return dihedral_bond_df, simulation_df, path_results, path_images


def clustering(simulation_df,
               dihedral_bond_df,
               path_results,
               path_images):
    """
    Function
    ----------
    Cluster the results obtained and stored in a data frame.

    Parameters
    ----------
    - simulation_df : pd.DataFrame
        Data frame with all the rotatable bonds' dihedral angle values of all the simulation
        with corresponding model, trajectory and epoch.
    - dihedral_bond_df : pd.DataFrame
        Data frame with rotatable bonds, atoms conforming it and the index assigned.
    - path_results : str 
        Path to the directory where the results will be stored. 

    """

    def binning(simulation_df,
                path_images):
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

        entropy_contribution = []

        rotatable_bonds = simulation_df['rotatable bond'].to_numpy()
        values = simulation_df['value'].to_numpy()

        results = defaultdict(list)

        for rot_bond, value in zip(rotatable_bonds, values):
            results[rot_bond].append(value)

        rot_bond_values = list(results.items())
     
        for rot_bond, values in rot_bond_values:

            bin_edges = np.histogram_bin_edges(values, bins=15)
            density, _ = np.histogram(
                values, bins=bin_edges, density=True)
            dense_bins = density[density != 0]

            entropy_contribution.append(
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
                                     str(rot_bond) + '_distribution.png'), format='png', transparent=True)
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

    binning(simulation_df,
                path_images)

    simulation_df.to_csv(os.path.join(path_results, 'data.csv'), index=False)
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
    print('* Dihedral Clustering *')
    print(' ')
    print(' -   Gathering information')
    #

    dihedral_bond_df, \
    simulation_df,    \
    path_results,     \
    path_images = dihedral_angles_retriever_main(input_folder=args.input_folder,
                                                 residue_name=args.residue_name,
                                                 input_file=args.input_file,
                                                 evolution_bool=args.evolution_bool)

    #
    print(' ')
    print(' -   Beginning of clustering')
    #

    clustering(simulation_df,
               dihedral_bond_df,
               path_results,
               path_images)

    #
    print(' ')
    print(' -   Results written in /dihedrals.\n')


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)

