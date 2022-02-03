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

    parsed_args = parser.parse_args(args)

    return parsed_args


def dihedral_angles_retriever_main(input_folder,
                                   residue_name):
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
        Data frame with all the rotatable bonds' dihedral angle values of all the simulation with corresponding model, trajectory and epoch.
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
        """

        path = str(pathlib.Path().absolute())
        path_template = os.path.join(path, 'DataLocal', 'LigandRotamerLibs')
        path_output = os.path.join(path, input_folder)

        return path_template, path_output

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
                       rotatable_bonds_dict):
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
        path_input = os.path.join(path, 'input')
        path_ligand = os.path.join(path_input, 'ligand.pdb')

        atoms = []

        with open(path_ligand) as filein:

            for line in filein:

                if residue_name in line:

                    atom = line.split()[2]
                    atoms.append(atom)

        m = Chem.rdmolfiles.MolFromPDBFile(path_ligand)

        heavy_atoms = [k for k in atoms if 'H' not in k]

        neighbours_dict = {}
        dihedral_bond_dict = {}
        atom_list = []

        for atom in range(len(heavy_atoms)):

            neighbours_dict[heavy_atoms[atom]] = [
                atoms[x.GetIdx()] for x in m.GetAtomWithIdx(atom).GetNeighbors()]

        for key in rotatable_bonds_dict:

            all_neighbours_atom_1 = neighbours_dict[rotatable_bonds_dict[key][0]]
            all_neighbours_atom_2 = neighbours_dict[rotatable_bonds_dict[key][1]]

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

                    model_cont = 0

                    for model in models_list:

                        atoms_positions_dict = {}

                        model_lines = model.splitlines(True)
                        ligand_lines = [
                            line for line in model_lines if 'LIG' in line]

                        for line in ligand_lines:

                            if any(ext in line for ext in atom_list):

                                line = line.split()
                                atoms_positions_dict[line[2]] = [float(line[5]),
                                                                 float(
                                                                     line[6]),
                                                                 float(line[7])]

                        dihedral_angles_trajectory[model_cont] = dihedral_angle_calculator(dihedral_bond_dict,
                                                                                           atoms_positions_dict)

                        model_cont += 1

                dihedral_angles_epoch[int(
                    trajectory_number)] = dihedral_angles_trajectory

            return dihedral_angles_epoch

        files = os.listdir(path_output)
        numeric_files = [s for s in files if s.isnumeric()]

        simulation_dict = {}

        if len(numeric_files) != 0:

            for document in numeric_files:

                new_directory = os.path.join(path_output, document)

                if os.path.isdir(new_directory) and document.isnumeric():

                    files = os.listdir(new_directory)

                    simulation_dict[int(document)] = trajectory_positions_retriever(new_directory,
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
                                 'model', 'rotable bond', 'value']

        return simulation_df

    path_template, path_output = path_definer(input_folder)
    rotatable_bonds_dict = template_info_retriever(path_template,
                                                   residue_name)
    dihedral_bond_dict, atom_list, dihedral_bond_df = atoms_to_track(residue_name,
                                                                     rotatable_bonds_dict)
    simulation_df = trajectory_positions(path_output,
                                         atom_list,
                                         dihedral_bond_dict)

    return dihedral_bond_df, simulation_df


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

    dihedral_angles_retriever_main(input_folder=args.input_folder,
                                   residue_name=args.residue_name)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
