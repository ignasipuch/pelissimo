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

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

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

    parsed_args = parser.parse_args(args)

    return parsed_args

def main_function(input_folder,
                  residue_name):

    def path_definer(input_folder):

        path = str(pathlib.Path().absolute())
        path_template = os.path.join(path,'DataLocal','LigandRotamerLibs')
        path_output = os.path.join(path,input_folder)

        return path_template, path_output
      
    def template_info_retriever(path_template,
                                residue_name):

        rotatable_bonds_dict = {}

        with open(os.path.join(path_template,residue_name + '.rot.assign')) as filein:

            cont = 1

            for line in filein:

                line = line.split()

                if line[0].startswith('sidelib'):
                    
                    rotatable_bonds_dict[cont] = tuple([line[2].split('_')[1],line[3].split('_')[1]])
                    cont += 1

        return rotatable_bonds_dict

    def atoms_to_track(residue_name,
                       rotatable_bonds_dict):

        path = str(pathlib.Path().absolute())
        path_input = os.path.join(path,'input')
        path_ligand = os.path.join(path_input,'ligand.pdb')

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

        for atom in range(len(heavy_atoms)):

            neighbours_dict[heavy_atoms[atom]] = [atoms[x.GetIdx()] for x in m.GetAtomWithIdx(atom).GetNeighbors()]

        for key in rotatable_bonds_dict:

            all_neighbours_atom_1 = neighbours_dict[rotatable_bonds_dict[key][0]]
            all_neighbours_atom_2 = neighbours_dict[rotatable_bonds_dict[key][1]]

            neighbours_atom_1 = [x for x in all_neighbours_atom_1 if rotatable_bonds_dict[key][1] not in x]
            neighbours_atom_2 = [x for x in all_neighbours_atom_2 if rotatable_bonds_dict[key][0] not in x]

            dihedral_bond_dict[key] = neighbours_atom_1[0],\
                                      rotatable_bonds_dict[key][0],\
                                      rotatable_bonds_dict[key][1],\
                                      neighbours_atom_2[0]

        return dihedral_bond_dict

    def trajectory_positions(path_output):
        
        def trajectory_positions_retriever(path,
                                           file_list):

            trajectory_dict = {}
            epoch_dict = {}

            file_list = [k for k in file_list if 'trajectory' in k]

            for file in file_list:

                if 'trajectory' in file:

                    trajectory_number = file.split('_')[-1].split('.pdb')[0]
                    atom_position = {}

                    with open(os.path.join(path,file)) as filein:

                        for line in filein:

                            if 'MODEL' in line:

                                line = line.split()
                                model = line[-1]

                            elif residue_name in line:

                                line = line.split()
                                atom_position[line[2]] = tuple([line[5],line[6],line[7]])

                            elif 'ENDMDL' in line:

                                trajectory_dict[model] = atom_position
                    
                    epoch_dict[trajectory_number] =  trajectory_dict

                
            return epoch_dict
        
        files = os.listdir(path_output)
        numeric_files = [s for s in files if s.isnumeric()]

        if len(numeric_files) != 0:

            simulation_dict = {}

            for document in numeric_files:

                new_directory = os.path.join(path_output, document)

                if os.path.isdir(new_directory) and document.isnumeric():

                    files = os.listdir(new_directory)       

                    simulation_dict[document] = trajectory_positions_retriever(new_directory,
                                                                          files)

        else:

            simulation_dict = trajectory_positions_retriever(path_output,
                                                        files)
        
        return simulation_dict

    path_template, path_output = path_definer(input_folder)
    rotatable_bonds_dict = template_info_retriever(path_template,
                                                   residue_name)
    dihedral_bond_dict = atoms_to_track(residue_name,
                                        rotatable_bonds_dict)
    simulation_dict = trajectory_positions(path_output)

    print('dihedral angles dict')
    print(rotatable_bonds_dict)
    print('dihedral bonds')
    print(dihedral_bond_dict)

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

    main_function(input_folder=args.input_folder,
                  residue_name=args.residue_name)

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)