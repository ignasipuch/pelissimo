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
from termios import NL1
import numpy as np
from itertools import groupby, chain

from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import AllChem

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
#        conf = AllChem.EmbedMolecule(m)
#        conf = m.GetConformer(0)

        heavy_atoms = [k for k in atoms if 'H' not in k]

        neighbours_dict = {}
        dihedral_bond_dict = {}

        for atom in range(len(heavy_atoms)):

            #print(atom,heavy_atoms[atom])

            neighbours_dict[heavy_atoms[atom]] = [atoms[x.GetIdx()] for x in m.GetAtomWithIdx(atom).GetNeighbors()]

        #print(rdMolTransforms.GetDihedralDeg(conf, 10,15,16,17))
        #print(rdMolTransforms.GetDihedralDeg(conf, 7,10,14,15))

        atom_list = []

        for key in rotatable_bonds_dict:

            all_neighbours_atom_1 = neighbours_dict[rotatable_bonds_dict[key][0]]
            all_neighbours_atom_2 = neighbours_dict[rotatable_bonds_dict[key][1]]

            neighbours_atom_1 = [x for x in all_neighbours_atom_1 if rotatable_bonds_dict[key][1] not in x]
            neighbours_atom_2 = [x for x in all_neighbours_atom_2 if rotatable_bonds_dict[key][0] not in x]

            dihedral_bond_dict[key] = neighbours_atom_1[0],\
                                      rotatable_bonds_dict[key][0],\
                                      rotatable_bonds_dict[key][1],\
                                      neighbours_atom_2[0]
            
            atom_list.append(neighbours_atom_1[0])
            atom_list.append(rotatable_bonds_dict[key][0])
            atom_list.append(rotatable_bonds_dict[key][1])
            atom_list.append(neighbours_atom_2[0])

        atom_list = sorted(set(atom_list))

        return dihedral_bond_dict, atom_list

    def trajectory_positions(path_output,
                             atom_list,
                             dihedral_bond_dict):
        
        def trajectory_positions_retriever(path,
                                           file_list,
                                           atom_list,
                                           dihedral_bond_dict):

            def dihedral_angle_calculator(dihedral_bond_dict,
                                          atoms_positions_dict):

                dihedral_angle_dict = {}

                for key in dihedral_bond_dict:

                    atom_name_a = dihedral_bond_dict[key][0]
                    atom_name_b = dihedral_bond_dict[key][1]
                    atom_name_c = dihedral_bond_dict[key][2]
                    atom_name_d = dihedral_bond_dict[key][3]

                    coordinates_atom_a = np.array(atoms_positions_dict[atom_name_a])
                    coordinates_atom_b = np.array(atoms_positions_dict[atom_name_b])
                    coordinates_atom_c = np.array(atoms_positions_dict[atom_name_c])
                    coordinates_atom_d = np.array(atoms_positions_dict[atom_name_d])

                    b1 = coordinates_atom_b - coordinates_atom_a
                    b2 = coordinates_atom_c - coordinates_atom_b; b2_u = b2/np.linalg.norm(b2)
                    b3 = coordinates_atom_d - coordinates_atom_c

                    n1 = np.cross(b1,b2)/np.linalg.norm(np.cross(b1,b2))
                    n2 = np.cross(b2,b3)/np.linalg.norm(np.cross(b2,b3))
                    m1 = np.cross(n1,b2_u)

                    y = np.dot(m1,n2) 
                    x = np.dot(n1,n2)
                    
                    dihedral_angle = np.arctan2(y,x)
                    dihedral_angle_dict[key] = dihedral_angle
                
                return dihedral_angle_dict

            dihedral_angles_epoch = {}
            file_list = [k for k in file_list if 'trajectory' in k]

            for file in file_list:

                if 'trajectory' in file:

                    dihedral_angles_trajectory = {}
                    trajectory_number = file.split('_')[-1].split('.pdb')[0]

                    filein = open(os.path.join(path,file), "r")
                    content = filein.read()
                    models_list = content.split('MODEL')
                    models_list = models_list[1:len(models_list)]

                    model_cont = 0

                    for model in models_list:

                        atoms_positions_dict = {}

                        model_lines = model.splitlines(True)
                        ligand_lines = [line for line in model_lines if 'LIG' in line] 

                        for line in ligand_lines:

                            if any(ext in line for ext in atom_list):

                                line = line.split()
                                atoms_positions_dict[line[2]] = [float(line[5]),
                                                                 float(line[6]),
                                                                 float(line[7])]

                        dihedral_angles_trajectory[model_cont] = dihedral_angle_calculator(dihedral_bond_dict,
                                                                    atoms_positions_dict)
                        
                        model_cont += 1

                dihedral_angles_epoch[trajectory_number] = dihedral_angles_trajectory
        
            return dihedral_angles_epoch
        
        files = os.listdir(path_output)
        numeric_files = [s for s in files if s.isnumeric()]

        if len(numeric_files) != 0:

            simulation_dict = {}

            for document in numeric_files:

                new_directory = os.path.join(path_output, document)

                if os.path.isdir(new_directory) and document.isnumeric():

                    files = os.listdir(new_directory)       

                    simulation_dict[document] = trajectory_positions_retriever(new_directory,
                                                                               files,
                                                                               atom_list,
                                                                               dihedral_bond_dict)

        else:

            simulation_dict = trajectory_positions_retriever(path_output,
                                                        files,
                                                        atom_list,
                                                        dihedral_bond_dict)
        
        return simulation_dict

    path_template, path_output = path_definer(input_folder)
    rotatable_bonds_dict = template_info_retriever(path_template,
                                                   residue_name)
    dihedral_bond_dict, atom_list = atoms_to_track(residue_name,
                                                   rotatable_bonds_dict)
    simulation_dict = trajectory_positions(path_output,
                                           atom_list,
                                           dihedral_bond_dict)

    print(simulation_dict)

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