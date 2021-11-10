# -*- coding: utf-8 -*-

# Imports
import sys
import os
import pathlib 
import argparse
import shutil
import time
from distutils.dir_util import copy_tree

# *********************************************************************** #
#                        PREPARING READING OF FILES                       #
# *********************************************************************** #

# Time
start_time = time.time()

# Parse command line
parser = argparse.ArgumentParser()

# Inputs
parser.add_argument("-d", "--directory", type=str, dest = "input_folder",\
    default = 'LIG_Pele', help="Name of the directory where the simulation\
    is located.")
parser.add_argument("-r", "--residue_name", type=str, dest = "residue_name",\
    default = 'LIG', help="Ligand's residue name.")
parser.add_argument("-cl", "--clusters_folder", type=str, dest = "clusters_folder",\
    default = 'results', help="Name of the directory containing the folder: clusters.")
parser.add_argument("-co", "--conf_file_name", type=str, dest = "conf_file_name",\
    default = 'pele', help="Name of the .conf file used for the simulation.")

args = parser.parse_args()

input_folder = args.input_folder
residue_name = args.residue_name
clusters_folder = args.clusters_folder
conf_file_name = args.conf_file_name

# Creating a variable with the paths that we are going to use
path = str(pathlib.Path().absolute())
path_simulation = path + '/' + input_folder
path_results = path_simulation + '/' + clusters_folder
path_clusters = path_results + '/clusters' 

# Checking existance
if os.path.isdir(path_simulation) == False:
    raise Exception('PathSimulationError: There is no folder with this name. Please check the path and the folder name.')

# Listing directories inside the directory
files = os.listdir(path_simulation)

# Creating a folder to store the info for the calculations
path_energies = path + '/' + residue_name + '_peleRes'
path_energies_input = path_energies + '/' + 'input'
path_energies_simulation = path_energies + '/' + 'simulation'

if  os.path.exists(path_energies) == False:
    os.mkdir(path_energies)

if  os.path.exists(path_energies_input) == False:
    os.mkdir(path_energies_input)

if  os.path.exists(path_energies_simulation) == False:
    os.mkdir(path_energies_simulation)

# *********************************************************************** #
#                    COPYING FILES TO NEW DIRECTORY                       #
# *********************************************************************** #

# Copying the clusters as input in the new folder and saving their names
cluster_files = []
labels = []

# Printing header
print(' ')
print('*******************************************************************')
print('*                        peleLigandResults                        *')
print('* --------------------------------------------------------------- *')
print('*      Ligand\'s internal energy from induced fit results          *')
print('*******************************************************************')
print(' ')

# Storing information and copying files
if os.path.isdir(path_clusters): 

    files = os.listdir(path_clusters)

    for document in files: 

        if 'cluster' in document and '.pdb' in document:

            cluster_files.append(os.path.join(path_energies_simulation,document))
            labels.append((document.split('cluster_'))[1].split('.pdb')[0])
            shutil.copy(os.path.join(path_clusters,document), path_energies_input)

# List of clusters' letters for the run_file
run_file_labels = ' '.join(labels)

#
print(' -   Number of clusters obtained in the simulation:',len(labels))
#

# Removing the protein from the pdbs.
clusters = os.listdir(path_energies_input)

for cluster in clusters:

    with open(os.path.join(path_energies_input,cluster)) as filein:

        lines = (l for l in filein if residue_name in l)
        new_path = path_energies_simulation + '/' + cluster.split('.pdb')[0]
        path_DataLocal = path_energies_simulation + '/DataLocal'

        if  os.path.exists(new_path) == False:
            os.mkdir(new_path)

        if  os.path.exists(path_DataLocal) == False:
            os.mkdir(path_DataLocal)

        copy_tree(os.path.join(path_simulation,'DataLocal'), path_DataLocal)
           
        with open(os.path.join(new_path,cluster), 'w') as fileout:

            fileout.writelines(lines)


# Copying .conf document and extracting information
files =  files = os.listdir(path_simulation)

#
print(' -   Extracting information from ' + conf_file_name + '.conf.' )
#

for document in files:

    if conf_file_name + '.conf' in document:

        with open (os.path.join(path_simulation,document)) as filein:

            for line in filein:

                if "ForceField" in line:

                    line = line.split(':')
                    line = line[1].split('     ')
                    line = line[0].split('"')
                    forcefield = line[1]

                elif "VDGBNP" in line:

                    solvent = 'VDGBNP'

                elif "OBC" in line:

                    solvent = 'OBC'

#
print('     -   Forcefield used:', forcefield + '.' )
print('     -   Solvent model used:', solvent + '.' )
#

# Generating all the necessary control files
cont = 0

#
print(' -   Generating control files for the energy calculation.' )
#

for label in labels:

    new_path = os.path.join(path_energies_simulation + '/cluster_' + label)

    with open (os.path.join(new_path,'energy' + label + '.conf'), 'w') as fileout:

        fileout.writelines(
        '{\n'
        '   "licenseDirectoryPath" : "/gpfs/projects/bsc72/PELE++/license",'
        '\n'
        '   "Initialization" : {\n'
        '      "Complex" : {\n'
        '         "files" : [\n'
        '            {\n'
        '               "path": "' + new_path + '/cluster_' + label + '.pdb"\n'
        '            }\n'
        '         ]\n'
        '      },\n'
        '      "ForceField" : "' + forcefield + '",\n'
        '      "Solvent" : {\n'
        '         "ionicStrength" : 0.250,\n'
        '         "solventType" : "' + solvent + '",\n'
        '         "useDebyeLength" : true\n'
        '      }\n'
        '   },\n'
        '   "commands" : [\n'
        '      {\n'
        '       "commandType":"energyComputation"\n'
        '      }\n'
        '   ]\n'
        '}\n'
        )
        
    cont += 1

with open (os.path.join(path_energies_simulation,'run'), 'w') as fileout:

    fileout.writelines(
    '#!/bin/bash\n'
    '#SBATCH -J PELEne\n'
    '#SBATCH --output=PELEne.out\n'
    '#SBATCH --error=PELEne.err\n'
    '#SBATCH --qos=debug\n'
    '#SBATCH --time=00:30:00\n'
    '\n'
    'module purge\n'
    'module load intel mkl impi gcc\n'
    'module load impi\n'
    'module load boost/1.64.0\n'
    '\n'
    'list="' + run_file_labels + '"\n'
    '\n'
    'for i in $list\n'
    'do\n'
    '\n'
    '    echo " --------------------------------------------------------------------"\n'
    '    echo "|                            CLUSTER $i                              |"\n'
    '    echo " --------------------------------------------------------------------"\n'
    '    /gpfs/projects/bsc72/PELE++/mniv/V1.7.1/bin/PELE-1.7.1_serial /gpfs/projects/bsc72/ignasi/PhD/strain/second_set/MTAP/OPLS/normal/1CB0/LIG_peleRes/simulation/cluster_${i}/energy${i}.conf\n'
    '    echo " "\n'
    '    echo "**********************************************************************"\n'
    '    echo "**********************************************************************"\n'
    '    echo " "\n'
    '\n'
    'done\n'
    )



#
print(' ')
print('                    --Duration of the execution--                   ')
print('                      %s seconds' % (time.time() - start_time)  )
print(' ')
print('*******************************************************************')
print(' -   To run the energy calculation for all the clusters:')
print('     :> cd ' + residue_name + '_peleRes/simulation')
print('     :> sbatch run')
print(' -   Results are stored in PELEne.out.')
#

                


                

