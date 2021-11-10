# -*- coding: utf-8 -*-

# Imports
import sys
import os
import pathlib 
import pandas as pd
import re
from collections import Counter
import argparse
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------- #
# Constants:
T = 298.
R = 1.985e-3
# ----------------------------------------------------------------------- #


# *********************************************************************** #
#                        PREPARING READING OF FILES                       #
# *********************************************************************** #

# Parse command line
parser = argparse.ArgumentParser()

# Inputs
parser.add_argument("-d", "--directory", type=str, dest = "input_folder",\
    default = 'output', help="Name of the directory where the files to\
    analyze are located.")
parser.add_argument("-fn", "--file_name", type=str, dest = "file_name",\
    default = 'report_', help="Name of the report files.")
parser.add_argument("-T", "--temperature", type=float, dest = "temperature",\
    default = 298., help="Temperature of the experiment.")
parser.add_argument("-pS", "--pele_Steps", type=int, dest = "peleSteps",\
    default = 20, help="Number of Pele Steps in the simulation.")

args = parser.parse_args()

input_folder = args.input_folder
file_name = args.file_name
T = args.temperature
peleSteps = args.peleSteps

# Creating a variable with the path of the folder
path = str(pathlib.Path().absolute())
folderpath = path + '/' + input_folder

# Checking existance
if os.path.isdir(folderpath) == False:
    raise Exception('FolderPathError: There is no folder with this name. Please check the path and the folder name.')

# Listing directories inside the directory
files = os.listdir(folderpath)

# *********************************************************************** #
#                              READING FILES                              #
# *********************************************************************** #

be = []
step = []

for document in files:

    # Checking if the folder exist
    new_directory = os.path.join(folderpath,document)

    if os.path.isdir(new_directory) and document.isnumeric():

        # Listing files inside
        files = os.listdir(new_directory)

        if file_name in files == False:
            raise Exception('FilePathError: There is no file containing ' + file_name + ' in it. \
            Please check the path to the files and the files name.')
        
        # Iterating over files
        for file in files:
                
            cont = 0

            file_path = os.path.join(new_directory,file)

            # Checking whether path exists and the file has the file_name selected
            if os.path.isfile(file_path) and file_name in file:

                with open(file_path, 'r') as filein:

                    for line in filein:

                        if cont != 0:

                            line = line.split('   ')
                            be.append(float(line[4]))
                            step.append(int(line[1]))

                        cont += 1


minimum_energy = min(be)
be = np.array(be)
be4 = be/4

# *********************************************************************** #
#                 CALCULATING BOLTZAMNN WEIGHTED ENERGY                   #
# *********************************************************************** #

exp_bz = np.exp(-be/(R*T)) 
nominator = be.dot(exp_bz)
denominator = np.sum(exp_bz)
ene_bz = nominator/denominator

# *********************************************************************** #
#           CALCULATING BOLTZAMNN WEIGHTED CORRECTED ENERGY               #
# *********************************************************************** #

exp_bz4 = np.exp(-be4/(R*T)) 
nominator4 = be4.dot(exp_bz4)
denominator4 = np.sum(exp_bz4)
ene_bz4 = nominator4/denominator4

# *********************************************************************** #
#                          CALCULATING STEP WEIGHT                        #
# *********************************************************************** #

num_steps = []

for i in range(len(step) - 1):
 
    if step[i] == 0 and step[i+1] != 0:

        num_steps.append(step[i+1] - step[i])

    elif step[i] == 0 and step[i+1] == 0:

        num_steps.append(peleSteps)
        
    elif step[i] != 0 and step[i+1] == 0:

        num_steps.append(peleSteps-step[i])

    else:

        num_steps.append(step[i+1] - step[i])

num_steps.append(peleSteps - step[len(step)-1])

num_steps = np.array(num_steps)   
numerator = be.dot(num_steps)
denominator = np.sum(num_steps)
ene_step = numerator/denominator


# *********************************************************************** #
#                 CALCULATING STEP-BOLTZAMNN WEIGHTED ENERGY              #
# *********************************************************************** #

step_bz = num_steps.dot(exp_bz)
nominator = np.sum(np.multiply.reduce((be, num_steps, exp_bz)))
denominator = num_steps.dot(exp_bz)
ene_step_bz = nominator/denominator

# *********************************************************************** #
#                           PRINTING RESULTS                              #
# *********************************************************************** #

print(' ')
print('**************************************************************')
print('*                        ENERGIES                            *')
print('**************************************************************')
print(' ')
print('--------------------------------------------------------------')
print('- Minimum Binding Energy:', minimum_energy)
print('- Average Binding Energy:', np.average(be))
print('- Boltzmann weighted Energy:', ene_bz)
print('- Step weighted Energy:', ene_step)
print('- Step-Boltzmann weighted Energy:', ene_step_bz)
print('- Boltzmann weighted corrected Energy:', ene_bz4)
print('-------------------------------------------------------------- ')
print(' ')

with open('energy.csv', 'w') as fileout:
    fileout.writelines(
    'Minimum,Average,Boltzmann weighted,Step weighted,Step-Boltzmann weighted,Boltzmann weighted corrected\n' 
    '' + str(minimum_energy) + ',' + str(np.average(be)) + ',' + str(ene_bz) + ',' + str(ene_step) + ',' + str(ene_step_bz) + ',' + str(ene_bz4) + '\n'
    )