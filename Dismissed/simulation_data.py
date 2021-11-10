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
#                 READING AND STORING INFORMATION FROM FILES              #
# ----------------------------------------------------------------------- #

# Parse command line
parser = argparse.ArgumentParser()

# Inputs
parser.add_argument("-d", "--directory", type=str, dest = "input_folder",\
    default = None, help="Name of the directory where the files to\
    analyze are located.")

args = parser.parse_args()

input_folder = args.input_folder

# Creating a variable with the path of the folder
path = str(pathlib.Path().absolute())
folderpath = path + '/' + input_folder

# Checking existance
if os.path.isdir(folderpath) == False:
    raise Exception('FolderPathError: There is no folder with this name. Please check the path and the folder name.')

# Listing files inside the directory
files = os.listdir(folderpath)

# Print ---

print(' ')
print(' ----------------------------------------------')
print('|                READING INPUT                 |')
print(' ----------------------------------------------')
print(' ')
print('- Folder to analyze:',folderpath)
print('- Number of files to analyze:', len(files))
print(' ')
print('---------------------------------------------------------')

cont = 0
cluster = []
population = []
mbe_min = []
mbe_5 = []
mbe_25 = []
mbe = []
mbe_75 = []
mbe_95 = []
mbe_max = []
mbe_sd = []
label = []

enthalpy = []
entropy = []
be = []

enthalpy_exp = [0.0,0.2,0.0,0.2,0.0,0.2,0.0,0.2]
entropy_exp = [-7.6,-12.5,-7.6,-12.5,-7.6,-12.5,-7.6,-12.5]
be_exp = [-7.6,-12.3,-7.6,-12.3,-7.6,-12.3,-7.6,-12.3]

for document in files:

    cont = 0

    # Checking if the file and folder exist
    if os.path.isfile(os.path.join(folderpath,document)):

        document_data = document.split('_')
                
        print('---------------------------------------------------------')
        print(' ')
        print('Forcefield:',document_data[1])
        print('Sampling:',document_data[2])
        print('System:',document_data[3])
        print(' ')
        print('---------------------------------------------------------')

        label.append(document_data[1] + ' ' + document_data[2] + ' ' + document_data[3])

        # Setting the document path for each document 
        documentpath = folderpath + '/' + document
        
        with open (documentpath, 'r') as filein:
        
            for line in filein:
            
                if cont != 0:
                
                    line = line.split(',')
                    cluster.append(int(line[0]))
                    population.append(float(line[1]))
                    mbe_min.append(float(line[11]))
                    mbe_5.append(float(line[12]))
                    mbe_25.append(float(line[13]))
                    mbe.append(float(line[14]))
                    mbe_75.append(float(line[15]))
                    mbe_95.append(float(line[16]))
                    mbe_max.append(float(line[17]))
                    mbe_sd.append(float(line[18]))

                cont += 1

        num_clusters = max(cluster)
        print(mbe)

        S = 0.

        for i in range(num_clusters):
            S += population[i]*np.log(population[i])

        S = -1.9858775e-3*S
        T = 298

        print('---------------------------------------------------------')
        print('Entropic term:',-T*S)
        print('Minimum Mean Binding Energy:',min(mbe))
        print('Sum of terms:',-T*S + min(mbe))
        print('---------------------------------------------------------')

        entropy.append(-T*S)
        enthalpy.append(min(mbe))
        be.append(-T*S + min(mbe))

print(enthalpy)
print(entropy)
print(be)


fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)

ax1.plot(label,abs(np.array(enthalpy)-np.array(enthalpy_exp)/np.array(enthalpy)))
ax1.set_title('Enthalpy')
ax1.set_xticklabels(label,rotation=90)

ax2.plot(label,abs(np.array(entropy)-np.array(entropy_exp)/np.array(entropy)))
ax2.set_title('Entropy')
ax2.set_xticklabels(label,rotation=90)

ax3.plot(label,abs(np.array(be)-np.array(be_exp)/np.array(be)))
ax3.set_title('Binding Energy')
ax3.set_xticklabels(label,rotation=90)

plt.tight_layout()
plt.savefig('energy_figure.pdf')


