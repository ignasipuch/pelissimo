import os
import pathlib 
import numpy as np

path = str(pathlib.Path().absolute())

lineclusters = []
linenergies = []

with open (os.path.join(path,'PELEne.out'),'r') as filein:

    for line in filein:

        if 'CLUSTER' in line:

            line = line.split('CLUSTER')
            line = line[1].strip()
            lineclusters.append(line[0])
        
        if 'ENERGY VACUUM + SGB + CONSTRAINTS + SELF + NON POLAR:' in line:

            line = line.split('ENERGY VACUUM + SGB + CONSTRAINTS + SELF + NON POLAR:')
            linenergies.append(float(line[1].strip()))

    sorted_clusters = [x for _, x in sorted(zip(linenergies, lineclusters))]
    sorted_energies = sorted(linenergies)

sorted_energies_corrected = np.array(sorted_energies) - min(linenergies)

with open('energy.csv', 'w') as fileout:
    fileout.writelines(
    'Cluster,Internal energy,Internal energy change\n' 
    )

    for i in range(len(lineclusters)):

        fileout.write(str(sorted_clusters[i]) + ',' + str(sorted_energies[i]) + ',' + str(sorted_energies_corrected[i]) + '\n')

