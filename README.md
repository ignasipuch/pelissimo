# Pelissimo

This repository contains certain python scripts to automatize specific [PELE](https://pele.bsc.es/pele.wt) simulations and some others to analyze PELE simulations.

## Contents


 * <b> 1. binding_energy:</b> Directory with the code related to calculating the binding energy of a simulation. Contains: 
  * **1.1. disc.py**: Script to: (1) Read all the information about a metric in a simulation and calculating different scoring functions.
 (2) Make plots of the evolution of the Boltzmann weighted average or minimum of a certain metric throughout a simulation.

  * **1.2. bootstrapping.py**: Script to generate bootstrap datasets and have statistics.

  * **1.3. pele_fetcher.py**: Script that fetches specific snapshot of a simulation.

 * <b> 2. corrector:</b> Directory containing the code involved in the correction of the energies of a simulation.

  * **2.1. correction.py**: Script that retrieves entropy and strain corrections and implements them in the reports of the > induced fit simulation.

 * <b> 3. entropy:</b> Directory with the code related to the Ligand Conformational Entropy (LiCE) calculation

  * **3.1 lice.py**: Script that calculates the entropy change upon binding. Requires pele platform analysis
 of an induced fit simulation and an isolated ligand simulation.
 
  * **3.2. dihedral_clustering.py**: Script that calculates dihedral angles for all the rotatable bonds in all the conformations reached 
 by the ligand in a simulation. This information will be used to cluster conformations. 

  * **3.3. entropy.py**: Script that ensembles the other two to make a complete entropy analysis of a simulation.

 * <b> 4. initial_files:</b> Directory with the files necessary to perform a full simulation of the protocol.

 * <b> 5. strain:</b> Directory with the code related to the calculation of the Ligand Internal Energy LInEn. Contains:
 
  * **5.1. linen.py**: Script to calculate LInEN with a PELE simulation of the ligand. 

 * <b> 6. structural:</b> Directory with the code related to the rmsd calculation of the protein-ligand induced fit simulation. Contains:

  * **6.1. rmsd.py**: Script that prepares all to run an rmsd calculation. 
