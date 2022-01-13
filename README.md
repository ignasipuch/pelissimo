# PhD
## Development of a energy correctors to PELE predictions.

---
### 1. alchem: 
#### Directory with all the code related to alchemPELE. Contains:

 **1.1. alchem_a.py**: Script to retrieve the reports related to AdaptivePELE's trajectories and perform an integration
of all of them with the Trapezoidal rule.

### 2. binding_energy: 
#### Directory with the code related to calculating the binding energy of a simulation. Contains: 

 **2.1. statistics.py**: Script to read all the energetic information of a simulation and calculating different scoring functions.

 **2.2. kTr.py**: Script to test whether PELE introduces bias to the energies calculated in the simulations. 

### 3. corrector:
#### Directory containing the code involved in the correction of the energies of a simulation.

 **3.1. correction.py**: Script that retrieves entropy and strain corrections and implements them in the reports of the 
induced fit simulation.


### 4. entropy:
#### Directory with the code related to the Ligand Conformational Entropy (LiCE) calculation

 **4.1 lice.py**: Script that calculates the entropy change upon binding. Requires pele platform analysis
of an induced fit simulation and an isolated ligand simulation.

 **4.2. LIG_PELE**: Directory tree with needed files to test the script.

 **4.3. LIG_linen_cry**: Directory tree with needed files to test the script.


### 5. strain: 
#### Directory with the code related to the calculation of the Ligand Internal Energy LInEn with different approaches (ConACry). Contains:

 5.1. analysis

   **5.1.a. linen_a.py**: Script to calculate LInEN with the platform's analysis of the PELE simulation.

   **5.1.b. LIG_Pele**: Directory to test the script.

   **5.1.c. LIG_linen_a**: Directory to test the script.

 5.2. conformer:
 
   **5.2.a. linen_con.py**: Script to calculate LInEN with conformation ensamble generation of the ligand.
   
 5.3. crystal:
 
   **5.3.a. linen_cry.py**: Script to calculate LInEN with a PELE simulation of the ligand. 

### 6. structural:
#### Directory with the code related to the rmsd calculation of the protein-ligand induced fit simulation.. Contains:

 **6.1. rmsd.py**: Script that prepares all to run an rmsd calculation. 
