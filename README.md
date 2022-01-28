# PhD
## Development of a energy correctors to PELE predictions.

---
> ### 1. alchem: 
> #### Directory with all the code related to alchemPELE. Contains:
>
> - **1.1. alchem_a.py**: Script to retrieve the reports related to AdaptivePELE's trajectories and perform an integration
of all of them with the Trapezoidal rule.

> ### 2. binding_energy: 
> #### Directory with the code related to calculating the binding energy of a simulation. Contains: 
>
> - **2.1. disc.py**: Script to: (1) Read all the information about a metric in a simulation and calculating different scoring functions.
> (2) Make plots of the evolution of the Boltzmann weighted average or minimum of a certain metric throughout a simulation.
>
> - **2.2. kTr.py**: Script to test whether PELE introduces bias to the energies calculated in the simulations. 

> ### 3. corrector:
> #### Directory containing the code involved in the correction of the energies of a simulation.
>
> - **3.1. correction.py**: Script that retrieves entropy and strain corrections and implements them in the reports of the 
> induced fit simulation.

> ### 4. entropy:
> #### Directory with the code related to the Ligand Conformational Entropy (LiCE) calculation
>
> - **4.1 lice.py**: Script that calculates the entropy change upon binding. Requires pele platform analysis
> of an induced fit simulation and an isolated ligand simulation.
>
> - **4.2. LIG_PELE**: Directory tree with needed files to test the script.
>
> - **4.3. LIG_linen_cry**: Directory tree with needed files to test the script.

> ### 5. initial_files:
> #### Directory with the files necessary to perform a full simulation of the protocol.

> ### 6. strain: 
> #### Directory with the code related to the calculation of the Ligand Internal Energy LInEn with different approaches (ConACry). Contains:
> 
> - **6.1. linen_cry.py**: Script to calculate LInEN with a PELE simulation of the ligand. 

> ### 7. structural:
> #### Directory with the code related to the rmsd calculation of the protein-ligand induced fit simulation. Contains:
>
> - **7.1. rmsd.py**: Script that prepares all to run an rmsd calculation. 
