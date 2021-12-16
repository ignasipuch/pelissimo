# PhD
## Development of a energy correctors to PELE predictions.

---
### alchem: 
#### Directory with all the code related to alchemPELE. Contains:

 1. alchem_a.py: Script to retrieve the reports related to AdaptivePELE's trajectories and perform an integration
of all of them with the Trapezoidal rule.

### binding_energy: 
#### Directory with the code related to calculating the binding energy of a simulation. Contains: 

 1. statistics.py: Script to read all the energetic information of a simulation and calculating different scoring functions.

### strain: 
#### Directory with the code related to the calculation of the LInEN (Ligand Internal Energy) with different approaches (ConACry). Contains:

 1. analysis

   a. linen_a.py: Script to calculate LInEN with the platform's analysis of the PELE simulation.

   b. LIG_Pele: Directory to test the script.

   c. LIG_linen_a: Directory to test the script.

 2. conformer:
 
   a. linen_con.py: Script to calculate LInEN with conformation ensamble generation of the ligand.
   
 3. crystal:
 
   a. linen_cry.py: Script to calculate LInEN with a PELE simulation of the ligand. 



