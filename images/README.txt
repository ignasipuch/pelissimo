
********************************************************************************************
Simulations performed with two different variables in three different targets and 8 systems
********************************************************************************************

	- Variables:
		1. Forcefield: OPLS2005 (OPLS) or Open Forcefield 1.2.0	(OFF)
		2. Protocol: 
			2.a) normal: 
				2.a.i) peleSteps: 25
				2.a.ii) iterations: 20
			2.b) fast:
				2.a.i) peleSteps: 10
				2.a.ii) iterations: 10

	- Targets:
		1. WT-HIV
		2. Galectin-3C
		3. MTAP

	- Systems PDB:
 		WT-HIV:
			1. 1HPV
			2. 1HSG
			3. 1MSM
			4. 1T3R
		Galectin-3C:
			5. 6QGE
			6. 6QGF
		MTAP:
			7. 1CB0
			8. 1K27

********************************************************************************************
Six scoring functions
********************************************************************************************

	1. Minimum Binding Energy -> min
	2. Average Binding Energy -> av
	3. Boltzmann weighted average -> bz
	4. Step weighted average -> st
	5. Step and Boltzamnn weighted average -> stbz
	6. Boltzmann weighted average on Binding energy corrected by a 1/4 factor -> bz_corr

********************************************************************************************
Two experimental values to compare to:
********************************************************************************************

	1. DG -> Experimental Binding free energy (Gibbs free energy)
	2. DH -> Experimental Binding enthalpy



				