from schrodinger import structure
from schrodinger.application.jaguar.input import JaguarInput
from schrodinger.structure import StructureReader
from schrodinger.job.jobcontrol import launch_job

import argparse
import sys
import os
import shutil

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

    parser.add_argument("-f", "--file", type=str, dest="input_file",
                        default=None, required=True, help="Name of the file corresponding to the isolated ligand with connectivity.")
    
    parsed_args = parser.parse_args(args)

    return parsed_args

def multiplicity(input_file):
    """
    Function
    ----------
    Calculate or retrieve the number of electrons and multiplicity.

    Parameters
    ----------
    - input_file : str
        Ligand file in pdb format.

    Returns
    ----------
    - electrons : int
        Number of electrons calculated with rdkit (without double bonds)  
    - parity_bool : bool
        Bool that indicates if the number of electrons is odd or even.
    """

    outfile = input_file.split('.')[0] + '_qm_charges.out'
    electrons = 0

    if os.path.isfile(outfile):

        print('\n -     Reading .out file.')

        error_bool = False

        with open(outfile) as filein:
            for line in filein:
                if 'ERROR' in line:
                    error_bool = True

        if error_bool:

            with open(outfile) as filein:
                for line in filein:
                    if 'Total' in line and 'number' in line and 'electrons' in line:
                        line = line.split()
                        electrons = int(line[4])

            if electrons%2 == 0:
                parity_bool = True
            else: parity_bool = False

    else:

        from rdkit import Chem

        print("\n -   This code assumes there are no double bonds.\n If that is not the case, please re-run this code.")

        m = Chem.MolFromPDBFile(input_file)
        m1 = Chem.AddHs(m)

        for atom in m1.GetAtoms():
            electrons += atom.GetAtomicNum()

        if electrons%2 == 0:
            parity_bool = True
        else: parity_bool = False

    return electrons, parity_bool

def jaguar_input(input_file,
                 electrons,
                 parity_bool):
    """
    Function
    ----------
    Prepare input for the schrödinger job.

    Parameters
    ----------
    - input_file : str
        Ligand file in pdb format.
    - electrons : int
        Number of electrons calculated with rdkit (without double bonds)  
    - parity_bool : bool
        Bool that indicates if the number of electrons is odd or even.

    Returns
    ----------
    - jaguar_input_file : str
        Name of the output file the geometry optimization will yield.
    """

    print(" -     Writing .in file.")

    if parity_bool:

        print("     -     Spin multiplicity set to 1 due to even number of electrons (%a)." %electrons)

        options = {
               "igeopt": "1",        
               "isolv" : "2",
               "maxitg" : "10",
               "basis" : "6-31g**",
               "igeopt" : "1",
               "dftname" : "B3LYP-D3",
               "icfit" : "1",
               "nogas" : "1"   
              }

    else: 

        print("     -     Spin multiplicity set to 2 due to odd number of electrons (%a)." %electrons)

        options = {
               "igeopt": "1",        
               "isolv" : "2",
               "maxitg" : "10",
               "basis" : "6-31g**",
               "igeopt" : "1",
               "dftname" : "B3LYP-D3",
               "icfit" : "1",
               "nogas" : "1",
               "multip" : "2"
              }

    jaguar_input_file = str(input_file.split('.')[0] + "_qm_charges.in")

    st = next(StructureReader(input_file))
    ji = JaguarInput(structure=st, genkeys=options)
    ji.saveAs(jaguar_input_file)

    return jaguar_input_file

def jaguar_job(jaguar_input):
    """
    Function
    ----------
    Launches the schrödinger job.

    Parameters
    ----------
    - jaguar_input : str
        Name of the job's input.
    """

    run_command = ["jaguar", "run", jaguar_input]
    print(" -     Running Jaguar optimization under jobcontrol...")
    job = launch_job(run_command)
    job.wait()

def jaguar_charges(optimization_file,
                   charges_file):
    """
    Function
    ----------
    Prepare input for charge calculation and launch schrödinger job.

    Parameters
    ----------
    - optimization_file : str
        Name of the out file of the schrödinger optimization.
    - charges_file : str
        Name of the input file used in the charge calculation.  
    """

    with open(optimization_file, 'r') as filein:
        with open(charges_file, 'w') as fileout:
            for line in filein:
                if 'igeopt=1' in line:
                    fileout.write('igeopt=0\n')
                else:
                    fileout.write(line)

    run_command = ["jaguar", "run", charges_file]
    print("\n -     Running Jaguar charges under jobcontrol...")
    job = launch_job(run_command)
    job.wait()

def jaguar_output(output_file):
    """
    Function
    ----------
    Prepare directory to parametrize the ligand with peleffy.

    Parameters
    ----------
    - output_file : str
        Name of the output file of the schrödinger optimization.
    """

    def file_generation(system,
                        pdb_file,
                        output_file):
        """
        Function
        ----------
        Write the run file for MN4.

        Parameters
        ----------
        - system : str
            Name of the PDB the ligand is extracted from.

        - pdb_file : str
            Name of the PDB we are going to use to parametrize the ligand.

        - output_file : str
            Name of the outputfile we are going to be using for the charges.

        """

        with open('run', 'w') as fileout:

            fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH -J ' + system + '_jag\n'
                '#SBATCH --output=jag.out\n'
                '#SBATCH --error=jag.err\n'
                '#SBATCH --ntasks=48\n'
                '#SBATCH --qos=debug\n'
                '#SBATCH --time=00-01:00:00\n'
                '\n'
                'module purge\n'
                'module load ANACONDA/2019.10\n'
                '\n'
                'eval "$(conda shell.bash hook)"\n'
                '\n'
                'conda activate /gpfs/projects/bsc72/conda_envs/platform/1.6.3\n'
                '\n'
                'python -m peleffy.main ' + pdb_file + ' -f OPLS2005 --charges_from_file ' + output_file + '\n'
            )

    system = output_file.split('_qm_charges_POP.out')[0]
    pdb_file = output_file.split('.out')[0] + '.01.pdb'
    mae_file = output_file.split('.out')[0] + '.01.mae'
    path_results = system + '_jag'

    print(' ')
    os.system('$SCHRODINGER/utilities/structconvert -imae ' + mae_file + ' -opdb ' + pdb_file)

    file_generation(system,pdb_file,output_file)

    if os.path.exists(path_results) == False:
        os.mkdir(path_results)

    print('\n -     Job prepared in ' + str(path_results) + '.\n')

    shutil.copyfile(output_file, os.path.join(path_results,output_file))    
    shutil.copyfile('run', os.path.join(path_results,'run'))    
    shutil.copyfile(pdb_file, os.path.join(path_results,pdb_file))

def assemble(input_file):
    """
    Function
    ----------
    Function that assembles all the other functions to easy-to-tead format.

    Parameters
    ----------
    - input_file : str
        Ligand file in pdb format.
    """

    optimization_file = input_file.split('.')[0] + '_qm_charges.01.in'
    charges_file = input_file.split('.')[0] + '_qm_charges_POP.in'
    output_file = input_file.split('.')[0] + '_qm_charges_POP.out'

    if not os.path.isfile(optimization_file):

        electrons, parity_bool = multiplicity(input_file)
        jaguar_input_file = jaguar_input(input_file = input_file,
                                     electrons = electrons,
                                     parity_bool = parity_bool)
        jaguar_job(jaguar_input_file)

    if not os.path.isfile(output_file):

        jaguar_charges(optimization_file,
                       charges_file)

    jaguar_output(output_file)

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

    assemble(args.input_file)

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)


