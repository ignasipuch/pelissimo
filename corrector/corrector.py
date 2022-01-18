# -*- coding: utf-8 -*-
"""
This module is designed to correct energies of a simulation.
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import os
import pathlib
import argparse
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

    parser.add_argument("-d", "--directory", type=str, dest="input_folder",
                        default='LIG_Pele', help="Name of the directory where the simulation\
        is located.")
    parser.add_argument("-r", "--residue_name", type=str, dest="residue_name",
                        default='LIG', help="Ligand's residue name.")
    parser.add_argument("-rn", "--report_name", type=str, dest="report_name",
                        default='report', help="Name of the report files used for the simulation.")
    parser.add_argument("-cbe", "--column_binding_energy", type=int, dest="column_binding_energy",
                        default=None, help="Column of the report where the binding energy\
                             metric is located.")
    parser.add_argument("-cie", "--column_internal_energy", type=int, dest="column_internal_energy",
                        default=None, help="Column of the report where the internal energy of the ligand\
                             metric is located.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def corrector(input_folder,
              residue_name,
              report_name,
              column_binding_energy,
              column_internal_energy):
    """
    Function
    ----------
    With all the corrections calculated, they are introduced in the reports.

    Parameters
    ----------
    - input_folder : str
        The path to the directory created by the induced fit simulation.
    - residue_name : str
        Residue name of the ligand in the pdb of each cluster
    -report_name : str
        Name of the report files we want to correct.
    - column_binding_energy : int
        Column of the report where the binding energy metric is located.
    - column_internal_energy : int
        Column of the report where the internal energy metric is located.
    """

    def path_definer(input_folder,
                     residue_name):
        """
        Function
        ----------
        Defines all the paths that are going to be used

        Parameters
        ----------
        - input_folder : str
            The path to the directory created by the induced fit simulation.
        - residue_name : str
            Residue name of the ligand in the pdb of each cluster
        - clusters_folder : str
            Name of the directory where the directory clusters is located (results/analysis).

        Returns
        ----------
        - path_pl_simulation: str
            The path to the protein-ligand simulation.
        - path_pl_output : str
            The path to the protein-ligand simulation output.
        - path_l_simulation : str
            The path to the ligand simulation.
        """

        path = str(pathlib.Path().absolute())

        path_pl_simulation = os.path.join(path, input_folder)
        path_pl_output = os.path.join(path_pl_simulation, 'output')

        path_l_simulation = os.path.join(path, residue_name + '_linen_cry')

        if os.path.isdir(path_pl_simulation) == False:
            raise Exception('PathError: There is no folder with this path: ' +
                            path_pl_simulation + '. Please check the path and the folder name.')

        return path_pl_simulation, path_pl_output, path_l_simulation

    def corrections_detector(path_pl_simulation,
                             path_l_simulation):
        """
        Function
        ----------
        Checks whether there is a correction to be implemented to the reports of the 
        protein ligand simulation.

        Parameters
        ----------
        - path_pl_simulation : str
            The path to the protein-ligand simulation.
        - path_l_simulation : str
            The path to the ligand simulation.

        Returns
        ----------
        - correction_number : int
            Number that indicates univocally what corrections are to be implemented.

        """

        correction_number = 0

        if os.path.isfile(os.path.join(path_l_simulation, 'energy.csv')) == True:

            print(' -   Strain correction found.')
            correction_number += 1

        else:

            print(' -   Strain correction couldn\'t be found.'
                  '     File containing energy information (energy.csv)\n'
                  '     has not been found in\n\n'
                  '    ' + path_l_simulation + '.\n')

        if os.path.isfile(os.path.join(path_pl_simulation, 'entropy.csv')) == True:

            print(' -   Enropy correction found.')
            correction_number += 2

        else:

            print(' -   Entropy correction couldn\'t be found.\n'
                  '     File containing entropy information (entropy.csv)\n'
                  '     has not been found in\n\n' + path_pl_simulation + '.\n')

        if correction_number == 0:

            raise Exception('NoCorrectionToImplement: Either there is no correction to apply to the data' +
                            ' or they have not been found in the paths they are supposed to. (1) ' + path_pl_simulation + ' (2) ' + path_l_simulation + '.')

        return correction_number

    def ligand_min(path_l_simulation):
        """
        Function
        ----------
        Retrieves the energetic scoring of the ligand (we have chosen).

        Parameters
        ----------
        - path_l_simulation : str
            The path to the ligand simulation.

        Returns
        ----------
        - ligand_min_energy : float
            Ligand's minimum energy of the scoring function chosen.
        """

        cont = 0

        with open(os.path.join(path_l_simulation, 'energy.csv')) as filein:

            for line in filein:

                if cont != 0:

                    # This could change depending on the user.
                    line = line.split(',')
                    ligand_min_energy = float(line[2])

                cont += 1

        return ligand_min_energy

    def entropy_correction(path_pl_simulation):
        """
        Function
        ----------
        Retrieves the entropic change of the ligand.

        Parameters
        ----------
        - path_pl_simulation : str
            The path to the protein-ligand simulation.

        Returns
        ----------
        - entropy_change : float
            Entropy change calculated from the protein-ligand and the ligand simulations.
        """

        cont = 0

        with open(os.path.join(path_pl_simulation, 'entropy.csv')) as filein:

            for line in filein:

                if cont != 0:

                    line = line.split(',')
                    entropy_change = float(line[3])

                cont += 1

        return entropy_change
       
    def column_retriever(file,
                         column_binding_energy,
                         column_internal_energy):
        """
        Function
        ----------
        Retrieves the position of the binding energy.

        Parameters
        ----------
        - file : list
            The path to the one report.
        - column : int 
            Value give to the --column flag. 

        Returns
        ----------
        - column : int 
            Column where the interesting data is located. 
        """

        cont = 0

        with open(file, 'r') as filein:

            for line in filein:

                if cont == 0:

                    line = line.split('    ')
                    column_binding_energy = line.index('BindingEnergy') + 1
                    column_internal_energy = line.index('InternalEnergy') + 1

                    #
                    print(
                        ' -   Column for the binding energy:', column_binding_energy - 1)
                    print(
                        ' -   Column for the internal energy:', column_internal_energy - 1)
                    #

                cont += 1

        return column_binding_energy, column_internal_energy

    def correction_implementer(column_binding_energy,
                               column_internal_energy,
                               path_pl_output,
                               report_name,
                               ligand_min_energy,
                               entropy_change):
        """
        Function
        ----------
        Correct the reports with the calculated corrections at disposal.

        Parameters
        ----------
        - column_binding_energy : int
            Column of the report where the binding energy metric is located.
        - column_internal_energy : int
            Column of the report where the internal energy metric is located.
        - path_pl_output : str
            The path to the protein-ligand simulation output.
        - report_name : str
            Name of the report files we want to correct.
        - ligand_min_energy : float
            Ligand's minimum energy of the scoring function chosen.
        - entropy_change : float
            Entropy change calculated from the protein-ligand and the ligand simulations.
        """

        def copier(path_pl_output,
                   report_name):
            """
            Function
            ----------
            Copies al the reports in place with the prefix mod_ added. It also 
            stores information for the correction.

            Parameters
            ----------
            - path_pl_output : str
                The path to the protein-ligand simulation output.
            - report_name : str
                Name of the report files we want to correct.

            Returns
            ----------
            - report_paths : list
                List of all the reports' paths.
            """

            cont_reports = 0
            report_paths = []

            # Copying reports
            if os.path.isdir(path_pl_output):

                files = os.listdir(path_pl_output)

                for folder in files:

                    if folder.isnumeric():

                        full_path = os.path.join(path_pl_output, folder)

                    files_subdir = os.listdir(full_path)

                    for report in files_subdir:

                        if report.startswith(report_name) and report.split(report_name + '_')[1].isnumeric():

                            shutil.copy(os.path.join(full_path, report),
                                        os.path.join(full_path, 'mod_' + report_name + '_' + report.split(report_name + '_')[1]))

                            report_paths.append(
                                os.path.join(full_path, report))
                            cont_reports += 1

                        else:
                            continue

                    if cont_reports == 0:

                        raise Exception('ReportNameError: No reports beginning with \"' + report_name + '\" were found in '
                                        + full_path + '.')

            return report_paths

        def corrector(column_binding_energy,
                      column_internal_energy,
                      report_paths,
                      report_name,
                      ligand_min_energy,
                      entropy_change):
            """
            Function
            ----------
            Correct the new reports with the calculated corrections.

            Parameters
            ----------
            - column_binding_energy : int
                Column of the report where the binding energy metric is located.
            - column_internal_energy : int
                Column of the report where the internal energy metric is located.
            - report_paths : list
                List of all the reports' paths.
            - report_name : str
                Name of the report files we want to correct.
            - ligand_min_energy : float
                Ligand's minimum energy of the scoring function chosen.
            - entropy_change : float
                Entropy change calculated from the protein-ligand and the ligand simulations.
            """

            for report_path in report_paths:

                path_modified_report = report_path.replace(
                    report_name, 'mod_' + report_name)

                with open(path_modified_report, 'w') as fileout:

                    cont = 0

                    with open(report_path) as filein:

                        for line in filein:

                            if cont == 0:

                                fileout.write(line)

                            else:

                                line = line.split()
                                strain_energy = float(
                                    line[column_internal_energy-1]) - ligand_min_energy
                                line[column_binding_energy-1] = \
                                    str(float(
                                        line[column_binding_energy-1]) + strain_energy + entropy_change)
                                fileout.write("     ".join(line) + '\n')

                            cont += 1

        # Copying  reports
        report_paths = copier(path_pl_output, report_name)

        # Correcting copied reports
        corrector(column_binding_energy,
                  column_internal_energy,
                  report_paths,
                  report_name,
                  ligand_min_energy,
                  entropy_change)

    def analysis_files_writer(column_binding_energy,
                              path_pl_simulation):
        """
        Function
        ----------
        Write files for further analysis with pele platform.

        Parameters
        ----------
        - column_binding_energy : int
            Column of the report where the binding energy metric is located.
        - path_pl_simulation: str
            The path to the protein-ligand simulation.
        """

        with open(os.path.join(path_pl_simulation, 'run_analysis'), 'w') as fileout:

            fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH --job-name=analysis\n'
                '#SBATCH --output=analysis.out\n'
                '#SBATCH --error=analysis.err\n'
                '#SBATCH --ntasks=48\n'
                '#SBATCH --qos=debug\n'
                '#SBATCH --time=00-00:30:00\n'
                '\n'
                'module load ANACONDA/2019.10\n'
                'module load intel mkl impi gcc # 2> /dev/null\n'
                'module load impi\n'
                'module load boost/1.64.0\n'
                '\n'
                'eval "$(conda shell.bash hook)"\n'
                'conda activate /gpfs/projects/bsc72/conda_envs/platform/1.6.2\n'
                '\n'
                'python script.py\n'
            )

        with open(os.path.join(path_pl_simulation, 'script.py'), 'w') as fileout:

            fileout.writelines(
                'from pele_platform.analysis import Analysis\n'
                '\n'
                'analysis = Analysis(resname="' + residue_name +
                '", chain="L", simulation_output="output", be_column = ' + str(column_binding_energy) + ', report="' +
                'mod_' + report_name + '", cpus=48)\n'
                'analysis.generate(path="analysis", clustering_type="meanshift")\n'
            )

    #
    print(' ')
    print('*******************************************************************')
    print('*                        peleCorrector                            *')
    print('* --------------------------------------------------------------- *')
    print('*                   Corrector of induced fit results              *')
    print('*******************************************************************')
    print(' ')
    #

    path_pl_simulation, path_pl_output, path_l_simulation = path_definer(
        input_folder, residue_name)

    correction_number = corrections_detector(
        path_pl_simulation, path_l_simulation)

    # Corrections values ---
    if correction_number == 1:

        ligand_min_energy = ligand_min(path_l_simulation)
        entropy_change = 0.

    elif correction_number == 2:

        ligand_min_energy = 0.
        entropy_change = entropy_correction(path_pl_simulation)

    elif correction_number == 3:

        ligand_min_energy = ligand_min(path_l_simulation)
        entropy_change = entropy_correction(path_pl_simulation)

    # Corrections column location ---
    file = os.path.join(path_pl_simulation,'output','0',report_name + '_1')

    if column_internal_energy == None and column_binding_energy == None:

        column_binding_energy, column_internal_energy = \
        column_retriever(file,
                         column_binding_energy,
                         column_internal_energy)

    elif column_binding_energy == None:
        
        column_binding_energy, _ = \
        column_retriever(file,
                         column_binding_energy,
                         column_internal_energy)

    elif column_internal_energy == None:

        _, column_internal_energy = \
        column_retriever(file,
                         column_binding_energy,
                         column_internal_energy)

    #
    print(' ')
    print(' -   Implementing corrections...')
    #

    correction_implementer(column_binding_energy,
                           column_internal_energy,
                           path_pl_output,
                           report_name,
                           ligand_min_energy,
                           entropy_change)

    #
    print(' -   Job finished succesfully. Energies corrected will be found in \n'
          '     mod_' + report_name + ' files.')
    #

    analysis_files_writer(column_binding_energy,
                          path_pl_simulation)

    #
    print(' -   run_analysis and script.py files have been generated.')
    print(' ')
    print('------------------------------ INFO -------------------------------')
    print(' ')
    print(' -   To perform a new analysis:')
    print('     :> cd ' + input_folder)
    print('     :> sbatch run_analysis')
    print(' -   The new plots and files will be in /analysis.')
    print(' ')
    print('-------------------------------------------------------------------')
    print(' ')


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

    corrector(input_folder=args.input_folder,
              residue_name=args.residue_name,
              report_name=args.report_name,
              column_binding_energy=args.column_binding_energy,
              column_internal_energy=args.column_internal_energy)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
