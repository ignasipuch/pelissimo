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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import OrderedDict


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
    parser.add_argument("-q", "--quantile", type=float, dest="quantile",
                    default=1, help="Percentage of data with lowest total energy snapshots you want\
        to keep to assess the strain energy")

    parsed_args = parser.parse_args(args)

    return parsed_args


def corrector(input_folder,
              residue_name,
              report_name,
              quantile):
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
        - path_pl_results : str
            The path to the folder where results will be stored.
        """

        path = str(pathlib.Path().absolute())

        path_pl_simulation = os.path.join(path, input_folder)
        path_pl_output = os.path.join(path_pl_simulation, 'output')
        path_pl_results = os.path.join(path_pl_simulation, 'strain')

        path_l_simulation = os.path.join(path, residue_name + '_linen_cry')

        if os.path.isdir(path_pl_simulation) == False:
            raise Exception('PathError: There is no folder with this path: ' +
                            path_pl_simulation + '. Please check the path and the folder name.')

        if os.path.isdir(path_pl_results) == False:
            os.mkdir(path_pl_results)

        return path_pl_simulation, path_pl_output, path_l_simulation, path_pl_results

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

            print(' -   Entropy correction found.')
            correction_number += 2

        else:

            print(' -   Entropy correction couldn\'t be found.\n'
                  '     File containing entropy information (entropy.csv)\n'
                  '     has not been found in\n')
            print('     ' + path_pl_simulation)
            print(' ')

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
                    ligand_min_energy = float(line[0])

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

    def column_retriever(file):
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
                    column_current_energy = line.index('currentEnergy') + 1
                    column_binding_energy = line.index('BindingEnergy') + 1
                    column_internal_energy = line.index('InternalEnergy') + 1

                    #
                    print(
                        ' -   Column for the binding energy:', column_binding_energy - 1)
                    print(
                        ' -   Column for the internal energy:', column_internal_energy - 1)
                    #

                cont += 1

        return column_current_energy, column_binding_energy, column_internal_energy

    def correction_implementer(column_current_energy,
                               column_binding_energy,
                               column_internal_energy,
                               path_pl_output,
                               report_name,
                               ligand_min_energy,
                               entropy_change):
        """
        Function
        ----------
        Correct the reports with the calculated corrections at disposal and 
        store simulation data.

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

        Returns
        ----------
        - strain_energy_list : list
            List of all the strain energies calculated.
        - simulation_df : pd.DataFrame
            Data frame with currentEnergy and strainEnergy values of the simulation.
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

            Returns
            ----------
            - strain_energy_list : list
                List with all the strain energies calculated in a simulation.
            """

            epochs = []
            trajectories = []
            models = []
            currentEnergies = []
            strainEnergies = []

            for report_path in report_paths:

                epoch = report_path.split('/')[-2]
                trajectory = report_path.split('/')[-1].split('_')[-1]

                path_modified_report = report_path.replace(
                    report_name, 'mod_' + report_name)

                with open(path_modified_report, 'w') as fileout:

                    cont = 0

                    with open(report_path) as filein:

                        for line in filein:

                            if cont == 0:

                                fileout.write(line)

                            else:

                                if '/0/' + report_name in report_path:

                                    if cont == 1:

                                        line = line.split()
                                        strain_energy = float(
                                            line[column_internal_energy-1]) - ligand_min_energy
                                        line[column_binding_energy-1] = \
                                            str(float(
                                                line[column_binding_energy-1]) + strain_energy + entropy_change)
                                        fileout.write(
                                            "     ".join(line) + '\n')

                                    else:

                                        line = line.split()

                                        currentEnergy = line[column_current_energy-1]
                                        model = line[2]
                                        strain_energy = float(
                                            line[column_internal_energy-1]) - ligand_min_energy

                                        line[column_binding_energy-1] = \
                                            str(float(
                                                line[column_binding_energy-1]) + strain_energy + entropy_change)
                                        fileout.write(
                                            "     ".join(line) + '\n')

                                        epochs.append(int(epoch))
                                        trajectories.append(int(trajectory))
                                        models.append(int(model))
                                        currentEnergies.append(
                                            float(currentEnergy))
                                        strainEnergies.append(
                                            float(strain_energy))

                                else:

                                    line = line.split()

                                    currentEnergy = line[column_current_energy-1]
                                    model = line[2]
                                    strain_energy = float(
                                        line[column_internal_energy-1]) - ligand_min_energy

                                    line[column_binding_energy-1] = \
                                        str(float(
                                            line[column_binding_energy-1]) + strain_energy + entropy_change)
                                    fileout.write("     ".join(line) + '\n')

                                    epochs.append(int(epoch))
                                    trajectories.append(int(trajectory))
                                    models.append(int(model))
                                    currentEnergies.append(
                                        float(currentEnergy))
                                    strainEnergies.append(float(strain_energy))

                            cont += 1

            simulation_df = pd.DataFrame(OrderedDict([('epoch', epochs),
                                                      ('trajectory', trajectories),
                                                      ('model', models),
                                                      ('currentEnergy',
                                                       currentEnergies),
                                                      ('strainEnergy', strainEnergies)]))

            simulation_df = simulation_df.sort_values(
                ['epoch', 'trajectory', 'model', 'currentEnergy', 'strainEnergy'], ascending=True).reset_index(drop=True)

            return strainEnergies, simulation_df

        # Copying  reports
        report_paths = copier(path_pl_output, report_name)

        # Correcting copied reports and storing information
        strain_energy_list, simulation_df = corrector(column_binding_energy,
                                                      column_internal_energy,
                                                      report_paths,
                                                      report_name,
                                                      ligand_min_energy,
                                                      entropy_change)

        return strain_energy_list, simulation_df

    def quantile_retriever(simulation_df,
                          quantile):
        """
        Function
        ----------
        Obtain a list of strain energies corresponding to those conformations
        with energies below a certain quantile.

        Parameters
        ----------
        - simulation_df : pd.DataFrame
            Data frame with current energies and strain information.
        - quantile : float
            Percentage of data with lowest total energy snapshots you want to keep 
        to assess the strain energy.

        Returns
        ----------
        - strain_energy_quantile : list
            List with the strain energies belonging to the quantile percentage of data 
        with lowest total energy.
        """

        quantile_value = simulation_df.currentEnergy.quantile(quantile)
        simulation_df_quantile = simulation_df[simulation_df.currentEnergy < quantile_value]
        strain_energy_quantile = list(np.array(simulation_df_quantile['strainEnergy']))
        
        return strain_energy_quantile

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

    def results_writer(strain_energy_list,
                       strain_energy_quantile,
                       path,
                       report_name,
                       quantile,
                       ligand_min_energy):
        """
        Function
        ----------
        Writes and plots result files of strain values.

        Parameters
        ----------
        - strain_energy_list : list
            List with all the strain values calculated from the entire simulation.
        - path : str
            The path to the protein-ligand strain results folder.
        """

        def histogram_function(strain_energy):
            """
            Function
            ----------
            Writes and plots result files of strain values.

            Parameters
            ----------
            - strain_energy_list : list
                List with all the strain values calculated from the entire simulation.
            """

            strain_energy = [
                item for item in strain_energy if item < 200]

            strain_energy_vector = np.array(strain_energy)

            bin_edges = np.histogram_bin_edges(strain_energy_vector, bins='fd')
            density, _ = np.histogram(strain_energy_vector, bins=bin_edges)

            hist_ene = 0.5*(bin_edges[np.argmax(density)] +
                    bin_edges[np.argmax(density) + 1])
            minimum_ene = min(strain_energy_vector)
            average_ene = np.average(strain_energy_vector)
            max_ene = max(strain_energy_vector)

            return strain_energy_vector, bin_edges, hist_ene, minimum_ene, average_ene, max_ene
        
        if quantile == 1.:
            
            strain_energy_vector, bin_edges, hist_ene, minimum_ene, average_ene, max_ene = histogram_function(strain_energy_list)

            if min(strain_energy_vector) < 0.:

                strain_energy_vector = strain_energy_vector - minimum_ene

                #
                print('\n'
                      '                              WARNING:                               \n'
                      '   Lower ligand energies were found in the induced fit simulations.  \n'
                      '   The results in mod_' + report_name + ' have been corrected with   \n'
                      '   Ligand minimum = ' + "{:3.3f}".format(ligand_min_energy) + ' > ' +
                      "{:3.3f}".format(minimum_ene + ligand_min_energy) +
                      ' = Induced fit minimum.\n'
                      '   Better sampling of the ligand is required.          \n'
                      )
                #


            with open(os.path.join(path, 'strain.csv'), 'w') as fileout:
                fileout.writelines(
                    'Minimum,Histogram max,Average,Maximum\n'
                    '' + str(minimum_ene) + ',' + str(hist_ene) + ',' +
                    str(average_ene) + ',' + str(max_ene) + '\n')

            # Plot
            plt.title('Strain distribution')
            plt.hist(strain_energy_vector, bins=bin_edges, density=True, color = "blue")
            plt.xlabel('Strain (kcal/mol)')
            plt.ylabel('Density')
            plt.savefig(os.path.join(path, 'density_strain.png'), format='png')
        
        else: 

            print(' ')
            print(' -   Quantile chosen for the results:', quantile)

            strain_energy_vector, \
            bin_edges, \
            hist_ene,\
            minimum_ene,\
            average_ene,\
            max_ene = histogram_function(strain_energy_list)

            strain_energy_vector_q,\
            bin_edges_q,\
            hist_ene_q,\
            minimum_ene_q,\
            average_ene_q,\
            max_ene_q = histogram_function(strain_energy_quantile)

            with open(os.path.join(path, 'strain.csv'), 'w') as fileout:
                fileout.writelines(
                    'Minimum,Histogram max,Average,Maximum\n'
                    '' + str(minimum_ene) + ',' + str(hist_ene) + ',' +
                    str(average_ene) + ',' + str(max_ene) + '\n'
                    '' + str(minimum_ene_q) + ',' + str(hist_ene_q) + ',' +
                    str(average_ene_q) + ',' + str(max_ene_q) + '\n'
                )

            # Plot
            plt.title('Strain distribution')
            plt.hist(strain_energy_vector, bins=bin_edges, density=True, color = "blue", label='original')
            plt.xlabel('Strain (kcal/mol)')
            plt.ylabel('Density')
            plt.savefig(os.path.join(path, 'density_strain.png'), format='png')

            plt.title('Strain distribution quantile ' + str(int(quantile*100)) + '%')
            plt.hist(strain_energy_vector_q, bins=bin_edges_q, density=True, color = "red", label='quantile')
            plt.xlabel('Strain (kcal/mol)')
            plt.ylabel('Density')
            plt.legend(loc='best')
            plt.savefig(os.path.join(path, 'density_strain_quantile.png'), format='png')

    #
    print(' ')
    print('*******************************************************************')
    print('*                        peleCorrector                            *')
    print('* --------------------------------------------------------------- *')
    print('*                   Corrector of induced fit results              *')
    print('*******************************************************************')
    print(' ')
    #

    path_pl_simulation,\
        path_pl_output,\
        path_l_simulation,\
        path_pl_results = path_definer(input_folder, residue_name)

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
    file = os.path.join(path_pl_simulation, 'output', '0', report_name + '_1')

    column_current_energy, column_binding_energy, column_internal_energy = \
        column_retriever(file)

    #
    print(' ')
    print(' -   Implementing corrections...')
    print(' ')
    #

    strain_energy_list, simulation_df = correction_implementer(column_current_energy,
                                                              column_binding_energy,
                                                              column_internal_energy,
                                                              path_pl_output,
                                                              report_name,
                                                              ligand_min_energy,
                                                              entropy_change)

    strain_energy_quantile = quantile_retriever(simulation_df,
                                                quantile)

    #
    print(' -   Job finished succesfully. Energies corrected will be found in \n'
          '     mod_' + report_name + ' files.')
    #

    results_writer(strain_energy_list,
                   strain_energy_quantile,
                   path_pl_results,
                   report_name,
                   quantile,
                   ligand_min_energy)

    #
    print(' -   Report about strain values can be found in ' +
          input_folder + '/strain')
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
              quantile=args.quantile)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
