# -*- coding: utf-8 -*-
"""
This module is designed to calculate different scoring methods
for the binding energy of a simulation. Different Scorings Calculator
(DiSC).
"""

# Imports
import sys
import os
import pathlib
import argparse
import numpy as np
import time
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------- #
# Constants:
T = 298.
R = 1.985e-3
# ----------------------------------------------------------------------- #


def parse_args(args):
    """
    Function
    ----------
    It parses the command-line arguments.

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
                        default='output', help="Name of the output directory where the files to\
        analyze are located.")
    parser.add_argument("-rn", "--report_name", type=str, dest="report_name",
                        default='report', help="Name of the report files.")
    parser.add_argument("-T", "--temperature", type=float, dest="temperature",
                        default=298., help="Temperature of the experiment.")
    parser.add_argument("-pS", "--pele_Steps", type=int, dest="pele_steps",
                        default=None, help="Number of Pele Steps in the simulation.")
    parser.add_argument("-c", "--column", type=int, dest="column",
                        default=None, help="Column of the report where the interesting\
                             metric is located.")
    parser.add_argument("-a", "--action", type=str, dest="action",
                        default='all', help="Function the user wants the script to do: all or evolution.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def statistics(input_folder,
               report_name,
               T,
               pele_steps,
               column,
               action):
    """
    Function
    ----------
    Reads the data from all the reports and calculates different scoring functions.
    The information is stored in energy.csv

    Parameters
    ----------
    - input_folder : str
        The path to the directory where the output of the PELE simulation is located.
    - report_name : str
        Name of the reports containing the energetic data of all the simulation.
    - T : float
        Temperature to perform the Boltzmann weights with.
    - pele_steps : int
        Number of pele_steps in the simulation.
    - column : int 
        Column where the interesting data is located.
    """

    def reader(files,
               folderpath,
               report_name,
               column):
        """
        Function
        ----------
        Reads the data from all the reports and stores it in lists.

        Parameters
        ----------
        - files : list
            The path to the directory where the output of the PELE simulation is located.
        - folderpath : str
            Path where the different epochs of the simulation are located.
        - report_name : str
            Name of the reports to obtain the data from
        - column : int 
            Column where the interesting data is located.


        Returns
        ----------
        - be : list
            Binding energies of all the simulation.
        - step : list
            Steps associated to poses for all the simulation.
        """

        def file_reader(files,
                        folderpath,
                        report_name,
                        column):
            """
            Function
            ----------
            Reads the data from all the reports in a list.

            Parameters
            ----------
            - files : list
                The path to the directory where the output of the PELE simulation is located.
            - folderpath : str
                Path where the different epochs of the simulation are located.
            - report_name : str
                Name of the reports to obtain the data from
            - column : int 
                Column where the interesting data is located. 

            Returns
            ----------
            - be : list
                Binding energies of all the simulation.
            - step : list
                Steps associated to poses for all the simulation.
            """

            new_directory = folderpath

            for file in files:

                cont = 0

                if file.startswith(report_name):

                    file_path = os.path.join(new_directory, file)

                    with open(file_path, 'r') as filein:

                        for line in filein:

                            line = line.split('    ')

                            if cont != 0:

                                be.append(float(line[column-1]))
                                step.append(int(line[1]))

                            cont += 1

            return be, step

        def column_retriever(file,
                             column):
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
                        column = line.index('BindingEnergy') + 1

                        #
                        print(
                            ' -   No column introduced, picking BindingEnergy column:', column - 1)
                        #

                    cont += 1

            return column

        be = []
        step = []
        numeric_files = [s for s in files if s.isnumeric()]

        if len(numeric_files) != 0:

            if column == None:

                column_file = os.path.join(
                    folderpath, numeric_files[0], report_name + '_1')
                column = column_retriever(column_file,
                                          column)

            for document in numeric_files:

                # Checking if the folder exist
                new_directory = os.path.join(folderpath, document)

                if os.path.isdir(new_directory) and document.isnumeric():

                    # Listing files inside
                    files = os.listdir(new_directory)

                    if report_name in files == False:
                        raise Exception('FilePathError: There is no file containing ' + report_name + ' in it. \
                        Please check the path to the files and the files name.')

                    be, step = file_reader(files,
                                           new_directory,
                                           report_name,
                                           column)

        else:

            if report_name in files == False:
                raise Exception('FilePathError: There is no file containing ' + report_name + ' in it. \
                Please check the path to the files and the files name.')

            if column == None:

                column_file = os.path.join(folderpath, report_name + '1')
                column = column_retriever(column_file,
                                          column)

            be, step = file_reader(files,
                                   folderpath,
                                   report_name,
                                   column)

        return be, step

    def boltzmann_weighted(be,
                           T,
                           steps=[]):
        """
        Function
        ----------
        Calculates boltzmann weighted energy.

        Parameters
        ----------
        - be : list
            Binding energies of all the simulation.
        - T : float
            Temperature to perform the Boltzmann weights with.
        - steps : list
            Steps associated to poses for all the simulation.

        Returns
        ----------
        - ene_bz : float
            Value of the boltzmann weighted energy.
        """
        exp_bz = np.exp(-be/(R*T))

        if not steps:

            nominator = be.dot(exp_bz)
            denominator = np.sum(exp_bz)
            ene_bz = nominator/denominator

        else:

            steps = np.array(steps)
            nominator = np.sum(np.multiply.reduce((be, steps, exp_bz)))
            denominator = steps.dot(exp_bz)
            ene_bz = nominator/denominator

        return ene_bz

    def pelesteps_retriever():
        """
        Function
        ----------
        Reads the pelesteps from adaptive.conf or pele.conf.

        Returns
        ----------
        - pele_steps : int
            Number of pele steps of the simulation(s).
        """

        path = str(pathlib.Path().absolute())
        adaptive_path = os.path.join(path, 'adaptive.conf')

        if os.path.isfile(adaptive_path) == True:

            #
            print('     -   Retrieving information from adaptive.conf.')
            #

            with open(adaptive_path) as filein:

                for line in filein:

                    if "peleSteps" in line:

                        peleSteps_string = line.split()[2]
                        pele_steps = int(peleSteps_string.split(',')[0])

        else:

            peleconf_path = os.path.join(path, 'pele.conf')

            #
            print('     -   Retrieving information from pele.conf.')
            #

            if os.path.isfile(peleconf_path) == True:

                with open(peleconf_path) as filein:

                    for line in filein:

                        if "numberOfPeleSteps" in line:

                            pele_steps = int(line.split()[-1])

            else:

                pele_steps = None

                #
                print('     -   No .conf was found.')
                print('     -   The step weighted scoring function will not \n'
                      '         be taken into account.')
                #

        return pele_steps

    def step_weighted(be,
                      step,
                      pele_steps):
        """
        Function
        ----------
        Calculates step weighted energy.

        Parameters
        ----------
        - be : list
            Binding energies of all the simulation.
        - step : list
            Steps associated to poses for all the simulation.
        - pele_steps : list
            Steps associated to poses for all the simulation.

        Returns
        ----------
        - ene_step : float
            Value of the step weighted energy.
        - num_step : numpy.array
            Array with all the steps of all the simulation.
        """

        num_steps = []

        for i in range(len(step) - 1):

            if step[i] == 0 and step[i+1] != 0:

                num_steps.append(step[i+1] - step[i])

            elif step[i] == 0 and step[i+1] == 0:

                num_steps.append(pele_steps)

            elif step[i] != 0 and step[i+1] == 0:

                num_steps.append(pele_steps-step[i])

            else:

                num_steps.append(step[i+1] - step[i])

        num_steps.append(pele_steps - step[len(step)-1])

        num_steps = np.array(num_steps)
        numerator = be.dot(num_steps)
        denominator = np.sum(num_steps)
        ene_step = numerator/denominator

        return ene_step, num_steps

    def simulation_evolution(files,
                             report_name,
                             pele_steps,
                             column,
                             folderpath):
        """
        Function
        ----------
        Calculates the evolution of the boltzmann weighted value  and the minimum value 
        of the metric chosen with column and stores data in two dictionaries.

        Parameters
        ----------
        - files : list
            The path to the directory where the output of the PELE simulation is located.
        - report_name : str
            Name of the reports to obtain the data from
        - pele_steps : int
            Number of pele steps in each epoch.
        - column : int 
            Column where the interesting data is located. 
        - folderpath : str
            Path where the different epochs of the simulation are located.

        Returns
        ----------
        - bz_dict : dict
            Dictionary with the acummulated boltzmann weighted values of the chosen metric 
            per epoch.
        - min_dict : dict
            Dictionary with the accumulated minium values of the chosen metric per epoch.
        - step_dict :dict
            Dictionary with the number of maximum accepted steps per epoch.
        """

        def step_data_reader(files,
                             report_name,
                             column,
                             folderpath,
                             step,
                             step_data):
            """
            Function
            ----------
            Opens all reports and extracts information from the column inputed by the user
            at a certain step.

            Parameters
            ----------
            - files : list
                The path to the directory where the output of the PELE simulation is located.
            - report_name : str
                Name of the reports to obtain the data from
            - column : int 
                Column where the interesting data is located. 
            - folderpath : str
                Path where the different epochs of the simulation are located.
            - step : int
                Accepted step we are retrieving data from.
            - step_data : list
                List with all the accumulated data of previous steps.

            Returns
            ----------
            - step : int
                Accepted step next iteration is going to retrieve data from.
            - step_data : list
                List with all the accumulated data of previous steps.
            """

            for file in files:

                cont = 0

                if file.startswith(report_name):

                    with open(os.path.join(folderpath, file)) as filein:

                        for line in filein:

                            if cont != 0:

                                line = line.split()

                                if step == int(line[2]):

                                    step_data.append(float(line[column - 1]))

                            cont += 1

                        else:

                            continue

            step += 1

            return step_data, step

        numeric_files = [s for s in files if s.isnumeric()]

        step_dict = {}
        bz_dict = {}
        min_dict = {}

        if len(numeric_files) != 0:

            #
            print(' -   Calculating...')
            #

            for document in numeric_files:

                step = 0

                step_data = []
                step_list = []
                minimum_energy_list = []
                bz_energy_list = []

                new_directory = os.path.join(folderpath, document)

                if os.path.isdir(new_directory) and document.isnumeric():

                    files = os.listdir(new_directory)

                    while step < pele_steps:

                        if step > 1 and bz_energy_list[-1] == bz_energy_list[-2]:

                            break

                        step_list.append(int(step))
                        step_data, step = step_data_reader(files,
                                                           report_name,
                                                           column,
                                                           new_directory,
                                                           step,
                                                           step_data)

                        ene_bz = boltzmann_weighted(np.array(step_data), T)
                        minimum_energy = min(np.array(step_data))

                        bz_energy_list.append(ene_bz)
                        minimum_energy_list.append(minimum_energy)

                step_dict[int(document)] = step_list
                bz_dict[int(document)] = bz_energy_list
                min_dict[int(document)] = minimum_energy_list

        else:

            step = 0

            minimum_energy_list = []
            bz_energy_list = []
            step_data = []
            step_list = []

            #
            print(' -   Calculating...')
            #

            while step < pele_steps:

                if step > 1 and bz_energy_list[-1] == bz_energy_list[-2]:

                    print(
                        ' -   Maximum number of acepted steps in any report:', step - 1)

                    break

                step_list.append(int(step))
                step_data, step = step_data_reader(files,
                                                   report_name,
                                                   column,
                                                   folderpath,
                                                   step,
                                                   step_data)

                ene_bz = boltzmann_weighted(np.array(step_data), T)
                minimum_energy = min(np.array(step_data))

                bz_energy_list.append(ene_bz)
                minimum_energy_list.append(minimum_energy)

            step_dict[0] = step_list
            bz_dict[0] = bz_energy_list
            min_dict[0] = minimum_energy_list

        return bz_dict, min_dict, step_dict

    def plotter(step_dict,
                bz_dict,
                min_dict,
                path,
                input_folder):
        """
        Function
        ----------
        Orders dictionaries by epochs and plots the data calculated to save the results
        in a newu directory /evolution.

        Parameters
        ----------
        - bz_dict : dict
            Dictionary with the acummulated boltzmann weighted values of the chosen metric 
            per epoch.
        - min_dict : dict
            Dictionary with the accumulated minium values of the chosen metric per epoch.
        - step_dict :dict
            Dictionary with the number of maximum accepted steps per epoch.
        - folderpath : str
            Path where the different epochs of the simulation are located.
        - path : str
            Path to the working directory.
        - input_folder : str
            Name of the output directory from which we have obtained the data 
            of the dictionaries.
        """

        path_plots = os.path.join(path, 'evolution')

        if os.path.exists(path_plots) == False:
            os.mkdir(path_plots)

        step_dict = dict(sorted(step_dict.items()))
        bz_dict = dict(sorted(bz_dict.items()))
        min_dict = dict(sorted(min_dict.items()))

        plt.figure(figsize=(7,3.5))

        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()

        for key in step_dict:

            ax1.plot(step_dict[key], bz_dict[key], label=key)
            ax1.legend(loc='best', fontsize=5)
            ax2.plot(step_dict[key], min_dict[key], label=key)
            ax2.legend(loc='best', fontsize=5)

        if len(step_dict[key]) > 100:

            left, bottom, width, height = [0.65, 0.3, 0.2, 0.2]
            ax1_1 = fig1.add_axes([left, bottom, width, height])
            ax1_1.plot(step_dict[key], bz_dict[key], label=key)
            ax1_1.set_xlim([0,20])

        fig1.savefig(os.path.join(
            path_plots, input_folder + '_bz_evolution.png'))
        fig2.savefig(os.path.join(
            path_plots, input_folder + '_min_evolution.png'))

    #
    print(' ')
    print('**************************************************************')
    print('*                           DiSC                             *')
    print('**************************************************************')
    print(' ')
    #

    path = str(pathlib.Path().absolute())
    folderpath = os.path.join(path, input_folder)

    files = os.listdir(folderpath)

    if action == 'all':

        be, step = reader(files,
                          folderpath,
                          report_name,
                          column)

        minimum_energy = min(be)
        be = np.array(be)
        be4 = be/4

        ene_bz = boltzmann_weighted(be, T)
        ene_bz4 = boltzmann_weighted(be4, T)
        ene_step_bz = boltzmann_weighted(be, T, step)

        if pele_steps == None:

            #
            print(' -   No information about pele_steps was given.')
            #

            pele_steps = pelesteps_retriever()

            if pele_steps == None:

                ene_step = 0.

            else:

                ene_step, _ = step_weighted(be, step, pele_steps)

            #
            print(' -   Pele steps:', pele_steps)
            #

        else:

            ene_step, _ = step_weighted(be, step, pele_steps)

        #
        print(' ')
        print(' RESULTS:')
        print(' -   Minimum Binding Energy:', minimum_energy)
        print(' -   Average Binding Energy:', np.average(be))
        print(' -   Boltzmann weighted Energy:', ene_bz)
        print(' -   Step weighted Energy:', ene_step)
        print(' -   Step-Boltzmann weighted Energy:', ene_step_bz)
        print(' -   Boltzmann weighted corrected Energy:', ene_bz4)
        print(' ')
        #

        with open('energy.csv', 'w') as fileout:
            fileout.writelines(
                'Minimum,Average,Boltzmann weighted,Step weighted,Step-Boltzmann weighted,Boltzmann weighted corrected\n'
                '' + str(minimum_energy) + ',' + str(np.average(be)) + ',' + str(ene_bz) +
                ',' + str(ene_step) + ',' + str(ene_step_bz) +
                ',' + str(ene_bz4) + '\n'
            )

    elif action == 'evolution':

        if pele_steps == None:

            #
            print(' -   No information about pele_steps was given.')
            #

            pele_steps = pelesteps_retriever()

            #
            print(' -   Pele steps:', pele_steps)
            print(' ')
            #

        bz_dict, min_dict, step_dict = simulation_evolution(files,
                                                            report_name,
                                                            pele_steps,
                                                            column,
                                                            folderpath)

        plotter(step_dict, bz_dict, min_dict, path, input_folder)

        #
        print(' -   Images have been stored in /evolution directory.')
        #

    else:

        raise Exception(
            'ActionError: The action you chose is not either all of evolution.')


def main(args):
    """
    Function
    ----------
    It runs statistics function. 

    Parameters
    ----------
    - args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    start_time = time.time()

    statistics(input_folder=args.input_folder,
               report_name=args.report_name,
               T=args.temperature,
               pele_steps=args.pele_steps,
               column=args.column,
               action=args.action)

    print(' ')
    print('                    --Duration of the execution--                   ')
    print('                      %s seconds' % (time.time() - start_time))
    print(' ')
    print('*******************************************************************')


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
