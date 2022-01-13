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
                        default='output', help="Name of the directory where the files to\
        analyze are located.")
    parser.add_argument("-rn", "--report_name", type=str, dest="report_name",
                        default='report_', help="Name of the report files.")
    parser.add_argument("-T", "--temperature", type=float, dest="temperature",
                        default=298., help="Temperature of the experiment.")
    parser.add_argument("-pS", "--pele_Steps", type=int, dest="pele_steps",
                        default=25, help="Number of Pele Steps in the simulation.")
    parser.add_argument("-c", "--column", type=int, dest="column",
                        default=5, help="Column of the report where the interesting\
                             metric is located.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def statistics(input_folder,
               report_name,
               T,
               pele_steps,
               column):
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
    """

    def reader(files, folderpath, column):
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

            for file in files:

                new_directory = folderpath

                cont = 0

                # Checking whether path exists and the file begins with the report_name selected
                if file.startswith(report_name):

                    file_path = os.path.join(new_directory, file)

                    with open(file_path, 'r') as filein:

                        for line in filein:
                            if cont != 0:
                                line = line.split('   ')
                                be.append(float(line[column-1]))
                                step.append(int(line[1]))
                            cont += 1
            return be, step

        be = []
        step = []
        numeric_files = [s for s in files if s.isnumeric()]

        if len(numeric_files) != 0:

            for document in numeric_files:

                # Checking if the folder exist
                new_directory = os.path.join(folderpath, document)

                if os.path.isdir(new_directory) and document.isnumeric():

                    # Listing files inside
                    files = os.listdir(new_directory)

                    if report_name in files == False:
                        raise Exception('FilePathError: There is no file containing ' + report_name + ' in it. \
                        Please check the path to the files and the files name.')

                    be, step = file_reader(
                        files, new_directory, report_name, column)

        else:

            if report_name in files == False:
                raise Exception('FilePathError: There is no file containing ' + report_name + ' in it. \
                Please check the path to the files and the files name.')

            be, step = file_reader(files, folderpath, report_name, column)

        return be, step

    def boltzmann_weighted(be, T, steps=[]):
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

    def step_weighted(be, step, pele_steps):
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

    path = str(pathlib.Path().absolute())
    folderpath = os.path.join(path, input_folder)

    if os.path.isdir(folderpath) == False:
        raise Exception(
            'FolderPathError: There is no folder with this name. Please check the path and the folder name.')

    files = os.listdir(folderpath)

    be, step = reader(files, folderpath, column)

    minimum_energy = min(be)
    be = np.array(be)
    be4 = be/4

    ene_bz = boltzmann_weighted(be, T)
    ene_bz4 = boltzmann_weighted(be4, T)
    ene_step_bz = boltzmann_weighted(be, T, step)
    ene_step, _ = step_weighted(be, step, pele_steps)

    #
    print(' ')
    print('**************************************************************')
    print('*                        ENERGIES                            *')
    print('**************************************************************')
    print(' ')
    print('- Minimum Binding Energy:', minimum_energy)
    print('- Average Binding Energy:', np.average(be))
    print('- Boltzmann weighted Energy:', ene_bz)
    print('- Step weighted Energy:', ene_step)
    print('- Step-Boltzmann weighted Energy:', ene_step_bz)
    print('- Boltzmann weighted corrected Energy:', ene_bz4)
    print(' ')
    #

    with open('energy.csv', 'w') as fileout:
        fileout.writelines(
            'Minimum,Average,Boltzmann weighted,Step weighted,Step-Boltzmann weighted,Boltzmann weighted corrected\n'
            '' + str(minimum_energy) + ',' + str(np.average(be)) + ',' + str(ene_bz) +
            ',' + str(ene_step) + ',' + str(ene_step_bz) +
            ',' + str(ene_bz4) + '\n'
        )


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
               column=args.column)

    print('                    --Duration of the execution--                   ')
    print('                      %s seconds' % (time.time() - start_time))
    print(' ')
    print('*******************************************************************')


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
