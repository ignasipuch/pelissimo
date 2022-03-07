# -*- coding: utf-8 -*-

# Imports
import sys
import os
import pathlib
import argparse
import numpy as np
import time
from scipy import stats
import matplotlib.pyplot as plt

# Constant
R = 1.985e-3

DG = {'1CB0': -7.6,
      '1K27': -12.3,
      '1HPV': -13.3,
      '1HSG': -13.1,
      '1MSM': -14.9,
      '1T3R': -15.0,
      '6QGE': -7.8,
      '6QGF': -8.3}


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
    parser.add_argument("-c", "--column", type=int, dest="column",
                        default=5, help="Column of the report where the interesting\
                             metric is located.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def statistics(input_folder,
               report_name,
               T,
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
                file_path = os.path.join(new_directory, file)

                # Checking whether path exists and the file has the report_name selected
                if os.path.isfile(file_path) and report_name in file:

                    with open(file_path, 'r') as filein:

                        for line in filein:
                            if cont != 0:
                                line = line.split('   ')
                                be.append(float(line[column-1]))
                                ene_t.append(float(line[3]))
                            cont += 1
            return be, ene_t

        be = []
        ene_t = []
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

                    be, ene_t = file_reader(
                        files, new_directory, report_name, column)

        else:

            if report_name in files == False:
                raise Exception('FilePathError: There is no file containing ' + report_name + ' in it. \
                Please check the path to the files and the files name.')

            be, ene_t = file_reader(files, folderpath, report_name, column)

        min_energy = min(ene_t)
        ene_t = np.array(ene_t) - min_energy
        be = np.array(be)

        return be, ene_t

    def boltzmann_weighted(be, ene_t, T):
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

        exp_bz = np.exp(-ene_t/(R*T))
        nominator = be.dot(exp_bz)
        denominator = np.sum(exp_bz)
        ene_bz = nominator/denominator

        return ene_bz

    def dict_to_list(dict_bz, dict_dG):
        """
        Function
        ----------
        Reads data from dictionaries and stores it in two lists.

        Parameters
        ----------
        - dict_bz : dict
            Dictionary containing the Boltzmann weighted energies
            calculated in the simulation.
        - dict_exp : dict
            Dictionary containing the experimental energies of the systems.

        Return
        ----------
        - dG_exp : list
            List with ordered results of experimental energies
        - dG_bz : list
            List with ordered results of Boltzmann weighted energies
        """

        dG_exp = []
        dG_bz = []

        for key in dict_bz:

            dG_exp.append(dict_dG[key])
            dG_bz.append(dict_bz[key])

        return dG_exp, dG_bz

    def correlation(dG_exp, dG_bz):
        """
        Function
        ----------
        Calculates the correlation between the calculated and the experimental
        data.

        Parameters
        ----------
        - dG_exp : list
            List with ordered results of experimental energies
        - dG_bz : list
            List with ordered results of Boltzmann weighted energies

        Returns
        ----------
        - r : float
            Correlation between experimental and calculated data.
        """

        _, _, r, _, _ = stats.linregress(dG_exp, dG_bz)

        return r

    path = str(pathlib.Path().absolute())
    outputspath = os.path.join(path, input_folder)

    if os.path.isdir(outputspath) == False:
        raise Exception(
            'FolderPathError: There is no folder with this name. Please check the path and the folder name.')

    outputpaths = os.listdir(outputspath)

    dict = {}

    for folder in outputpaths:

        system = folder.split('_')[4]
        folderpath = os.path.join(outputspath, folder)
        files = os.listdir(folderpath)

        be, ene_t = reader(files, folderpath, column)
        be = np.array(be, dtype=np.float128)
        ene_bz = boltzmann_weighted(be, ene_t, T)

        dict[system] = float(ene_bz)

    #
    print(dict)
    print(DG)
    #

    dG_exp, dG_bz = dict_to_list(dict,DG)
    r = correlation(dG_exp, dG_bz)

    return r


def plotter(kT, r):

    plt.title('kT vs r')
    plt.xlabel('kT (kcal/mol)')
    plt.ylabel('r')
    plt.xscale('log')
    plt.axvline(x=R*298., color='g')
    plt.plot(kT, r)
    plt.savefig(str(pathlib.Path().absolute())+'/energyplot.png')


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

    T1 = np.linspace(2, 10, 40, endpoint=False)
    T2 = np.linspace(10, 100, 20)
    T_list = list(np.concatenate((T1, T2))/R)

    kT = []
    r = []

    for T in T_list:

        r_val = statistics(input_folder=args.input_folder,
                           report_name=args.report_name,
                           T=T,
                           column=args.column)

        kT.append(R*T)
        r.append(r_val)

        print(' ')
        print('kT:', R*T)
        print('r:', r_val)
        print('R:', r_val**2)
        print(' ')
        print(str((100.*T_list.index(T) + 1)/60) + ' %')
        print('---------------------------------------------')

    plotter(kT, r)

    print(' ')
    print('                    --Duration of the execution--                   ')
    print('                      %s seconds' % (time.time() - start_time))
    print(' ')
    print('*******************************************************************')


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
