# -*- coding: utf-8 -*-
"""
This module is designed to perform bootstrapping to a PELE output.
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import os
import pathlib
import argparse

from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from seaborn import displot

# Constants
R = 1.985e-3
T = 298.


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
                        default='output', help="Name of the output directory where the simulation\
        is located.")
    parser.add_argument("-rn", "--report_name", type=str, dest="report_name",
                        default='report', help="Name of the report files used for the simulation.")
    parser.add_argument("-ns", "--number_of_samples", type=int, dest="number_of_samples",
                        default=10, help="Number of bootstrap datasets to generate.")
    parser.add_argument("-m", "--metric", type=str, dest="metric",
                        default='BindingEnergy', help="Name of the metric you are interested in.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def bootstrapping(input_folder,
                  report_name,
                  number_of_samples,
                  metric):
    """
    Function
    ----------
    Bootstrap function that performs all the calculations.

    Parameters
    ----------
    - input_folder : str
        Name of the output directory where the simulation
        is located.
    - report_name : str
        Name of the report files used for the simulation.
    - number_of_samples : int
        Number of bootstrap datasets to generate.
    - metric : str
        Name of the metric you are interested in.
    """

    def path_definer(input_folder):
        """
        Function
        ----------
        Defines the important paths that are going to be used throughout the simulation.

        Parameters
        ----------
        - input_folder : str
            Name of the folder where the output of the simulation is located.

        Returns
        ----------
        - path_output : str
            Path to the output folder of the simulation.
        - path_results : str 
            Path to the results folder of the bootstrap analysis.
        """

        path = str(pathlib.Path().absolute())
        path_output = os.path.join(path, input_folder)
        path_results = os.path.join(path, 'bootstrap')

        if os.path.isdir(path_results) is False:
            os.mkdir(path_results)

        return path_output, path_results

    def data_retriever(report_name,
                       path_output):
        """
        Function
        ----------
        Reads all the results of a simulation and stores them in a dataframe.

        Parameters
        ----------
        - report_name : str
            Name of the report files where the results are located
        - path_output : str
            Path to the output directory.

        Return
        ----------
        - original_data : pd.DataFrame
            Data frame with all the data from the reports.
        - snapshots : int
            Total number of snapshots in the simulation.
        """

        def metric_retirever(path_epoch):
            """
            Function
            ----------
            Assigns to each metric an index to store it in a dictionary.

            Parameters
            ----------
            - path_epoch : str
                Path to an epoch directory.

            Return
            ----------
            - metrics_dict : dict
                Dictionary with different metrics of the simulation and their index in the report.
            """

            report_paths = [os.path.join(path_epoch, report) for report in os.listdir(
                path_epoch) if report.startswith(report_name)]

            with open(report_paths[0]) as filein:

                cont = 0

                for line in filein:

                    if cont == 0:

                        metrics = line.split()[1:]

                    cont += 1

            metrics_dict = {k: v for k, v in enumerate(metrics)}

            return metrics_dict

        def epoch_retriever(report_name,
                            path_epoch,
                            metric_dict):
            """
            Function
            ----------
            Reads all the results of an epoch and stores them in a dataframe.

            Parameters
            ----------
            - input_folder : str
                Name of the output folder where the results are located
            - report_name : str
                Name of the report files where the results are located

            Return
            ----------
            - original_data : pd.DataFrame
                Data frame with all the data from the reports.
            """

            epoch_dict = {}
            report_paths = [os.path.join(path_epoch, report) for report in os.listdir(
                path_epoch) if report.startswith(report_name)]

            for report_path in report_paths:

                trajectory = report_path.split('/')[-1].split('_')[1]
                trajectory_dict = {}

                with open(report_path) as filein:

                    cont = 0

                    for line in filein:

                        model_dict = {}

                        if cont != 0:

                            line = line.split()[1:]
                            cont_metric = 0

                            for metric in line:

                                if metric_dict[cont_metric] == 'Step':

                                    model_dict[metric_dict[cont_metric]] = int(
                                        metric)

                                elif metric_dict[cont_metric] == 'numberOfAcceptedPeleSteps':

                                    model_dict[metric_dict[cont_metric]] = int(
                                        metric)

                                else:

                                    model_dict[metric_dict[cont_metric]] = float(
                                        metric)

                                cont_metric += 1

                        trajectory_dict[cont-1] = model_dict

                        cont += 1

                epoch_dict[trajectory] = trajectory_dict

            return epoch_dict

        def dict_to_dataframe(data_dict):
            """
            Function
            ----------
            Transforms a dictionary of dictionaries into a dataframe.

            Parameters
            ----------
            - data_dict : dict
                Dictionary with all the simulation data.

            Return
            ----------
            - original_df : pd.DataFrame
                Data frame with all the data from the reports.
            """

            original_df = pd.DataFrame([(epoch, trajectory, model, met, value)
                                        for epoch, traj_mod_met_val in data_dict.items()
                                        for trajectory, mod_met_val in traj_mod_met_val.items()
                                        for model, met_val in mod_met_val.items()
                                        for met, value in met_val.items()])

            original_df.columns = ['epoch', 'trajectory', 'model',
                                   'metric', 'value']

            original_df['epoch'] = original_df['epoch'].astype(int)
            original_df['trajectory'] = original_df['trajectory'].astype(int)
            original_df['model'] = original_df['model'].astype(int)
            original_df['value'] = original_df['value'].astype(float)

            original_df = original_df.sort_values(
                ['epoch', 'trajectory', 'model'], ascending=[True, True, True])
            original_df = original_df.reset_index()
            original_df = original_df.drop('index', 1)

            return original_df

        data_dict = {}
        epoch_paths = [os.path.join(path_output, epoch) for epoch in os.listdir(
            path_output) if epoch.isnumeric()]

        if len(epoch_paths) != 0:

            metrics_dict = metric_retirever(epoch_paths[0])

            for epoch_path in epoch_paths:

                epoch = epoch_path.split('/')[-1]
                data_dict[epoch] = epoch_retriever(report_name,
                                                   epoch_path,
                                                   metrics_dict)

        else:

            metrics_dict = metric_retirever(path_output)
            data_dict[0] = epoch_retriever(report_name,
                                           path_output,
                                           metrics_dict)

        original_df = dict_to_dataframe(data_dict)
        snapshots = int(len(original_df)/6)

        return original_df, snapshots

    def bootstrap_dataset(original_df,
                          snapshots):
        """
        Function
        ----------
        Generates a bootstrap datset from the data retrieved and random numbers.

        Parameters
        ----------
        - original_df : pd.DataFrame
            Data frame with all the data from the reports.
        - snapshots : int
            Number of total snapshots in a simulation.

        Return
        ----------
        - bootstrap_df : pd.DataFrame
            Data frame with all the data of the bootstrap dataset.
        """

        def random_numbers_generator(snapshots):
            """
            Function
            ----------
            Generates a random array of length: snapshots, of numbers
            ranging from 0 to: snapshots.

            Parameters
            ----------
            - snapshots : int
                Number of snaphots in the simulation.

            Return
            ----------
            - random_numbers : np.array
                Array of length: snapshots, of numbers ranging [0,snapshots].
            """

            seed = int(''.join(c for c in str(
                datetime.now().time()) if c.isdigit())[-5:])
            np.random.seed(seed)
            random_numbers = np.random.randint(0, snapshots, snapshots)

            return random_numbers

        def random_to_data(original_df,
                           random_numbers):
            """
            Function
            ----------
            Generates a dataframe from the original dataframe and the 
            random values of the array.

            Parameters
            ----------
            - original_df : pd.DataFrame
                Data frame with all the data from the reports.
            - random_numbers : np.array
                Array of length: snapshots, of numbers ranging [0,snapshots].

            Return
            ----------
            - bootstrap_df : pd.DataFrame
                Data frame with all the data from original dataframe and the 
                random numbers.
            """

            indexes_list = [[number*6, number*6 + 1, number*6 + 2, number *
                             6 + 3, number*6 + 4, number*6 + 5] for number in random_numbers]
            indexes = [item for sublist in sorted(
                indexes_list) for item in sublist]
            bootstrap_df = original_df.iloc[indexes]

            return bootstrap_df

        random_numbers = random_numbers_generator(snapshots)
        bootstrap_df = random_to_data(original_df,
                                      random_numbers)

        return bootstrap_df

    def metric_calculator(metric,
                          bootstrap_df):
        """
        Function
        ----------
        Calculates all the scores for a given dataframe and metric.

        Parameters
        ----------
        - metric : str
            Metric the user is interested in.
        - bootstrap_df : pd.DataFrame
            Data frame with all the data from original dataframe and the 
            random numbers.

        Return
        ----------
        - average : float
            The dataframe's minimum value for the specific metric.
        - maximum : float
            The dataframe's maximum value for the specific metric.
        - average : float
            The dataframe's average value for the specific metric.
        - boltzmann : float
            The dataframe's boltzmann average value for the specific metric.
        """

        def from_dataframe_to_vector(metric,
                                     bootstrap_df):
            """
            Function
            ----------
            From the dataframe obtains the values of the interesting metric and 
            the total energy.

            Parameters
            ----------
            - metric : str
                Metric the user is interested in.
            - bootstrap_df : pd.DataFrame
                Data frame with all the data from original dataframe and the 
                random numbers.

            Return
            ----------
            - metric_vector: np.array
                Array with all the values of the dataframe corresponding to the 
                metric of interest.
            - total_ene_vector : np.array
                Array with all the values of the dataframe corresponding to the 
                total energy.
            """

            metric_vector = bootstrap_df[bootstrap_df['metric']
                                         == metric]['value'].to_numpy()
            total_ene_vector = bootstrap_df[bootstrap_df['metric']
                                            == 'currentEnergy']['value'].to_numpy()

            return metric_vector, total_ene_vector

        def minimum_score(metric_vector):
            """
            Function
            ----------
            Compute minimum of a vector.

            Parameters
            ----------
            - metric_vector: np.array
                Array with all the values of the dataframe corresponding to the 
                metric of interest.

            Return
            ----------
            - minimum : float
                Minimum value of the inputed vector.
            """

            minimum = min(list(metric_vector))

            return minimum

        def maximum_score(metric_vector):
            """
            Function
            ----------
            Compute maximum of a vector.

            Parameters
            ----------
            - metric_vector: np.array
                Array with all the values of the dataframe corresponding to the 
                metric of interest.

            Return
            ----------
            - maximum : float
                Maximum value of the inputed vector.
            """

            maximum = max(list(metric_vector))

            return maximum

        def average_score(metric_vector):
            """
            Function
            ----------
            Compute average of a vector.

            Parameters
            ----------
            - metric_vector: np.array
                Array with all the values of the dataframe corresponding to the 
                metric of interest.

            Return
            ----------
            - average : float
                Average value of the inputed vector.
            """

            average = np.average(metric_vector)

            return average

        def boltzmann_score(metric_vector, total_ene_vector, T):
            """
            Function
            ----------
            Compute boltzmann average of a vector.

            Parameters
            ----------
            - metric_vector : np.array
                Array with all the values of the dataframe corresponding to the 
                metric of interest.
            - T : float
                Temperature to perform the Boltzmann weights with.
            - total_ene_vector : np.array
                Steps associated to poses for all the simulation.

            Returns
            ----------
            - boltzmann : float
                Boltzmann average value of the inputed vectors.
            """

            total_ene_vector -= min(list(total_ene_vector))

            exp_bz = np.exp(-total_ene_vector/(R*T))
            nominator = metric_vector.dot(exp_bz)
            denominator = np.sum(exp_bz)
            boltzmann = nominator/denominator

            return boltzmann

        metric_vector, total_ene_vector = from_dataframe_to_vector(
            metric, bootstrap_df)

        minimum = minimum_score(metric_vector)
        maximum = maximum_score(metric_vector)
        average = average_score(metric_vector)
        boltzmann = boltzmann_score(metric_vector, total_ene_vector, T)

        return minimum, maximum, average, boltzmann

    def iterator(number_of_samples,
                 metric,
                 original_df,
                 snapshots):
        """
        Function
        ----------
        Iterates number_of_samples times to obtain the different scores for
        the specified metric.

        Parameters
        ----------
        - number_of_samples : int
            Number of times we want to iterate the process of creating
            a bootsrap dataset and calculate scores.
        - metric : str
            Metric the user is interested in.
        - original_df : pd.DataFrame
            Data frame with all the data from the reports.
        - snapshots : int
            Number of total snapshots in a simulation.

        Return
        ----------
        - average : float
            The dataframe's minimum value for the specific metric.
        - maximum : float
            The dataframe's maximum value for the specific metric.
        - average : float
            The dataframe's average value for the specific metric.
        - boltzmann : float
            The dataframe's boltzmann average value for the specific metric.
        """

        minimum_list = []
        maximum_list = []
        average_list = []
        boltzmann_list = []

        for i in range(number_of_samples):

            bootstrap_df = bootstrap_dataset(original_df, snapshots)
            minimum, maximum, average, boltzmann = metric_calculator(
                metric, bootstrap_df)

            minimum_list.append(minimum)
            maximum_list.append(maximum)
            average_list.append(average)
            boltzmann_list.append(boltzmann)

        return minimum_list, maximum_list, average_list, boltzmann_list

    def statistics(metric_list):
        """
        Function
        ----------
        Calculates average and error of the inputed list.

        Parameters
        ----------
        - metric_list : list
            List of values we want to know the average and the standard
            deviation of.

        Return
        ----------
        - average : float
            Avergage of the inputed list.
        - standard_dev : float
            Standard deviation of the inputed list.
        """

        vector = np.array(metric_list)
        average = np.average(vector)
        standard_dev = np.std(vector)

        return average, standard_dev

    def plotter(path_results,
                data,
                average,
                error,
                scoring):

        """
        Function
        ----------
        Plots results obtained with bootstrap.

        Parameters
        ----------
        - path_results : str 
            Path to the results folder of the bootstrap analysis.
        - data : list 
            List of data points calculated from bootstrap we want to plot.
        - average : float
            Avergage of the inputed list.
        - error : float
            Standard deviation of the inputed list.
        - scoring : string 
            Name of the scoring function used.
        """

        displot(data, kind="kde", color='black', label='KDE plot')
        plt.title(scoring + ' Distribution')
        plt.axvline(x=average, color='red', label='Average = ' + str("{:.3f}".format(average)))
        plt.axvline(x=average + error, color='green', label='Error = ' + str("{:.3f}".format(error)))
        plt.axvline(x=average - error, color='green')
        plt.legend(loc="best")
        plt.xlabel(scoring)
        plt.ylabel('Density')
        plt.tight_layout()
        plt.savefig(os.path.join(path_results, scoring + '_distribution.png'), format='png', bbox_inches="tight")

    #
    print(' ')
    print('*******************************************************************')
    print('*                        peleBootstrap                            *')
    print('*******************************************************************')
    print(' ')
    #

    path_output, path_results = path_definer(input_folder)

    #
    print(' -   Retrieving data')
    #

    original_df, snapshots = data_retriever(report_name, path_output)

    #
    print(' -   Generating ' + str(number_of_samples) + ' datasets...')
    #

    minimum_list,\
        maximum_list,\
        average_list,\
        boltzmann_list = iterator(number_of_samples,
                                  metric,
                                  original_df,
                                  snapshots)

    average_minimum, standard_dev_minimum = statistics(minimum_list)
    average_maximum, standard_dev_maximum = statistics(maximum_list)
    average_average, standard_dev_average = statistics(average_list)
    average_boltzmann, standard_dev_boltzmann = statistics(boltzmann_list)

    data = {
        'minimum': [average_minimum, standard_dev_minimum],
        'maximum': [average_maximum, standard_dev_maximum],
        'average': [average_average, standard_dev_average],
        'boltzmann': [average_boltzmann, standard_dev_boltzmann],
    }

    #
    print(' -   Writing files\n')
    #

    results_df = pd.DataFrame(data, index=['average', 'error'])
    results_df.to_csv(os.path.join(path_results, 'results.csv'))

    plotter(path_results,
            minimum_list,
            average_minimum,
            standard_dev_minimum,
            'Minimum')

    plotter(path_results,
            maximum_list,
            average_maximum,
            standard_dev_maximum,
            'Maximum')

    plotter(path_results,
            average_list,
            average_average,
            standard_dev_average,
            'Average')

    plotter(path_results,
            boltzmann_list,
            average_boltzmann,
            standard_dev_boltzmann,
            'Boltzmann')


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

    bootstrapping(input_folder=args.input_folder,
                  report_name=args.report_name,
                  number_of_samples=args.number_of_samples,
                  metric=args.metric)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
