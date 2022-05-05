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

    # Analysis flags
    parser.add_argument("--bootstrap", dest='bootstrap_bool', default=False, action='store_true',
                        help="Flag to choose if bootstrap analysis is required.")
    parser.add_argument("--filter", dest='filter_bool', default=False, action='store_true',
                        help="Flag to choose if filter analysis is required.")
    parser.add_argument("--analyze", dest='analyzer_bool', default=False, action='store_true',
                        help="Flag to choose if analysis is required.")


    # General flags
    parser.add_argument("-d", "--directory", type=str, dest="input_folder",
                        default='output', help="Name of the output directory where the simulation\
        is located.")
    parser.add_argument("-rn", "--report_name", type=str, dest="report_name",
                        default='report', help="Name of the report files used for the simulation.")

    # Filter flags
    parser.add_argument("-es", "--equilibration_steps", type=int, dest="equilibration_steps",
                        default=None, help="Number of steps from first reports we want to omit.")
    parser.add_argument("-mt", "--metric_threshold", type=str, dest="metric_threshold",
                        default=None, help="List of [metric,[min,max]] where metric is str and min and max, float")
    parser.add_argument("-q", "--quantile", type=str, dest="quantile_flag",
                        default=None, help="Quantile we are interested in. List of [metric,value] where metric is\
                         str and values is float")

    # Bootstrap flags
    parser.add_argument("-ns", "--number_of_samples", type=int, dest="number_of_samples",
                        default=10, help="Number of bootstrap datasets to generate.")

    # Bootstrap and Analyzer flag                   
    parser.add_argument("-m", "--metric", type=str, dest="metric",
                        default='BindingEnergy', help="Name of the metric you are interested in.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def data(input_folder,
         report_name):
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
        Defines the important paths that are going to be used in this function.

        Parameters
        ----------
        - input_folder : str
            Name of the folder where the output of the simulation is located.

        Returns
        ----------
        - path_output : str
            Path to the output folder of the simulation.
        """

        path = str(pathlib.Path().absolute())
        path_output = os.path.join(path, input_folder)

        return path_output

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

                trajectory = report_path.split('/')[-1].split('_')[-1]
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

    #
    print(' ')
    print('*******************************************************************')
    print('*                        peleAnalysis                             *')
    print('*******************************************************************')
    print(' ')
    #

    path_output = path_definer(input_folder)

    #
    print(' *   Retrieving data')
    #

    original_df, snapshots = data_retriever(report_name, path_output)

    return original_df, snapshots


def filter(original_df,
           equilibration_steps,
           metric_threshold,
           quantile_flag):

    def equilibration_steps_remover(original_df,
                                    equilibration_steps):
        """
        Function
        ----------
        Obtain a data frame without equilibration steps.

        Parameters
        ----------
        - original_df : pd.DataFrame
            Data frame with reports' information.
        - equilibration_steps : int
            Number of steps we want to omit.

        Returns
        ----------
        - original_df : list
            Data frame with equilibration steps omited filtered.
        """

        if equilibration_steps > 5:
            raise Exception(
                'TooManyEquilibrationSteps: Maximum equilibration steps are 5.')

        print('     -   Deleting ' +
              str(equilibration_steps) + ' equilibration steps')
        equilibration_df = original_df.drop(
            original_df[original_df.epoch == 0][original_df[original_df.epoch == 0].model <= equilibration_steps].index)

        return equilibration_df

    def quantile_retriever(original_df,
                           quantile_flag):
        """
        Function
        ----------
        Obtain reports' metrics for reports with metric below a certain quantile.

        Parameters
        ----------
        - original_df : pd.DataFrame
            Data frame with reports' information.
        - quantile_flag : str
            Flag with metric and quantile information. Syntax: [metric,quantile]

        Returns
        ----------
        - quantile_df : list
            Data frame with quantile filtered.
        """

        quantile_flag = quantile_flag.strip('][').split(',')
        modified_df = original_df.loc[original_df['metric']
                                      == quantile_flag[0]]
        quantile_value = modified_df['value'].quantile(float(quantile_flag[1]))
        deletion_df = modified_df[modified_df['value'] > quantile_value]

        print('     -   Deleting all values above ' + str(quantile_value)
              + ' of the ' + quantile_flag[0] + '.')

        for _, row in deletion_df.iterrows():

            indexes = original_df[
                (original_df['epoch'] == row[0]) &
                (original_df['trajectory'] == row[1]) &
                (original_df['model'] == row[2])].index

            original_df = original_df.drop(indexes)

        quantile_df = original_df

        return quantile_df

    def metric_remover(original_df,
                       metric_threshold):
        """
        Function
        ----------
        Obtain reports' metrics for reports with metric withina range of a metric.

        Parameters
        ----------
        - original_df : pd.DataFrame
            Data frame with reports' information.
        - metric_threshold : str
            List with metric we want to apply the metric and threshold on and value 
            of quantile.

        Returns
        ----------
        - original_df : list
            Data frame with metric filtered.
        """

        metric_flag = metric_threshold.strip('][').split(',')
        modified_df = original_df.loc[original_df['metric']
                                      == metric_flag[0]]
        modified_df_1 = modified_df[modified_df['value'] > float(metric_flag[2])]
        modified_df_2 = modified_df[modified_df['value'] <= float(metric_flag[1])]
        deletion_df = pd.concat([modified_df_1, modified_df_2])

        print('     -   Deleting all values not belonging to [' \
            + metric_flag[1] + ',' + metric_flag[2] + ') of the ' + metric_flag[0] + '.')

        for _, row in deletion_df.iterrows():

            indexes = original_df[
                (original_df['epoch'] == row[0]) &
                (original_df['trajectory'] == row[1]) &
                (original_df['model'] == row[2])].index

            original_df = original_df.drop(indexes)

        metric_df = original_df

        return metric_df

    print('\n *   Filtering data.')

    cont = 0

    if equilibration_steps is None:

        cont += 1
        pass

    else:
        
        equilibration_df = equilibration_steps_remover(original_df,
                                                       equilibration_steps)

    if quantile_flag is not None and metric_threshold is not None:

        quantile_df = quantile_retriever(equilibration_df,
                                         quantile_flag)

        output_df = metric_remover(quantile_df,
                                   metric_threshold)
    
    elif quantile_flag is not None and metric_threshold is None:

        output_df = quantile_retriever(equilibration_df,
                                       quantile_flag)

    elif quantile_flag is None and metric_threshold is not None:

        output_df = metric_remover(equilibration_df,
                                   metric_threshold)

    else:

        cont += 1
        pass

    if cont == 2:

        raise Exception('NoCorrectionImplemented: The --filter flag has been written but no condition to filter\
        has been chosen.')

    return output_df


def bootstrapping(number_of_samples,
                  metric,
                  original_df,
                  snapshots):
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

    def path_definer():
        """
        Function
        ----------
        Defines the important paths that are going to be used in this function.

        Parameters
        ----------

        Returns
        ----------
        - path_results : str 
            Path to the results folder of the analysis.
        """

        path = str(pathlib.Path().absolute())
        path_results = os.path.join(path, 'analyzer','bootstrap')

        if os.path.isdir(path_results) is False:
            os.mkdir(path_results)

        return path_results

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
        - minimum : float
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

        from seaborn import displot

        displot(data, kind="kde", color='black', label='KDE plot')
        plt.title(scoring + ' Distribution')
        plt.axvline(x=average, color='red', label='Average = ' +
                    str("{:.3f}".format(average)))
        plt.axvline(x=average + error, color='green',
                    label='Error = ' + str("{:.3f}".format(error)))
        plt.axvline(x=average - error, color='green')
        plt.legend(loc="best")
        plt.xlabel(scoring)
        plt.ylabel('Density')
        plt.tight_layout()
        plt.savefig(os.path.join(path_results, scoring +
                    '_distribution.png'), format='png', bbox_inches="tight")

    path_results = path_definer()

    #
    print('\n *   Bootstrapping')
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
    print(' -   Writing files')
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


def dataframe_trimming(output_df):
    """
    Function
    ----------
    Trims non-important columns out of the dataframe and reshapes it.
    Parameters
    ----------
    output_df : pd.DataFrame
        Dataframe with all the information once it has been filtered.
    Returns
    ----------
    - output_df : pd.DataFrame 
        Dataframe without columns: Steps and nonAcceptedSteps
    """
    output_df = output_df[['metric','value']]
    output_df = output_df.set_index('metric').stack().reset_index(level=1, drop=True).to_frame()
    output_df['new_col'] = output_df.groupby(level='metric').cumcount()
    output_df = output_df.pivot(columns='new_col', values=0)
    output_df = output_df.reset_index().rename_axis(None, axis=1)
    output_df = output_df.T
    output_df.columns = output_df.iloc[0] 
    output_df = output_df[1:]
    output_df.head()
    output_df.drop('Step', inplace=True, axis=1)
    output_df.drop('numberOfAcceptedPeleSteps', inplace=True, axis=1)

    return output_df


def plot_function(output_df):
    """
    Function
    ----------
    Plot all the metrics vs all the metrics.

    Parameters
    ----------
    output_df : pd.DataFrame
        Dataframe with all the information once it has been filtered.
    """

    def path_definer():
        """
        Function
        ----------
        Defines the important paths that are going to be used in this function.

        Returns
        ----------
        - path_results : str 
            Path to the results folder of the analysis.
        """

        path = str(pathlib.Path().absolute())
        path_results = os.path.join(path, 'analyzer', 'images')

        if os.path.isdir(path_results) is False:
            os.mkdir(path_results)

        return path_results

    path_results = path_definer()

    metrics = list(output_df.columns.values.tolist())

    for metric_x in metrics:
        for metric_y in metrics:
            if metrics.index(metric_y) != metrics.index(metric_x):

                x = np.array(output_df[metric_x])
                y = np.array(output_df[metric_y])

                plt.scatter(x,y,marker = 'o')
                plt.xlabel(str(metric_x))
                plt.ylabel(str(metric_y))
                plt.margins()
                plt.title(str(metric_y) + ' vs ' + str(metric_x))
                plt.savefig(os.path.join(path_results,metric_x + '_' + metric_y + '.png'))
                plt.close()


def analyzer(df,
             metric_analyze):
    """
    Function
    ----------
    Calculate different scores fot the metric we want to analyze.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with all the information once it has been filtered.
    metric_analyze : str
        Metric we want to obtain information of.
    """
    def path_definer():
        """
        Function
        ----------
        Defines the important paths that are going to be used in this function.

        Parameters
        ----------

        Returns
        ----------
        - path_results : str 
            Path to the results folder of the analysis.
        """

        path = str(pathlib.Path().absolute())
        path_results = os.path.join(path, 'analyzer')

        if os.path.isdir(path_results) is False:
            os.mkdir(path_results)

        return path_results

    def boltzmann_weighted(be,
                           te,
                           T):
        """
        Function
        ----------
        Calculates boltzmann weighted energy.

        Parameters
        ----------
        - be : list
            Binding energies of all the simulation.
        - te : list
            Total energies of all the simulation.
        - T : float
            Temperature to perform the Boltzmann weights with.
        - steps : list
            Steps associated to poses for all the simulation.

        Returns
        ----------
        - ene_bz : float
            Value of the boltzmann weighted energy.
        """

        exp_bz = np.exp(-te/(R*T))
        nominator = be.dot(exp_bz)
        denominator = np.sum(exp_bz)
        ene_bz = nominator/denominator

        return ene_bz

    print('\n *   Data analyzer.')

    path_results = path_definer()
    vector = np.array(df[metric_analyze]).astype(float)

    te = np.array(df['currentEnergy']).astype(float)
    min_energy = np.min(te)
    te = np.array(te) - min_energy

    min = np.min(vector)
    ave = np.average(vector)
    bz = boltzmann_weighted(vector, te, T)

    with open(os.path.join(path_results,'energy.csv'), 'w') as fileout:
        fileout.writelines(
            'Minimum,Average,Boltzmann weighted\n'
            '' + str(min) + ',' +
            str(ave) + ',' + str(bz) + '\n'
        )

    return min, ave, bz


def ensambler(bootstrap_bool,
              filter_bool,
              analyzer_bool,
              input_folder,
              report_name,
              equilibration_steps,
              metric_threshold,
              quantile_flag,
              number_of_samples,
              metric):
    """
    Function
    ----------
    Function that joins all the other functions.

    Parameters
    ----------
    bootstrap_bool : bool
        Boolean to know whether bootstrap analysis is wanted.
    filter_bool : bool
        Boolean to know whether filter analysis is wanted.
    input_folder : str
        Name of the output directory where the simulation is located.
    report_name : str
        Name of the report files used for the simulation.
    equilibration_steps : int
        Number of steps from first reports we want to omit.
    metric_threshold : str
        List of [metric,[min,max]] where metric is str and min and max, float.
    quantile_flag : str
        Quantile we are interested in. List of [metric,value] where metric is str and values is float
    number_of_samples : int
        Number of bootstrap datasets to generate.
    metric : str
        Name of the metric you are interested in both for analyzer and bootstrap.
    """

    original_df, snapshots = data(input_folder,
                                  report_name)

    cont = 0

    if bootstrap_bool and not filter_bool:

        bootstrapping(number_of_samples,
                      metric,
                      original_df,
                      snapshots)

    elif not bootstrap_bool and filter_bool:

        output_df = filter(original_df,
                           equilibration_steps,
                           metric_threshold,
                           quantile_flag)

    elif bootstrap_bool and filter_bool:

        bootstrapping(number_of_samples,
                      metric,
                      original_df,
                      snapshots)

        output_df = filter(original_df,
                           equilibration_steps,
                           metric_threshold,
                           quantile_flag)
    
    else: 

        cont += 1

    output_df = dataframe_trimming(output_df)

    if analyzer_bool:

        analyzer(output_df,
                 metric)

    else:

        cont += 1

        if cont == 2:

            raise Exception('NoActionError: No action has been chosen: --bootstrap and/or --filter and/or --analyze')

    plot_function(output_df)

    print('\n -   All results stored inside analyzer.\n')


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

    ensambler(args.bootstrap_bool,
              args.filter_bool,
              args.analyzer_bool,
              args.input_folder,
              args.report_name,
              args.equilibration_steps,
              args.metric_threshold,
              args.quantile_flag,
              args.number_of_samples,
              args.metric)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
