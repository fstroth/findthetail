"""
.. module:: ftt
:platform: Unix, Windows
:synopsis: Module for Paper XY.
.. moduleauthor:: Frederik Strothmann <frstrothmann@gmail.com>

"""

import os
from multiprocessing import Pool

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import genpareto
from jinja2 import Environment, FileSystemLoader, select_autoescape

from .teststatistics import au2, cramer_von_mises, anderson_darling


class Ftt:
    def __init__(self, data, data_name, mc_steps=1000, threads=1):
        self.data = data
        self.data_name = data_name

        self.mc_steps = mc_steps
        self.mc_steps_run = 0
        self.mc_error = None
        self.mc_counter_au2 = 0
        self.mc_counter_a2 = 0
        self.mc_counter_w2 = 0
        self.p_value_au2 = None
        self.p_value_a2 = None
        self.p_value_w2 = None
        self.q = None
        self.cond = None
        self.significant_digit_of_data = None

        self.au_2_data = None
        self.cramer_data = None
        self.anderson_data = None

        self.optimal_tail_index = None
        self.optimal_tail = None
        self.rv_list = []
        self.cdf_list = []

        self.threads = threads

        # setup matplotlib parameters
        plt.rcParams['figure.figsize'] = 16, 9

        # make sure a diretory for the report is generated
        if not os.path.exists(os.getcwd() + '/reports'):
            os.mkdir(os.getcwd() + '/reports')
        if not os.path.exists(os.getcwd() + '/reports/' + self.data_name):
            os.mkdir(os.getcwd() + '/reports/' + self.data_name)

        # plot the data before is is sorted and prepared
        self.plot_data()
        self.perpare_data()

    @staticmethod
    def get_significant_digit(number):
        """Retrurns the first non zero digit after decimal point."""

        latter_number_part = str(number).split('.')[1]
        if latter_number_part == '0':
            return 0
        else:
            # leading zeros get removed automatically for integers
            return len(latter_number_part[::-1])

    def perpare_data(self):
        """
        This function prepares the data for the processing. The given data will be sorted in descending order and 
        a small linear increasing value will be added, so that two values won't be the same. Becuase this would 
        cause problems in the calculation of the test statistics.
        """
        # check if any values is peresent more than once
        if not np.unique(self.data).size == self.data.size:
            self.significant_digit_of_data = max(
                [Ftt.get_significant_digit(number) for number in self.data.astype('float')])
            # add random noise below the significant digit to make sure no two values are the same
            while np.unique(self.data).size != self.data.size:
                self.data += np.random.normal(size=self.data.size) / 10 ** (self.significant_digit_of_data + 6)
        self.data[::-1].sort()

    def generate_tails(self, data):
        """
        Generates the tails of a given data set. And and transforms it, so the location of the pareto distribution
        for the returned tail is 0.

        Args:
            data (numpy.ndarray): Data to generate the tail for.

        Yields:
            ndarray: The next tail
        """
        for i in range(1, data.size - 1):
            yield data[:i] - data[i]

    @staticmethod
    def fit_tail(tail):
        """
        Fitting the tail using scipys genpareto and calculating the cdf of the tail for the fitted distribution
        Args:
            tail (numpy.ndarray): tail to fit

        Returns:
            numpy.ndarray, tuple: Cdf of the data for the fitted tail, fit parameters (c, loc, scale).
        """
        # floc is set to zero because the data is expected to be transformed, so the location of the pareto distribution
        #  is 0. Check generate_tails for further information.
        fit_out = genpareto.fit(tail, floc=0)
        # generate distribution with the fitted parameters
        estimated_distribution = genpareto(c=fit_out[0], loc=fit_out[1], scale=fit_out[2])
        # calculate the cdf of the estimated distribution in ascending order
        cdf_of_tail = estimated_distribution.cdf(tail)
        cdf_of_tail.sort()
        return cdf_of_tail, fit_out

    def find_optimal_tail(self):
        """
        The function fits all tails and saves the generated fit information. After all tails have been fitted
        the tail with the minimal AU2 test statistic and the index of the tail are saved.
        
        Returns:
            None
        """
        # make sure all lists are cleaned up
        self.cdf_list = []
        self.rv_list = []
        # fit the tails
        for index, tail in enumerate(self.generate_tails(self.data)):
            print("\t" + str(index) + "/" + str(self.data.size), end='\r', flush=True)
            cdf, fit_out = self.fit_tail(tail)
            self.cdf_list.append(cdf)
            # save rv's
            rv = genpareto(c=fit_out[0],
                           loc=fit_out[1],
                           scale=fit_out[2])
            self.rv_list.append(rv)

        # calculate the test statitics
        self.au_2_data = np.array([au2(tail) for tail in self.cdf_list])
        self.cramer_data = np.array([cramer_von_mises(tail) for tail in self.cdf_list])
        self.anderson_data = np.array([anderson_darling(tail) for tail in self.cdf_list])

        self.optimal_tail_index = self.au_2_data.argmin()
        self.optimal_tail = self.cdf_list[self.au_2_data.argmin()]

    def montecarlo_simulation(self, mc_steps=None):
        """
        Runs Monte Carlo simulation for the optimal position.
        
        Args:
            mc_steps: number of Monte Carlo steps to run.

        Returns:
            float: p-value for the AU2 test statistic
            float: p-value for the Anderson-Darling test statistic
            float: p-value for the Cramér-von Mises test statistic
            int: number of montecarlo steps
            
        Raises:
            RuntimeError is the function gets called, when the fit for the optimal tail start has not been run before.
        """
        if (self.optimal_tail_index is None or
                self.rv_list is None or
                self.cdf_list is None):
            raise RuntimeError("Fits have to run before the Monte Carlo simulation")
        if mc_steps is None:
            mc_steps = self.mc_steps
        # generate mc points
        mc_counter_au2 = 0
        mc_counter_a2 = 0
        mc_counter_w2 = 0

        # make sure every thread has a different seed
        random_state = np.random.RandomState(np.random.seed())

        random_variates = self.rv_list[self.optimal_tail_index].rvs(size=(mc_steps, self.optimal_tail.size), random_state=random_state)
        for index, random_variate in enumerate(random_variates):
            print("\t" + str(index) + "/" + str(mc_steps), end='\r', flush=True)
            fit_out = genpareto.fit(np.sort(random_variate)[::-1], floc=0)
            my_pareto = genpareto(c=fit_out[0], loc=fit_out[1], scale=fit_out[2])
            cdf_of_tail = np.sort(my_pareto.cdf(random_variate))
            if au2(cdf_of_tail) > self.au_2_data[self.optimal_tail_index]:
                mc_counter_au2 += 1
            if anderson_darling(cdf_of_tail) > self.anderson_data[self.optimal_tail_index]:
                mc_counter_a2 += 1
            if cramer_von_mises(cdf_of_tail) > self.cramer_data[self.optimal_tail_index]:
                mc_counter_w2 += 1

        return mc_counter_au2, mc_counter_a2, mc_counter_w2, mc_steps

    def run_montecarlo_simulation(self, mc_steps=None, threads=None):
        """
        Runs the montecarlo simulation and saves the results in class variables.
        
        Args:
            mc_steps: Number of montecarlo steps per thread.
            threads: Number of threads to use.

        Returns:
            None
        """
        if mc_steps is None:
            mc_steps = self.mc_steps
        if threads is None:
            threads = self.threads

        with Pool(threads) as p:
            results = p.map(self.montecarlo_simulation, [mc_steps] * threads)
            for result in results:
                self.save_montecarlo_information(*result)

    def save_montecarlo_information(self, mc_counter_au2, mc_counter_a2, mc_counter_w2, mc_steps):
        self.mc_steps_run += mc_steps
        self.mc_counter_au2 += mc_counter_au2
        self.mc_counter_a2 += mc_counter_a2
        self.mc_counter_w2 += mc_counter_w2
        self.p_value_au2 = self.mc_counter_au2 / self.mc_steps_run
        self.p_value_a2 = self.mc_counter_a2 / self.mc_steps_run
        self.p_value_w2 = self.mc_counter_w2 / self.mc_steps_run
        self.mc_error = 1 / self.mc_steps_run ** 0.5

    def quantil_and_cvar(self, p_values=np.array([0.95, 0.97, 0.99, 0.999])):
        """
        Calculates the quantiles for given p-values.
        Args:
            p_values (np.ndarray): p-values to calculate the quantils and the conditional value at risk for.
                                   Defaults to [0.95, 0.97, 0.99, 0.999]

        Returns:
            None
        
        Raises:
            RuntimeError if the Monte Carlo simulation has not been run befor.
        """
        if self.mc_steps_run == 0:
            raise RuntimeError('Quantil and cvar can only be calculated after the montecarlo simulation')
        else:
            sigma = self.rv_list[self.optimal_tail_index].kwds['scale']
            xi = self.rv_list[self.optimal_tail_index].kwds['c']
            quantile = self.data[self.optimal_tail_index] + sigma / xi * (
                    (self.data.size / (self.optimal_tail_index + 1) * (1 - p_values)) ** -xi - 1)
            self.q = [(p, round(q, self.significant_digit_of_data)) for p, q in zip(p_values, quantile)]
            cond = quantile + (sigma + xi * quantile) / (1 - xi)
            self.cond = [(p, round(c, self.significant_digit_of_data)) for p, c in zip(p_values, cond)]

    def plot_statistics(self):
        """Plots the three test statistics and saves the plot"""
        fig, ax = plt.subplots(1, 1, figsize=(16, 9))
        ax.plot(self.au_2_data, label='AU2')
        ax.plot(self.anderson_data, label='Anderson-Darling')
        ax.plot(self.cramer_data, label='Cramér-von Mises')
        ax.grid()
        ax.set_xlabel(r"Sorted Data Index $k$", fontsize=24)
        ax.set_ylabel("Statistics", fontsize=24)
        ax.tick_params(labelsize=20)
        ax.set_yscale('log')
        ax.legend(loc='best', fontsize=24)
        fig.savefig(os.getcwd() + '/reports/' + self.data_name + '/test_statistics.png')
        plt.close(fig)

    def plot_data(self):
        """Plots the data and saves the plot"""
        fig, ax = plt.subplots(1, 1, figsize=(16, 9))
        ax.plot(np.arange(self.data.size), self.data / self.data.size, label=self.data_name)
        ax.set_xlabel(r"Data Index $i$", fontsize=24)
        ax.set_ylabel(r"$X_i$", fontsize=24)
        ax.tick_params(labelsize=20)
        ax.grid()
        ax.legend(loc='best', fontsize=24)
        fig.savefig(os.getcwd() + '/reports/' + self.data_name + '/data.png')
        plt.close(fig)

    def plot_empirical_distribution(self, closeup=True, save=True):
        """
        Plots the empirical distribution and saves the plot.
        Args:
            closeup (bool): If True the p-value range is set, so that the values of p > 0.95 are shown. This parameter
                            is used to have a closeup in the plot of the empirical distiribution.
            save (bool): If True the plot is saved else, the function just returns None. This is used for the picture in
                        picture of the close up.

        Returns: 
            None
        """
        fig, ax = plt.subplots(1, 1, figsize=(16, 9))
        ax.plot(self.data[::-1], np.arange(1, self.data.size + 1) / self.data.size, '+r', label='Data')
        x = np.arange(self.data[self.optimal_tail_index],
                      self.data[0] * 1.5,
                      (self.data[0] * 1.5 - self.data[self.optimal_tail_index]) / 100)
        tw = (self.optimal_tail_index + 1) / self.data.size
        ax.plot(x, (1 - tw) + self.rv_list[self.optimal_tail_index].cdf(x - x.min()) * tw,
                label="Generalized Pareto distribution")
        ax.legend(loc='upper center', bbox_to_anchor=(0.7, 0.9), fontsize=16)

        if closeup:
            axins = ax.inset_axes([.35, .1, .6, .6])
            axins.plot(x, (1 - tw) + self.rv_list[self.optimal_tail_index].cdf(x - x.min()) * tw,
                       label="Generalized Pareto distribution")
            axins.plot(self.data[::-1], np.arange(1, self.data.size + 1) / self.data.size, '+r', label='Data')
            axins.hlines((0.95, 0.97, 0.99, 0.999), self.data.min(), self.data.max() * 2, alpha=0.3)

            axins.set_yticks([0.95, 0.97, 0.99, 0.999])
            axins.set_xlim([self.data[self.optimal_tail_index] - self.data.std(), self.data[0] * 1.5])
            axins.set_ylim((0.94, 1.005))

        fig.savefig(os.getcwd() + '/reports/' + self.data_name + '/data_empirical.png')
        plt.close(fig)

    def report_html(self):
        """
        Generates a html report with.
        
        Returns:
            None
        """
        # TODO: only push the class dict to the stream function and chnage the html form
        env = Environment(
            loader=FileSystemLoader(os.path.dirname(os.path.realpath(__file__)) + '/templates'),
            autoescape=select_autoescape(['html', 'xml'])
        )
        template = env.get_template('report_base.html')
        template.stream(
            {'data_name': self.data_name,
             'data': self.data,
             'p_value_au2': self.p_value_au2,
             'au2_value': self.au_2_data[self.optimal_tail_index],
             'optimal_point_position': self.optimal_tail_index + 1,
             'montecarlo_steps': self.mc_steps_run,
             'fit_parameter': self.rv_list[self.optimal_tail_index].kwds,
             'quantile': self.q,
             'cond': self.cond,
             'optimal_point': round(self.data[self.optimal_tail_index], self.significant_digit_of_data),
             'w2_value': self.cramer_data[self.optimal_tail_index],
             'a2_value': self.anderson_data[self.optimal_tail_index],
             'data_size': self.data.size,
             'p_value_a2': self.p_value_a2,
             'p_value_w2': self.p_value_w2}
        ).dump(os.getcwd() + '/reports/' + self.data_name + '/' + self.data_name + '.html')

    def run_analysis(self, p_values=np.array([0.95, 0.97, 0.99, 0.999])):
        """
        Runs a complete analysis.
        Args:
            p_values: p values to calculate the quantiles and cvar for.

        Returns:
            None
        """
        print('Runnig fit')
        self.find_optimal_tail()
        print('Running Montecarlo simulation')
        self.run_montecarlo_simulation()
        print('Calculating q and cvar')
        self.quantil_and_cvar(p_values=p_values)
        print('Generating plots')
        self.plot_statistics()
        self.plot_empirical_distribution(save='save')
