"""Vector Borne Disease Model Class.

Contains class for specific diseases that inherit from the
class VectorBorneDiseaseModel. Also contains methods to read
in root configuration file and instantiate loggers.

    Typical usage example:

    logger = VBDM.create_logger(logger_name, config_file_name)
    disease = VBDM.DengueSEIRModel(config_file_name, simulation_duration)
"""

import numpy as np
# from abc import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import yaml
import logging
import os
import sys
from datetime import datetime
import argparse
import pandas as pd


def create_arg_parser():
    """Configures command line argument parser.

    Checks if the argument is a valid file.

    Returns:
        parser object.
    """
    def is_valid_file(parser, arg):
        if not os.path.isfile(arg):
            parser.error(f'File {arg} not found.')
        else:
            return arg

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--config_file', action='store',
                        type=lambda x: is_valid_file(parser, x),
                        default='config/config.yaml')

    return parser


def create_logger(name, config_file):
    """Configures and instantiates logger object.

    Creates logger that prints DEBUG statements and more
    serious statements to log files. Prints ERROR statements
    to console.

    Args:
        name: Name of logger to be created. Typically module name. String.
        config_file: Absolute path to root config file. String.

    Returns:
        logging object.
    """

    with open(config_file, 'r') as in_file:
        logfile_path = yaml.safe_load(in_file)['LOGFILE_PATH']

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('[%(asctime)s] %(name)s \
                %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    file_handler = logging.FileHandler(os.path.join(logfile_path,
                                                    datetime.now().strftime(f'{name}_%Y-%m-%d.log')))
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.ERROR)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


class VectorBorneDiseaseModel():
    """Defines a general vector borne disease model.

    Reads root config file and sets parameters.

    Attributes:
        params: ODE system parameters.
        initial_states: initial states for population sizes.

    """

    def __init__(self, config_file, config_name, days):
        # self.success = True
        try:
            self._read_config(config_file, config_name)
        except Exception:
            logger.exception(f'Exception occured opening config file for \
                             {config_name}')
            # self.success = False
            sys.exit(1)
        else:
            logger.info('Config file read in')

        self.t = np.linspace(0, days, days*500)

    def _read_config(self, config_file, config_name):
        """Reads root configuration file"""
        with open(config_file, 'r') as in_file:
            # self.config_dict = yaml.safe_load(in_file)[config_name]
            self.config_dict = yaml.safe_load(in_file)
            self.mosq = np.array(pd.read_csv(self.config_dict['MOSQUITOES_FILE_PATH'])['mosq'])
            self.params = self.config_dict[config_name]['PARAMETERS']
            self.initial_states = self.config_dict[config_name]['INITIAL_STATES']

    def save_model(self):
        pd.DataFrame(self.model_output).to_csv(os.path.join(self.config_dict['OUTPUT_DIR'],
                                                            'human_cases.csv'))
        # print(self.model_output)


class DengueSEIRModel(VectorBorneDiseaseModel):

    """Models the spread of dengue.

    Inherits from the VectorBorneDiseaseModel class. Solves ODE system
    of equations and plots the resulting curves.

    """

    def __init__(self, config_file, days):
        super().__init__(config_file, 'DENGUE', days)
        self.initial_states['Sv'] = self.mosq[0]
        print("DENGUE INITIAL", self.initial_states)

    def run_model(self):
        """Runs ODE solver to generate model output"""
        y0 = self.initial_states['Sh'], self.initial_states['Eh'], \
            self.initial_states['Iha'], self.initial_states['Ihs'], \
            self.initial_states['Rh'], self.initial_states['Sv'], \
            self.initial_states['Ev'], self.initial_states['Iv']

        try:
            self.model_output = odeint(self._model_dengue, y0,
                                       self.t, args=(self,))
        except Exception:
            logger.exception('Exception occured running dengue model')
            # self.success = False
            sys.exit(1)
        else:
            logger.info('dengue model run complete')

    def _model_dengue(self, y, t, p):
        """Defines system of ODEs for dengue model"""
        # States and population
        Sh, Eh, Iha, Ihs, Rh, Sv, Ev, Iv = y
        N_h = sum([Sh, Eh, Iha, Ihs, Rh])
        N_v = sum([Sv, Ev, Iv])

        # Biting rate
        b = self.params['sigma_h'] * self.params['sigma_v'] / \
            (self.params['sigma_h'] * N_h + self.params['sigma_v'] * N_v)
        b_h = b * N_v
        b_v = b * N_h

        # Force of infecton
        lambda_h = b_h * self.params['beta_h'] * Iv / N_v
        lambda_v = b_v * self.params['beta_v'] * (Iha + Ihs) / N_h

        # System of equations
        dSh = -lambda_h * Sh
        dEh = lambda_h * Sh - self.params['nu_h'] * Eh
        dIha = self.params['psi'] * self.params['nu_h'] * Eh - \
            self.params['gamma_h'] * Iha
        dIhs = (1 - self.params['psi']) * self.params['nu_h'] * \
            Eh - self.params['gamma_h'] * Ihs
        dRh = self.params['gamma_h'] * (Iha + Ihs)
        dSv = -lambda_v * Sh
        dEv = lambda_v * Sh - self.params['nu_v'] * Ev
        dIv = self.params['nu_v'] * Ev - self.params['mu_v'] * Iv

        return dSh, dEh, dIha, dIhs, dRh, dSv, dEv, dIv

    def graph_model(self):
        """Plots output of dengue model"""
        Sh, Eh, Iha, Ihs, Rh, Sv, Ev, Iv = self.model_output.T

        fig = plt.figure(facecolor='w', figsize=[2*6.4, 2*4.8])
        ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Cases')
        ax.yaxis.set_tick_params(length=0)
        ax.xaxis.set_tick_params(length=0)
        ax.grid(b=True, which='major', c='w', lw=2, ls='-')

        for spine in ('top', 'right', 'bottom', 'left'):
            ax.spines[spine].set_visible(False)

        ax.plot(self.t, Sh, 'b', alpha=0.5, lw=2, label='Susceptible Humans')
        ax.plot(self.t, Iha+Ihs, 'r', alpha=0.5, lw=2, label='Infected Humans')
        ax.plot(self.t, Iv, 'k', alpha=0.5, lw=2, label='Infected Vectors')
        legend = ax.legend()
        legend.get_frame().set_alpha(0.5)
        plt.title("Dengue Incidence")

        plt.show()


class WNVSEIRModel(VectorBorneDiseaseModel):

    """Models the spread of WNV.

    Inherits from the VectorBorneDiseaseModel class. Solves ODE system
    of equations and plots the resulting curves.

    """

    def __init__(self, config_file, days):
        super().__init__(config_file, 'WNV', days)
        self.initial_states['Sv'] = self.mosq[0]
        print("WNV initial", self.initial_states)

    def run_model(self):
        """Runs ODE solver to generate model output"""
        y0 = self.initial_states['Sh'], self.initial_states['Eh'], \
            self.initial_states['Iha'], self.initial_states['Ihs'], \
            self.initial_states['Rh'], self.initial_states['Sv'], \
            self.initial_states['Ev'], self.initial_states['Iv']

        try:
            self.model_output = odeint(self._model_WNV, y0,
                                       self.t, args=(self,))
        except Exception:
            logger.exception('Exception occured running WNV model')
            sys.exit(1)
            # self.success = False
        else:
            logger.info('WNV model run complete')

    def _model_WNV(self, y, t, p):
        """Defines system of ODEs for WNV model"""
        # States and population
        Sh, Eh, Iha, Ihs, Rh, Sv, Ev, Iv = y
        N_h = sum([Sh, Eh, Iha, Ihs, Rh])
        N_v = sum([Sv, Ev, Iv])

        # Biting rate
        b = self.params['sigma_h'] * self.params['sigma_v'] / \
            (self.params['sigma_h'] * N_h + self.params['sigma_v'] * N_v)
        b_h = b * N_v
        b_v = b * N_h

        # Force of infecton
        lambda_h = b_h * self.params['beta_h'] * Iv / N_v
        lambda_v = b_v * self.params['beta_v'] * (Iha + Ihs) / N_h

        # System of equations
        dSh = -lambda_h * Sh
        dEh = lambda_h * Sh - self.params['nu_h'] * Eh
        dIha = self.params['psi'] * self.params['nu_h'] * Eh - \
            self.params['gamma_h'] * Iha
        dIhs = (1 - self.params['psi']) * self.params['nu_h'] * \
            Eh - self.params['gamma_h'] * Ihs
        dRh = self.params['gamma_h'] * (Iha + Ihs)
        dSv = -lambda_v * Sh
        dEv = lambda_v * Sh - self.params['nu_v'] * Ev
        dIv = self.params['nu_v'] * Ev - self.params['mu_v'] * Iv

        return dSh, dEh, dIha, dIhs, dRh, dSv, dEv, dIv

    def graph_model(self):
        """Plots output of WNV model"""
        Sh, Eh, Iha, Ihs, Rh, Sv, Ev, Iv = self.model_output.T

        fig = plt.figure(facecolor='w', figsize=[2*6.4, 2*4.8])
        ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Cases')
        ax.yaxis.set_tick_params(length=0)
        ax.xaxis.set_tick_params(length=0)
        ax.grid(b=True, which='major', c='w', lw=2, ls='-')

        for spine in ('top', 'right', 'bottom', 'left'):
            ax.spines[spine].set_visible(False)

        ax.plot(self.t, Sh, 'b', alpha=0.5, lw=2, label='Susceptible Humans')
        ax.plot(self.t, Iha+Ihs, 'r', alpha=0.5, lw=2, label='Infected Humans')
        ax.plot(self.t, Iv, 'k', alpha=0.5, lw=2, label='Infected Vectors')
        legend = ax.legend()
        legend.get_frame().set_alpha(0.5)
        plt.title("WNV Incidence")

        plt.show()


if 'sphinx-build' in sys.argv[0]:
    logger = create_logger(__name__, '_default_config.yaml')
    # TODO phase out usage of sys.argv. Problem now is if program is run with
    # sys.argv = ['/Users/jkeithley/opt/anaconda3/bin/sphinx-build', '-M', 'html', '.', '_build']
else:
    parser = create_arg_parser()
    args = parser.parse_args()
    logger = create_logger(__name__, args.config_file)
