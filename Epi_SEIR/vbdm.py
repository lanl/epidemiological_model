"""Vector Borne Disease Model Class.

Contains class for specific diseases that inherit from the
class VectorBorneDiseaseModel. Also contains methods to read
in root configuration file.

"""

import numpy as np
import yaml
import os
import sys
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from utils import create_arg_parser, timer
from scipy.integrate import odeint
from abc import ABC, abstractmethod


class VectorBorneDiseaseModel(ABC):
    """Defines a general vector borne disease model.

    Reads root config file and sets parameters.

    Attributes:
        config_dict: All configurations read in from config file.
        params: ODE system parameters.
        initial_states: initial states for population sizes.
        mosq: Daily mosquito population data.

    """

    def __init__(self, config_file, disease_name):
        self._read_config(config_file, disease_name)

        # Read parameters
        self.params = self.config_dict[disease_name]['PARAMETERS']
        try:
            if not all(_ >= 0 for _ in self.params.values()):
                raise ValueError('Model parameters must be positive')
        except ValueError as e:
            self.logger.exception('Model parameters must be positive')
            raise e

        # Read initial states
        try:
            with open(self.config_dict['INITIAL_STATES_FILE_PATH'], 'r') as in_file:
                self.initial_states = yaml.safe_load(in_file)[disease_name]['INITIAL_STATES']
        except FileNotFoundError as e:
            self.logger.exception('Initial states input file not found.')
            raise e
        self.logger.info('Initial states data successfully opened')

        try:
            if not all(_ >= 0 for _ in self.initial_states.values()):
                raise ValueError('Model initial states must be positive')
        except ValueError as e:
            self.logger.exception('Model initial states must be positive')
            raise e

        # Read mosquito initial states
        try:
            self.mosq = np.array(pd.read_csv(self.config_dict['MOSQUITOES_FILE_PATH'])['Sv'])
        except FileNotFoundError as e:
            self.logger.exception('Mosquito population input file not found.')
            raise e
        else:
            self.logger.info('Mosquito population data successfully opened'
                            f' with {len(self.mosq)} days available')

        try:
            if not all(i >= 0 for i in self.mosq):
                raise ValueError('Mosquito initial states must be positive')
        except ValueError as e:
            self.logger.exception('Mosquito initial states must be positive')
            raise e

        # Check duration
        try:
            if not self.config_dict['DURATION'] > 0:
                raise ValueError('Simulation duration must be positive')
        except ValueError as e:
            self.logger.exception('Simulation duration must be positive')
            raise e

        try:
            if self.config_dict['DURATION'] > len(self.mosq):
                raise ValueError('Simulation duration exceeds days of'
                                 ' available mosquito population data.')
        except ValueError as e:
            self.logger.exception('Simulation duration exceeds days of'
                                  ' available mosquito population data.')
            raise e

        # Check resolution
        try:
            if not self.config_dict['RESOLUTION'] > 0:
                raise ValueError('Simulation resolution must be positive')
        except ValueError as e:
            self.logger.exception('Simulation resolution must be positive')
            raise e

    def _read_config(self, config_file, disease_name):
        """Reads root configuration file"""
        try:
            with open(config_file, 'r') as in_file:
                self.config_dict = yaml.safe_load(in_file)
        except OSError as e:
            self.logger.exception('Exception occured opening configuration file')
            raise e
        else:
            self.logger.info('Configuration file successfully opened')

    def __str__(self):
        return f"Running model for {self.config_dict['DURATION']} days at" \
               f" {self.config_dict['RESOLUTION']} points per day"

    @abstractmethod
    def set_y0(self):
        pass

    @abstractmethod
    def model_func(self, y, t):
        pass

    @timer
    def run_model(self):
        """Runs ODE solver to generate model output"""
        keys = ['Sh', 'Eh', 'Iha', 'Ihs', 'Rh', 'Sv', 'Ev', 'Iv']
        self.model_output = np.empty([0, len(keys)])

        t = np.linspace(0, 1, self.config_dict['RESOLUTION'] + 1)

        for i in range(0, self.config_dict['DURATION']):
            self.initial_states['Sv'] = self.mosq[i]
            y0 = self.set_y0()

            try:
                out = odeint(self.model_func, y0, t)
            except Exception as e:
                self.logger.exception('Exception occured running dengue model')
                raise e

            self.initial_states = dict(zip(keys, out[-1]))

            out = out[:-1]

            self.model_output = np.concatenate((self.model_output, out))

    def save_output(self, disease_name):
        """Save output to file"""
        keys = ['Susceptible Humans', 'Exposed Humans',
                'Asymptomatic Infected Humans', 'Symptomatic Infected Humans',
                'Recovered Humans', 'Susceptible Vectors',
                'Exposed Vectors', 'Infected Vectors']
        df = pd.DataFrame(dict(zip(keys, self.model_output.T)))

        if self.config_dict['OUTPUT_TYPE'].lower() == 'csv':
            output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                       f'{disease_name}_model_output.csv')
            df.to_csv(output_path)
        else:
            if self.config_dict['OUTPUT_TYPE'].lower() != 'parquet':
                self.logger.error(f'Output file type'
                                  ' {self.config_dict["OUTPUT_TYPE"].lower()}'
                                  ' not recognized. Output will be .parquet file')

            output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                       f'{disease_name}_model_output.parquet')
            pq.write_table(pa.Table.from_pandas(df), output_path)

        self.logger.info(f'Output saved to {output_path}')


if not sys.argv[0].endswith('sphinx-build'):
    parser = create_arg_parser()
    args = parser.parse_args()
