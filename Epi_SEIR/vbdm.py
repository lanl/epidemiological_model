"""Vector Borne Disease Model Class.

Contains class for specific diseases that inherit from the
class VectorBorneDiseaseModel. Also contains methods to read
in root configuration file.

"""

import numpy as np
import yaml
import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from utils import timer
from scipy.integrate import solve_ivp
from abc import ABC, abstractmethod


class VectorBorneDiseaseModel(ABC):
    """Defines a general vector borne disease model.

    Reads root config file and sets parameters.

    Attributes:
        config_dict: All configurations read in from config file.\n
        params: ODE system parameters.\n

        initial_states: initial states for population sizes.\n
        mosq: Daily mosquito population data.\n

    """

    def __init__(self, config_file, disease_name):
        self._read_config(config_file, disease_name)

        self.error_check_output_type()

        # Read parameters
        self.params = self.config_dict[disease_name]['PARAMETERS']
        
        self.param_of_interest = list(self.config_dict[disease_name]['PARAMETER_OF_INTEREST'].keys())[0]
        
        self.start_value = self.config_dict[disease_name]['PARAMETER_OF_INTEREST'][self.param_of_interest]['start']
        
        self.stop_value = self.config_dict[disease_name]['PARAMETER_OF_INTEREST'][self.param_of_interest]['stop']
        
        self.step_value = self.config_dict[disease_name]['PARAMETER_OF_INTEREST'][self.param_of_interest]['step']

        self.logger.info(f"\n\nParameters for model: {self.params}\n")

        # Read initial states
        try:
            with open(self.config_dict['INITIAL_STATES_FILE_PATH'], 'r') as in_file:
                self.initial_states = yaml.safe_load(in_file)[disease_name]
        except FileNotFoundError as e:
            self.logger.exception('Initial states input file not found.')
            raise e
        else:
            self.logger.info('Initial states data successfully opened')

        # SORT entire initial states dictionary
        self.initial_states = dict(sorted(self.initial_states.items(), key=lambda x: x[1]['position']))

        # EXTRACT initial states order names
        self.state_names_order = dict(zip(self.initial_states.keys(),
                                          [list(self.initial_states.values())[i]['name']
                                          for i in range(len(self.initial_states.values()))]))

        self.error_check_state_names()

        self.error_check_positions()

        # EXTRACT values, set as initial states
        self.initial_states = dict(zip(self.initial_states.keys(),
                                       [list(self.initial_states.values())[i]['value']
                                       for i in range(len(self.initial_states.values()))]))

        self.logger.info(f"\n\nInitial states for model: {self.initial_states}\n")

        self.error_check_initial_states()

        # Read mosquito initial states
        try:
            self.mosq = np.array(pd.read_csv(self.config_dict['MOSQUITOES_FILE_PATH'])['Sv'])
        except FileNotFoundError as e:
            self.logger.exception('Mosquito population input file not found.')
            raise e
        else:
            self.logger.info('Mosquito population data successfully opened'
                             f' with {len(self.mosq)} days available')

        self.error_check_mosq_initial_states()

    def _read_config(self, config_file, disease_name):
        """Reads root configuration file"""
        try:
            with open(config_file, 'r') as in_file:
                self.config_dict = yaml.safe_load(in_file)
        except OSError as e:
            self.logger.exception('Exception occurred opening configuration file')
            raise e
        else:
            self.logger.info('Configuration file successfully opened')

    def __str__(self):
        return f"Running model for {self.config_dict['DURATION']} days at" \
               f" {self.config_dict['RESOLUTION']} points per day"

    @abstractmethod
    def model_func(self, t, y):
        pass

    @timer
    def run_model(self, disease_name):
        """Runs ODE solver to generate model output"""
        keys = list(self.initial_states.keys())
        self.model_output = np.empty([0, len(keys)])
            
        t = (0, self.config_dict['DURATION'])
        t_eval = np.linspace(0, self.config_dict['DURATION'], self.config_dict['DURATION']*self.config_dict['RESOLUTION'])
            
        try:
            sol = solve_ivp(self.model_func, t, list(self.initial_states.values()), t_eval=t_eval)
            out = sol.y.T
        except Exception as e:
            self.logger.exception('Exception occurred running model')
            raise e
        self.model_output = out

    def save_output(self, disease_name, sim_labels = 'F'):
        """Save output to file"""
        df = pd.DataFrame(dict(zip(list(self.state_names_order.values()), self.model_output.T)))

        if self.config_dict['OUTPUT_TYPE'] == 'csv':
            if sim_labels == 'T':
                output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                           f'{disease_name}_{self.param_of_interest}_{self.params[self.param_of_interest]}_model_output.csv')
                df.to_csv(output_path, index=False)
            else:
                output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                           f'{disease_name}_model_output.csv')
                df.to_csv(output_path, index=False)
        else:
            if sim_labels == 'T':
                output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                           f'{disease_name}_{self.param_of_interest}_{self.params[self.param_of_interest]}_model_output.parquet')
                pq.write_table(pa.Table.from_pandas(df), output_path)
            else:
                output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                       f'{disease_name}_model_output.parquet')
                pq.write_table(pa.Table.from_pandas(df), output_path)

        self.logger.info(f'Output saved to {output_path}')


    def error_check_output_type(self):
        """check if output type is a string.

        check if output type is .csv or .parquet

        """

        try:
            if not isinstance(self.config_dict['OUTPUT_TYPE'], str):
                raise TypeError('Output type must be a string')
        except TypeError as e:
            self.logger.exception('Output type must be a string')
            raise e

        self.config_dict['OUTPUT_TYPE'] = self.config_dict['OUTPUT_TYPE'].strip('.').lower()

        try:
            if (self.config_dict['OUTPUT_TYPE'] != 'csv') and (self.config_dict['OUTPUT_TYPE'] != 'parquet'):
                raise ValueError('Output type must be .csv or .parquet')
        except ValueError as e:
            self.logger.exception('Output type must be .csv or .parquet')
            raise e

    def error_check_state_names(self):
        """check if compartment names field is string type"""

        try:
            if not all(isinstance(_, str) for _ in self.state_names_order.values()):
                raise TypeError('Initial state names must be strings')
        except TypeError as e:
            self.logger.exception('Initial state names must be strings')
            raise e

    def error_check_positions(self):
        """check if position array is unique.

        check if position array is int type.

        check if positions array is positive.

        """

        # EXTRACT list of position
        positions = [list(self.initial_states.values())[i]['position'] for i in
                     range(len(self.initial_states.values()))]

        try:
            if len(positions) != len(np.unique(positions)):
                raise ValueError('Position values must be unique')
        except ValueError as e:
            self.logger.exception('Position values must be unique')
            raise e

        try:
            if not all(isinstance(_, int) for _ in positions):
                raise TypeError('Position values must be integers')
        except TypeError as e:
            self.logger.exception('Position values must be integers')
            raise e

        try:
            if not all(_ >= 0 for _ in positions):
                raise ValueError('Position values must be positive')
        except ValueError as e:
            self.logger.exception('Position values must be positive')
            raise e

    def error_check_initial_states(self):
        """check if initial states are numerical values.

        check if initial states are positive.

        """
        try:
            if not all(isinstance(_, (int, float, np.int64)) for _ in self.initial_states.values()):
                raise TypeError('Initial states must be numerical values.'
                                ' Initialize all initial states.')
        except TypeError as e:
            self.logger.exception('Initial states must be numerical values.'
                                  ' Initialize all initial states.')
            raise e

        try:
            if not all(_ >= 0 for _ in self.initial_states.values()):
                raise ValueError('Model initial states must be positive')
        except ValueError as e:
            self.logger.exception('Model initial states must be positive')
            raise e

    def error_check_mosq_initial_states(self):
        """check if mosquito initial states are numerical values.

        check if mosquito initial states are positive.

        check if duration is positive.

        check if a long enough mosquito population vector is supplied.

        check if resolution is positive.

        """
        try:
            if not all(isinstance(_, (int, float, np.int64)) for _ in self.mosq):
                raise TypeError('Mosquito initial states must be numerical values.')
        except TypeError as e:
            self.logger.exception('Mosquito initial states must be numerical values.')
            raise e

        try:
            if not all(_ >= 0 for _ in self.mosq):
                raise ValueError('Mosquito initial states must be positive')
        except ValueError as e:
            self.logger.exception('Mosquito initial states must be positive')
            raise e

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

        try:
            if not self.config_dict['RESOLUTION'] > 0:
                raise ValueError('Simulation resolution must be positive')
        except ValueError as e:
            self.logger.exception('Simulation resolution must be positive')
            raise e
