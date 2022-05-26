"""Vector Borne Disease Model Class.

Contains class for specific diseases that inherit from the
class VectorBorneDiseaseModel. Also contains methods to read
in root configuration file.

"""

import math
import numpy as np
import yaml
import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from utils import timer
from scipy.integrate import solve_ivp
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt


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
        # Comment out for now since mosquito file and model loop is removed
#         try:
#             self.mosq = np.array(pd.read_csv(self.config_dict['MOSQUITOES_FILE_PATH'])['Sv'])
#         except FileNotFoundError as e:
#             self.logger.exception('Mosquito population input file not found.')
#             raise e
#         else:
#             self.logger.info('Mosquito population data successfully opened'
#                              f' with {len(self.mosq)} days available')

        self.error_check_mosq_initial_states()
    
        #EXTRACT and CALCULATE model run times (moving this here for now instead of model_func() for new fitting method)
        self.t = (0, self.config_dict['DURATION'])
        #need to add the +1 to get the correct step size
        self.t_eval = np.linspace(0, self.config_dict['DURATION'], self.config_dict['DURATION']*self.config_dict['RESOLUTION'] + 1)
        
    @classmethod
    def param_dict(cls, config_file, disease_name, param_dict):
        """Takes in a dictionary file that edits the parameter values for the model"""
        obj = cls(config_file, disease_name)
        obj.error_check_param_dict_names(param_dict)
        dict_keys = list(param_dict.keys())
        param_keys = [i for i in dict_keys if i in list(obj.params.keys())]
        init_keys = [i for i in dict_keys if i in list(obj.initial_states.keys())]
        dispersion = [i for i in dict_keys if i == 'dispersion']
        for k in param_keys:
            obj.params[k] = param_dict[k]
        for j in init_keys:
            obj.initial_states[j] = param_dict[j]
        if len(dispersion) != 0:
            obj.dispersion = param_dict['dispersion']
            
            
        obj.logger.info(f"\n\nParameters and/or initial states changed: {param_dict}\n")
        return obj

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

    def calc_Ih_wnv(self, df, verbose = True):
        """Calcualtees Ih compartment using Poisson distribution for WNV"""
        rng = np.random.default_rng()
        if verbose == True:
            try:
                df['Infected Humans'] = rng.poisson(lam=self.params['eta'] * df['Infected Vectors'])
            except ValueError:
                self.logger.exception(f"Used Normal distribution, lam = {np.trunc(self.params['eta']*df['Infected Vectors'])}")
                #added an absolute value around the sd in the normal so we never have sd < 0
                df['Infected Humans'] = [math.trunc(k) for k in rng.normal(loc = self.params['eta']*df['Infected Vectors'], scale = np.sqrt(abs(self.params['eta']*df['Infected Vectors'])))]
        elif verbose == False:
            try:
                df['Ih'] = rng.poisson(lam=self.params['eta'] * df['Iv'])
            except ValueError:
                self.logger.exception(f"Used Normal distribution, lam = {np.trunc(self.params['eta']*df['Iv'])}")
                #added an absolute value around the sd in the normal so we never have sd < 0
                df['Ih'] = [math.trunc(k) for k in rng.normal(loc = self.params['eta']*df['Iv'], scale = np.sqrt(abs(self.params['eta']*df['Iv'])))]
        return df
    
    def calc_Dh(self, df):
        #note this dependent on the fact that there is a cummulative human (Ch) column\
        #currently not adding a verbose option because this is just for model fitting purposes at the moment
        df['Dh'] = df['Ch'].diff()
        return df
    
    @timer
    def run_model(self, disease_name, verbose = True, calc_daily = False):
        """Runs ODE solver to generate model output"""
        keys = list(self.initial_states.keys())
        self.model_output = np.empty([0, len(keys)])
            
        try:
            sol = solve_ivp(self.model_func, self.t, list(self.initial_states.values()), t_eval=self.t_eval)
            out = sol.y.T
        except Exception as e:
            self.logger.exception('Exception occurred running model')
            raise e
        self.model_output = out
        
        if verbose == True:
            self.df = pd.DataFrame(dict(zip(list(self.state_names_order.values()), self.model_output.T)))
        elif verbose == False:
            self.df = pd.DataFrame(dict(zip(list(self.state_names_order.keys()), self.model_output.T)))
            
        if disease_name == 'wnv':
            self.df = self.calc_Ih_wnv(self.df, verbose = verbose)
        if calc_daily == True:
            self.df = self.calc_Dh(self.df)

    
    def save_output(self, disease_name, sim_labels = False, data = None):
        """Save output to file"""
#         self.df = pd.DataFrame(dict(zip(list(self.state_names_order.values()), self.model_output.T)))
#         if disease_name == 'wnv':
#             self.df = self.calc_Ih_wnv(self.df)
#         self.df['Time'] = self.t_eval
        self.df['Time'] = self.t_eval
    
        if sim_labels == True:
            dict_keys = data.columns
            param_keys = [i for i in dict_keys if i in list(self.params.keys())]
            init_keys = [i for i in dict_keys if i in list(self.initial_states.keys())]
            keys = []
            values = []
            for k in param_keys:
                keys.append(k)
                values.append(str(round(self.params[k],4)))
            for j in init_keys:
                keys.append(j)
                values.append(str(self.initial_states[j]))

            param_all = []
            for i in range(0,len(keys)):
                param_all.append(keys[i])
                param_all.append(values[i])

            self.output_names = '_'.join(param_all)

        if self.config_dict['OUTPUT_TYPE'] == 'csv':
            if sim_labels == True:
                output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                           f'{disease_name}_{self.output_names}_model_output.csv')
                self.df.to_csv(output_path, index=False)
            else:
                output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                           f'{disease_name}_model_output.csv')
                self.df.to_csv(output_path, index=False)
        else:
            if sim_labels == True:
                output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                           f'{disease_name}_{self.output_names}_model_output.parquet')
                pq.write_table(pa.Table.from_pandas(self.df), output_path)
            else:
                output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                       f'{disease_name}_model_output.parquet')
                pq.write_table(pa.Table.from_pandas(self.df), output_path)

        self.logger.info(f'Output saved to {output_path}')     
        
    def plot_output(self, disease_name, sim_labels = False, save_figure = False):
        human_vec = [x for x in self.df.columns if "Human" in x or "Time" in x]
        vector_vec = [x for x in self.df.columns if "Vector" in x or "Time" in x]
    
        if disease_name.lower() == "wnv":
            bird_vec = [x for x in self.df.columns if "Bird" in x or "Time" in x]


        human = self.df[human_vec]
        vector = self.df[vector_vec]

        if disease_name.lower() == "wnv":
            bird = self.df[bird_vec]
            data = {'human': human, 'vector': vector, 'bird': bird}
        elif disease_name.lower() == "dengue":
            data = {'human': human, 'vector': vector}
            
        for i in range(0, len(data.keys())):
            k = data[list(data.keys())[i]]
            n_plot = len(k.columns) -1
            k.plot(x='Time',subplots=True, figsize=(7.5,n_plot*2.5))
            if save_figure == False:
                plt.show()
                plt.close()
            elif save_figure == True:
                if sim_labels == True:
                    plt.savefig(f'plots/{disease_name}_{self.output_names}_{list(data.keys())[i]}.png')
                    plt.close()
                else:
                    plt.savefig(f'plots/{disease_name}_{list(data.keys())[i]}.png') 
                    plt.close()

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

    #commenting out portions refering to mosquito file since model loop is removed
    def error_check_mosq_initial_states(self):
        """check if mosquito initial states are numerical values.

        check if mosquito initial states are positive.

        check if duration is positive.

        check if a long enough mosquito population vector is supplied.

        check if resolution is positive.

        """
#         try:
#             if not all(isinstance(_, (int, float, np.int64)) for _ in self.mosq):
#                 raise TypeError('Mosquito initial states must be numerical values.')
#         except TypeError as e:
#             self.logger.exception('Mosquito initial states must be numerical values.')
#             raise e

#         try:
#             if not all(_ >= 0 for _ in self.mosq):
#                 raise ValueError('Mosquito initial states must be positive')
#         except ValueError as e:
#             self.logger.exception('Mosquito initial states must be positive')
#             raise e

        try:
            if not self.config_dict['DURATION'] > 0:
                raise ValueError('Simulation duration must be positive')
        except ValueError as e:
            self.logger.exception('Simulation duration must be positive')
            raise e

#         try:
#             if self.config_dict['DURATION'] > len(self.mosq):
#                 raise ValueError('Simulation duration exceeds days of'
#                                  ' available mosquito population data.')
#         except ValueError as e:
#             self.logger.exception('Simulation duration exceeds days of'
#                                   ' available mosquito population data.')
#             raise e

        try:
            if not self.config_dict['RESOLUTION'] > 0:
                raise ValueError('Simulation resolution must be positive')
        except ValueError as e:
            self.logger.exception('Simulation resolution must be positive')
            raise e
            
            
    def error_check_param_dict_names(self, param_dict):
        """check parameter dictionary names match model parameter names
        """
        constants = {**self.params, **self.initial_states}
        
        try:
            if len([x for x in list(param_dict.keys()) if ((x in list(constants.keys())) | (x == 'dispersion'))]) != len(list(param_dict.keys())):
                raise ValueError('Parameter dictionary names do not match model parameter names')
        except ValueError as e:
            self.logger.exception('Parameter dictionary names do not match model parameter names')
            raise e
