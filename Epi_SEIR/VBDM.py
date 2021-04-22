"""Vector Borne Disease Model Class.

Contains class for specific diseases that inherit from the
class VectorBorneDiseaseModel. Also contains methods to read
in root configuration file and instantiate loggers.

"""

import numpy as np
import yaml
import os
import sys
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from utils import create_arg_parser
from scipy.integrate import odeint
import abc


class VectorBorneDiseaseModel():
    """Defines a general vector borne disease model.

    Reads root config file and sets parameters.

    Attributes:
        params: ODE system parameters.
        initial_states: initial states for population sizes.

    """

    def __init__(self, config_file, disease_name):
        self._read_config(config_file, disease_name)

        # Read parameters
        self.params = self.config_dict[disease_name]['PARAMETERS']
        try:
            if not all(_ >= 0 for _ in self.params.values()):
                raise ValueError("Model parameters must be positive")
        except ValueError:
            self.logger.exception("Model parameters must be positive")

        # Read initial states
        self.initial_states = self.config_dict[disease_name]['INITIAL_STATES']
        try:
            if not all(_ >= 0 for _ in self.initial_states.values()):
                raise ValueError("Model initial states must be positive")
        except ValueError:
            self.logger.exception("Model initial states must be positive")

        # Read mosquito initial states
        self.mosq = np.array(pd.read_csv(self.config_dict['MOSQUITOES_FILE_PATH'])['mosq'])
        try:
            if not all(i >= 0 for i in self.mosq):
                raise ValueError("Mosquito initial states must be positive")
        except ValueError:
            self.logger.exception("Mosquito initial states must be positive")

        # TODO error check duration and resolution

    def _read_config(self, config_file, disease_name):
        """Reads root configuration file"""
        try:
            with open(config_file, 'r') as in_file:
                self.config_dict = yaml.safe_load(in_file)
        except OSError:
            self.logger.exception('Exception occured opening configuration file')
        else:
            self.logger.info("Configuration file successfully opened")

    @abc.abstractmethod
    def set_y0(self):
        pass
        # TODO raise exception if not set?

    @abc.abstractmethod
    def model_func(self, y, t, p):
        pass

    def run_model(self):
        # TODO move to VBDM module
        """Runs ODE solver to generate model output"""
        # TODO run_model STILL UNDER CONSTRUCTION FOR TIME DEPENDENT MOSQUITO
        y0 = self.set_y0()

        t = np.linspace(0, 1, self.config_dict['RESOLUTION'] + 1)

        try:
            self.model_output = odeint(self.model_func, y0,
                                       t, args=(self,))
        # TODO work on reraising exceptions
        except Exception as e:
            self.logger.exception('Exception occured running dengue model')
            raise e

        y0 = tuple(self.model_output[-1])
        self.model_output = self.model_output[:-1]

        for i in range(1, self.config_dict['DURATION']):
            self.initial_states['Sv'] = self.mosq[i]

            try:
                out = odeint(self.model_func, y0,
                             t, args=(self,))
            except Exception:
                self.logger.exception('Exception occured running dengue model')
                sys.exit(1)

            y0 = tuple(out[-1])
            out = out[:-1]

            self.model_output = np.concatenate((self.model_output, out))

    def save_model(self, disease_name):
        """Save output to file"""
        keys = ['Sh', 'Eh', 'Iha', 'Ihs', 'Rh', 'Sv', 'Ev', 'Iv']
        df = pd.DataFrame(dict(zip(keys, self.model_output.T)))

        if self.config_dict['OUTPUT_TYPE'].lower() == 'csv':
            output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                       f'{disease_name}_model_output.csv')
            df.to_csv(output_path)
        else:
            if self.config_dict['OUTPUT_TYPE'].lower() != 'parquet':
                self.logger.error(f'Output file type {self.config_dict["OUTPUT_TYPE"].lower()} not recognized. Output will be .parquet file')

            output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
                                       f'{disease_name}_model_output.parquet')
            pq.write_table(pa.Table.from_pandas(df), output_path)

        # print(self.model_output)


if not sys.argv[0].endswith('sphinx-build'):
    # print("FLAG----------------- entered docsys VBDM module")
    # logger = create_logger(__name__, '_default_config.yaml')
    # TODO phase out usage of sys.argv. Problem now is if program is run with
    # sys.argv = ['/Users/jkeithley/opt/anaconda3/bin/sphinx-build', '-M', 'html', '.', '_build']
    # else:
    parser = create_arg_parser()
    args = parser.parse_args()
    # logger = create_logger(__name__, args.config_file)
