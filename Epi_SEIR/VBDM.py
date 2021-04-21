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
from utils import create_arg_parser, create_logger


class VectorBorneDiseaseModel():
    """Defines a general vector borne disease model.

    Reads root config file and sets parameters.

    Attributes:
        params: ODE system parameters.
        initial_states: initial states for population sizes.

    """

    def __init__(self, config_file, disease_name, days):
        # self.success = True
        try:  # TODO get rid  of this try except block
            self._read_config(config_file, disease_name)
        except Exception:
            logger.exception(f'Exception occured opening config file for \
                             {disease_name}')
            # self.success = False
            sys.exit(1)
        else:
            logger.info('Config file read in')

        self.t = np.linspace(0, days, days*500)

    def _read_config(self, config_file, disease_name):
        """Reads root configuration file"""
        with open(config_file, 'r') as in_file:
            self.config_dict = yaml.safe_load(in_file)

        # Read parameters
        self.params = self.config_dict[disease_name]['PARAMETERS']
        try:
            if not all(_ >= 0 for _ in self.params.values()):
                raise ValueError("Model parameters must be positive")
        except ValueError:
            logger.exception("Model parameters must be positive")

        # Read initial states
        self.initial_states = self.config_dict[disease_name]['INITIAL_STATES']
        try:
            if not all(_ >= 0 for _ in self.initial_states.values()):
                raise ValueError("Model initial states must be positive")
        except ValueError:
            logger.exception("Model initial states must be positive")

        # Read mosquito initial states
        self.mosq = np.array(pd.read_csv(self.config_dict['MOSQUITOES_FILE_PATH'])['mosq'])
        try:
            if not all(i >= 0 for i in self.mosq):
                raise ValueError("Mosquito initial states must be positive")
        except ValueError:
            logger.exception("Mosquito initial states must be positive")

    def save_model(self):
        """Save output to file"""
        Sh, Eh, Iha, Ihs, Rh, Sv, Ev, Iv = self.model_output.T
        pd.DataFrame({'Sh': Sh, 'Eh': Eh, 'Iha': Iha, 'Ihs': Ihs, 'Rh': Rh,
                      'Sv': Sv, 'Ev': Ev, 'Iv': Iv}).to_csv(os.path.join(self.config_dict['OUTPUT_DIR'],
                                                            'human_cases.csv'))

        # print(self.model_output)


if sys.argv[0].endswith('sphinx-build'):
    logger = create_logger(__name__, '_default_config.yaml')
    # TODO phase out usage of sys.argv. Problem now is if program is run with
    # sys.argv = ['/Users/jkeithley/opt/anaconda3/bin/sphinx-build', '-M', 'html', '.', '_build']
else:
    parser = create_arg_parser()
    args = parser.parse_args()
    logger = create_logger(__name__, args.config_file)
