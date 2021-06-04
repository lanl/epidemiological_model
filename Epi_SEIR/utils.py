"""Utilities script to create argument parser and loggers

    Typical usage example:

    logger = VBDM.create_logger(logger_name, config_file_name)
    disease = VBDM.DengueSEIRModel(config_file_name, simulation_duration)

"""

import numpy as np
import logging
import os
import argparse
import yaml
from datetime import datetime
from functools import wraps
import time


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

    def is_disease(parser, arg):
        if arg.lower() not in ['wnv', 'dengue']:
            parser.error('Specify [wnv] or [dengue]')
        else:
            return arg

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', action='store',
                        type=lambda x: is_valid_file(parser, x),
                        default='config/config.yaml')
    parser.add_argument('-d', '--disease_name', action='store',
                        type=lambda x: is_disease(parser, x))

    return parser


def create_logger(name, config_file):
    """Configures and instantiates logger object.

    Creates logger that prints DEBUG statements and more
    serious statements to log files. Prints ERROR statements
    to console.

    Args:
        name: Name of logger to be created. Typically module name. String.
        config_file: Absolute path to root config file. String.

    Return:
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

def timer(func):
    """Print the runtime of the decorated function"""
    @wraps(func)
    def wrapper_timer(self, *args, **kwargs):
        start_time = time.perf_counter()
        value = func(self, *args, **kwargs)
        run_time = time.perf_counter() - start_time
        self.logger.info(f'Finished {func.__name__!r} in {run_time:.4f} seconds')
        return value
    return wrapper_timer

def error_check_state_names(self):
    try:
        if not all(isinstance(_, str) for _ in self.state_names_order):
            raise TypeError('Initial state names must be strings')
    except TypeError as e:
        self.logger.exception('Initial state names must be strings')
        raise e

def error_check_positions(self):
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
    try:
        if not all(_ >= 0 for _ in self.initial_states.values()):
            raise ValueError('Model initial states must be positive')
    except ValueError as e:
        self.logger.exception('Model initial states must be positive')
        raise e
    except TypeError as e:
        self.logger.exception('Initial states must be numerical values.'
                              ' Initialize all initial states.')
        raise e

def error_check_mosq_initial_states(self):
    try:
        if not all(i >= 0 for i in self.mosq):
            raise ValueError('Mosquito initial states must be positive')
    except ValueError as e:
        self.logger.exception('Mosquito initial states must be positive')
        raise e
    except TypeError as e:
        self.logger.exception('Mosquito initial states must be numerical values.'
                              ' Initialize all mosquito initial states.')
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
