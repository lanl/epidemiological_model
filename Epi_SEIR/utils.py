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
    """Configures command line argument parser for models_main.py

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
                        type=lambda x: is_valid_file(parser, x))
                        # default='config/local_test_config.yaml')
    parser.add_argument('-d', '--disease_name', action='store',
                        type=lambda x: is_disease(parser, x))
                        # default='dengue')
    parser.add_argument('-f', '--figure', dest='figure', action='store_true')
    parser.set_defaults(figure=False)
    parser.add_argument('-sf', '--save_figure', dest='save_figure', action='store_true')
    parser.set_defaults(save_figure=False)

    return parser

def create_arg_parser_LLM():
    """Configures command line argument parser for models_main.py

    Checks if the argument is a valid file.

    Returns:
        parser object.
    """

    def is_valid_file(parser, arg):
        if not os.path.isfile(arg):
            parser.error(f'File {arg} not found.')
        else:
            return arg
        
    def is_valid_dir(parser, arg):
        if not os.path.isdir(arg):
            parser.error(f'Directory {arg} not found.')
        else:
            return arg

    def is_disease(parser, arg):
        if arg.lower() not in ['wnv', 'dengue']:
            parser.error('Specify [wnv] or [dengue]')
        else:
            return arg

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', action='store',
                        type=lambda x: is_valid_file(parser, x))
                        # default='config/local_test_config.yaml')
    parser.add_argument('-d', '--disease_name', action='store',
                        type=lambda x: is_disease(parser, x))
    parser.add_argument('-m', '--mosquito_dir', action = 'store',
                        type=lambda x: is_valid_dir(parser, x))

    return parser

def create_arg_parser_exp():
    """Configures command line argument parser for models_params.py

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
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-d', '--disease_name', action='store',
                        type=lambda x: is_disease(parser, x))
    parser.add_argument('-l', '--sim_labels', dest='sim_labels', action='store_true')
    parser.set_defaults(sim_labels=False)
    parser.add_argument('-g', '--generate_params', dest='generate_params', action='store_true')
    parser.set_defaults(generate_params=False)
    parser.add_argument('-p', '--param_data_file', action='store', 
                        type=lambda x: is_valid_file(parser, x))
    #giving -f option but only for sequencing through parameters for fitting data
    #ultimately -sf will only do something for models_params and -f will only do something for models_param_guess
    parser.add_argument('-f', '--figure', dest='figure', action='store_true')
    parser.set_defaults(figure=False)
    parser.add_argument('-sf', '--save_figure', dest='save_figure', action='store_true')
    parser.set_defaults(save_figure=False)

    return parser

def create_arg_parser_plot():
    """Configures command line argument parser for models_params.py

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
    parser.add_argument('-of', '--output_file', action='store',
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-d', '--disease_name', action='store',
                        type=lambda x: is_disease(parser, x))
    parser.add_argument('-sf', '--save_figure', dest='save_figure', action='store_true')
    parser.set_defaults(save_figure=False)

    return parser

def create_arg_parser_fit():
    """Configures command line argument parser for models_main.py

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
                        type=lambda x: is_valid_file(parser, x))
                        # default='config/local_test_config.yaml')
    parser.add_argument('-d', '--disease_name', action='store',
                        type=lambda x: is_disease(parser, x))
                        # default='dengue')
    parser.add_argument('-p', '--best_guess_params', dest='best_guess_params', action='store_true')
    parser.set_defaults(best_guess_params=False)
    parser.add_argument('-rm', '--run_fit_model', dest='run_fit_model', action='store_true')
    parser.set_defaults(run_fit_model=False)
    parser.add_argument('-f', '--figure', dest='figure', action='store_true')
    parser.set_defaults(figure=False)
    parser.add_argument('-sf', '--save_figure', dest='save_figure', action='store_true')
    parser.set_defaults(save_figure=False)

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

    logger.handlers = []

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
