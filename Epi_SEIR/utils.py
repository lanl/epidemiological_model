"""Utilities script to create argument parser and loggers

    Typical usage example:

    logger = VBDM.create_logger(logger_name, config_file_name)
    disease = VBDM.DengueSEIRModel(config_file_name, simulation_duration)

"""

import logging
import os
import argparse
import yaml
from datetime import datetime


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
