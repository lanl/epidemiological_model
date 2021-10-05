"""Script to fit vector borne disease models to data
Instantiates VBDM (Vector Borne Disease Model)
and fitting process given paramters from local_param_config.yaml and config.yaml configuration files.

    Typical usage example:

    python models_fit.py -c <absolute path to config.yaml> -d <disease name>
"""
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser
from generate_param_values import param_data
import sys
import numpy as np
import pandas as pd
import yaml


def fit():
    """Runs vector borne disease simulations.

    Instaniates each disease from configuration file
    """
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser()
        args, unknown = parser.parse_known_args()

    disease_name = args.disease_name.lower()

    if disease_name == 'dengue':
        disease = DengueSEIRModel(args.config_file)
    elif disease_name == 'wnv':
        disease = WNVSEIRModel(args.config_file)

    disease.logger.info(disease)
    disease.fit_constants()
    disease.save_fit_output(disease_name)
    disease.logger.info('SUCCESS')
    
    if args.run_fit_model == True:
        fit_params = pd.read_csv(f'param_fit_output/{disease_name}_parameter_values.csv)
        fit_disease = DengueSEIRModel.param_dict(config_file = args.config_file, param_dict =  fit_params)
        disease.logger.info(disease)
        disease.run_model(disease_name)
        disease.save_output(disease_name, args.sim_labels, fit_params)
        disease.logger.info('SUCCESS')
    
    if args.plot == True:
        plot = 1

if __name__ == "__main__":
    fit()
