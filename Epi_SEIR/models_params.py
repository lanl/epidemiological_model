"""Script to execute vector borne disease models simulation through a series of parameter values
Instantiates VBDM (Vector Borne Disease Model)
given paramters from local_param_config.yaml and config.yaml configuration files.

    Typical usage example:

    python models_params.py -c <absolute path to config.yaml> -d <disease name
    [dengue][wnv]> -l <rename labels> -g <generate parameter data>
    -p <absolute path to parameter data file, only if no -g>
"""
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser_exp
from generate_param_values import param_data
import sys
import numpy as np
import pandas as pd
import yaml


def main():
    """Runs vector borne disease simulations.

    Instaniates each disease from configuration file
    """
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser_exp()
        args, unknown = parser.parse_known_args()

    disease_name = args.disease_name.lower()
    if args.generate_params == True:
        param_data(disease_name)
        param_values = pd.read_csv(f'parameters/{disease_name}_param_values.csv')
    elif args.generate_params == False:
        param_values = pd.read_csv(args.param_data_file)
    
    for k in range(0,len(param_values.index)):
        params_k = param_values.iloc[k,].to_dict()
            
        if disease_name == 'dengue':
             disease = DengueSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
        elif disease_name == 'wnv':
            disease = WNVSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
        disease.logger.info(disease)
        disease.run_model(disease_name)
        disease.save_output(disease_name, args.sim_labels, param_values)
        disease.logger.info('SUCCESS')


if __name__ == "__main__":
    main()
