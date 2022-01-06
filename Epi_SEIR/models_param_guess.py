"""Script to execute vector borne disease models simulation through a series of parameter values
Instantiates VBDM (Vector Borne Disease Model)
given paramters from local_param_config.yaml and config.yaml configuration files.

    Typical usage example:

    python models_param_guess.py -c <absolute path to config.yaml> -d <disease name
    [dengue][wnv]>  -g <generate parameter data> -p <absolute path to parameter data file, only if no -g> 
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
    """Runs vector borne disease simulations through a series of parameters, calculates cost function, determines parameter combination that minimizes cost

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
    
    min_value_list = list()
    for k in range(0,len(param_values.index)):
        params_k = param_values.iloc[k,].to_dict()
            
        if disease_name == 'dengue':
            disease = DengueSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
            #disease = DengueSEIRModel.param_dict(config_file = config_file, param_dict =  params_k)
        elif disease_name == 'wnv':
            disease = WNVSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
        #disease.logger.info(disease)
        if disease.fit_method == 'res':
            min_value = sum(abs(disease.fit_objective({**disease.params, **disease.initial_states})))
        else:
            min_value = disease.fit_objective({**disease.params, **disease.initial_states})
        min_value_list.append(min_value)
       # disease.logger.info('SUCCESS')
        
    param_values['min_value'] = min_value_list
    best_param_values = param_values[param_values['min_value'] == param_values['min_value'].min()]
    best_param_values.drop('min_value', axis = 1).to_csv(f'param_fit_output/{disease_name}_best_param_guess.csv', index = False)
    

        
if __name__ == "__main__":
    main()
