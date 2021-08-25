"""Script to execute vector borne disease models simulation through a series of parameter values
Instantiates VBDM (Vector Borne Disease Model)
given paramters from local_param_config.yaml and config.yaml configuration files.

    Typical usage example:

    python models_params.py -c <absolute path to config.yaml> -d <disease name
    [dengue][wnv]> -l <rename labels [T][F]> -g <generate parameter data [T][F]>
    -p <absolute path to parameter data file, only if -g is [F]>
"""
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser_exp
from generate_param_values import param_data
import sys
import numpy as np
import pandas as pd
import yaml

# TODO remove this comment

def main():
    """Runs vector borne disease simulations.

    Instaniates each disease from configuration file
    """
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser_exp()
        args, unknown = parser.parse_known_args()

    disease_name = args.disease_name.lower()
    if args.generate_params == 'T':
        param_data(disease_name)
        param_values = pd.read_csv(f'parameters/{disease_name}_param_values.csv')
    elif args.generate_params == 'F':
        param_values = pd.read_csv(args.param_data_file)
    
    for k in range(0,len(param_values.index)):
        with open(args.config_file) as f:
             list_doc = yaml.safe_load(f)
        
        for n in param_values.columns:
            list_doc[disease_name.upper()]['PARAMETERS'][n] = param_values.iloc[k,param_values.columns.get_loc(n)].item()
        
        #saving as a different file so I don't keep ruining the format of the normal one
        with open('config/local_test_config_change.yaml', "w") as f:
            yaml.safe_dump(list_doc, f, default_flow_style=False)
            
        if disease_name == 'dengue':
            disease = DengueSEIRModel('config/local_test_config_change.yaml')
        elif disease_name == 'wnv':
            disease = WNVSEIRModel('config/local_test_config_change.yaml')
        disease.logger.info(disease)
        disease.run_model(disease_name)
        disease.save_output(disease_name, param_values, args.sim_labels)
        disease.logger.info('SUCCESS')


if __name__ == "__main__":
    main()
