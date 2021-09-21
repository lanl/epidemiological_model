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
    
    data_list = list()
    for k in range(0,len(param_values.index)):
        params_k = param_values.iloc[k,].to_dict()
            
        if disease_name == 'dengue':
             disease = DengueSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
        elif disease_name == 'wnv':
            disease = WNVSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
        disease.logger.info(disease)
        disease.run_model(disease_name)
        disease.logger.info('SUCCESS')
        
       # out = disease.get_data()
       # final_size = out['Recovered Humans'][len(out['Recovered Humans']) -1]
       # R0 = disease.params['nu_v'] / (disease.params['mu_v'] + disease.params['nu_v']) * \
        #    (((disease.params['a_v']**2) * disease.params['beta_v'] * disease.params['beta_h'] * disease.params['K_v'] * final_size) / \
         #   (disease.params['gamma_h'] * disease.params['mu_v'] * (100000**2)))
       # data_list.append(pd.DataFrame(np.array([[R0]]).T, columns = ['R0'], index = [k]))
        
        disease.save_output(disease_name, args.sim_labels, param_values)
        
   # final_data = pd.concat(data_list)
   # final_data.to_csv('/Users/mbarnard/Documents/Dengue_Paper/global_sense_9_21_2021.csv', index=False)

if __name__ == "__main__":
    main()
