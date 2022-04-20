"""Script to execute vector borne disease models simulation through a series of parameter values
Instantiates VBDM (Vector Borne Disease Model)
given paramters from local_param_config.yaml and config.yaml configuration files.

    Typical usage example:

    python models_param_guess.py -c <absolute path to config.yaml> -d <disease name
    [dengue][wnv]>  -g < generate parameter data> -p <absolute path to parameter data file, only if no -g> 
    -f <include if a plot of the data and fitted model with the best guess params is desired>
"""
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser_exp
from generate_param_values import param_data
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
            disease = DengueSEIRModel(args.config_file)
        elif disease_name == 'wnv':
            disease = WNVSEIRModel(args.config_file)
        min_value = disease.fit_objective(params_k, disease_name)
        min_value_list.append(min_value)
        
    param_values['min_value'] = min_value_list
    best_param_values = param_values[param_values['min_value'] == param_values['min_value'].min()]
    best_param_values.drop('min_value', axis = 1).to_csv(f'param_fit_output/{disease_name}_best_param_guess.csv', index = False)
    
    if args.figure == True:
        best_params = best_param_values.drop('min_value', axis = 1).iloc[0,:].to_dict()
        if 'dispersion' in list(best_params.keys()):
            del best_params['dispersion']
            
        fit_disease = DengueSEIRModel.param_dict(config_file = args.config_file, param_dict =  best_params)
        
        if 'Dh' in fit_disease.fit_data_compartment:
            fit_disease.run_model(disease_name, False, True)
        else:
            fit_disease.run_model(disease_name, False)
            
        fit_disease.df['Time'] = fit_disease.t_eval
        model_df = fit_disease.df.copy()
        for i in range(0, len(disease.fit_data)):
            data = fit_disease.fit_data[i]
            compartment = fit_disease.fit_data_compartment[i]
            model_df = model_df[model_df[compartment].isna() == False].reset_index(drop = True)
            if fit_disease.fit_data_res[i] == 'weekly':
                data['Time'] = fit_disease.t_eval[::7]
                if compartment == 'Dh':
                    time = list(model_df['Time'][::7])
                    out = model_df.groupby(model_df.index // 7).sum()[compartment]
                    model_df = pd.DataFrame({'Time':time, compartment:out})
                    
            elif fit_disease.fit_data_res[i] == 'daily':
                data['Time'] = fit_disease.t_eval
                       
            plt.plot(data['Time'], data[compartment], 'ro', label=f'{compartment} data')
            plt.plot(model_df['Time'], model_df[compartment], 'b-', label= f'{compartment} fit')
            plt.legend(loc='best')
            plt.show()
            plt.close()
    
    

        
if __name__ == "__main__":
    main()
