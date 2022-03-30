"""Script to fit vector borne disease models to data
Instantiates VBDM (Vector Borne Disease Model)
and fitting process given paramters from local_param_config.yaml and config.yaml configuration files.

    Typical usage example:

    python models_fit.py -c <absolute path to config.yaml> -d <disease name> -p <if you want to use the best guess params> -rm <if you want to run model with fit params>
    -f <if you want to plot model output against data> -sf <if you want to save that plot>
"""
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser_fit
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
import numdifftools


def fit():
    """Runs vector borne disease simulations.

    Instaniates each disease from configuration file
    """
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser_fit()
        args, unknown = parser.parse_known_args()

    disease_name = args.disease_name.lower()

    if disease_name == 'dengue':
        if args.best_guess_params == True:
            disease = DengueSEIRModel.param_dict(args.config_file, pd.read_csv(f'param_fit_output/{disease_name}_best_param_guess.csv').iloc[0,:].to_dict())
        else:
            disease = DengueSEIRModel(args.config_file)
    elif disease_name == 'wnv':
        if args.best_guess_params == True:
            disease = WNVSEIRModel(args.config_file, pd.read_csv(f'param_fit_output/{disease_name}_best_param_guess.csv').to_dict())
        else:
            disease = WNVSEIRModel(args.config_file)

    disease.logger.info(disease)
    disease.fit_constants()
    disease.proflike()
    if disease.calc_ci_bool == True:
        disease.calc_ci()
    disease.save_fit_output(disease_name)
    disease.logger.info('SUCCESS')
    
    if args.run_fit_model == True:
        fit_params = pd.read_csv(f'param_fit_output/{disease_name}_parameter_values.csv').iloc[0,].to_dict()
        fit_disease = DengueSEIRModel.param_dict(config_file = args.config_file, param_dict =  fit_params)
        fit_disease.logger.info(fit_disease)
        fit_disease.run_model(disease_name)
        disease.logger.info('SUCCESS')

    if args.figure == True:
        model_df = pd.DataFrame(dict(zip(list(fit_disease.initial_states.keys()), fit_disease.model_output.T)))
        model_df['Time'] = fit_disease.t_eval
        for i in range(0, len(disease.fit_data)):
            data = disease.fit_data[i]
            if disease.fit_data_res[i] == 'weekly':
                data['Time'] = fit_disease.t_eval[::7]
            elif disease.fit_data_res[i] == 'daily':
                data['Time'] = fit_disease.t_eval
                       
            compartment = disease.fit_data_compartment[i]
            plt.plot(data['Time'], data[compartment], 'ro', label=f'{compartment} data')
            plt.plot(model_df['Time'], model_df[compartment], 'b-', label= f'{compartment} fit')
            plt.legend(loc='best')
            if args.save_figure == False:
                plt.show()
                plt.close()
            if args.save_figure == True:
                plt.savefig(f'plots/{disease_name}_fit_{compartment}.png')
                plt.close()

if __name__ == "__main__":
    fit()
