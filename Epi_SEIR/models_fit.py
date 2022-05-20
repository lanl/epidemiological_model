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
    disease.fit_constants(disease_name)
    disease.proflike(disease_name)
    disease.plot_proflike()
    if disease.calc_ci_bool == True:
        disease.calc_ci()
    disease.save_fit_output(disease_name)
    disease.logger.info('SUCCESS')
    
    if args.run_fit_model == True:
        fit_params_data = pd.read_csv(f'param_fit_output/{disease_name}_parameter_values.csv')
        fit_params = fit_params_data.iloc[0,].to_dict()
        if 'dispersion' in list(fit_params.keys()):
            del fit_params['dispersion']
        fit_disease = DengueSEIRModel.param_dict(config_file = args.config_file, param_dict =  fit_params)
        if 'Dh' in fit_disease.fit_data_compartment:
            fit_disease.run_model(disease_name, False, True)
        else:
            fit_disease.run_model(disease_name, False)
        fit_disease.save_output(disease_name, True, fit_params_data)

    if args.figure == True:
        fit_disease.df['Time'] = fit_disease.t_eval
        model_df = fit_disease.df.copy()
        for i in range(0, len(disease.fit_data)):
            data = disease.fit_data[i]
            compartment = disease.fit_data_compartment[i]
            model_df = model_df[model_df[compartment].isna() == False].reset_index(drop = True)
            if disease.fit_data_res[i] == 'weekly':
                data['Time'] = fit_disease.t_eval[::7]
                if compartment == 'Dh':
                    time = list(model_df['Time'][::7])
                    out = model_df.groupby(model_df.index // 7).sum()[compartment]
                    model_df = pd.DataFrame({'Time':time, compartment:out})
                    
            elif disease.fit_data_res[i] == 'daily':
                data['Time'] = fit_disease.t_eval
                       
            plt.plot(data['Time'], data[compartment], 'ro', label=f'{compartment} data')
            plt.plot(model_df['Time'], model_df[compartment], 'b-', label= f'{compartment} fit')
            plt.legend(loc='best')
            if args.figure == True:
                plt.show()
                plt.close()
            if args.save_figure == True:
                plt.savefig(f'plots/{disease_name}_fit_{compartment}.png')
                plt.close()

if __name__ == "__main__":
    fit()
