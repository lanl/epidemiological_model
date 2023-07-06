 """Main script to execute vector borne disease models simulation
Run from run_main.sh script. Instantiates VBDM (Vector Borne Disease Model)
given paramters from config.yaml configuration file.

    Typical usage example:

    python models_main.py -c <absolute path to config.yaml> -d <disease name
    [dengue][wnv]> -m <absolute path to mosquito data directory>
"""
import sys
import os

import numpy as np
import pandas as pd

from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from LLM import LogisticLinkModel
from utils import create_arg_parser_LLM


def main():
    """Runs vector borne disease simulations.

    Instaniates each disease from configuration file
    """
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser_LLM()
        args, unknown = parser.parse_known_args()
    
    #im thinking i'm going to want the mosquito directory as an argument
    #sort the filenames so the naming tag will match up with the LLM output
    mosq_data = sorted(os.listdir(args.mosquito_dir))
    #pulls the HPU tag
    naming_tag = [k[12:-4] for k in mosq_data]
    mosq_data_paths = [os.path.join(args.mosquito_dir, k) for k in mosq_data]
    
    for i in range(0, len(mosq_data_paths)):
        print(i)
        mosq = LogisticLinkModel(args.config_file, mosq_data_paths[i])
        #mosq = LogisticLinkModel('config/run_LLM_human.yaml', mosq_data_paths[i])
        mosq.run_model()
        mosq.save_output(naming_tag[i])
    
    #sort the filenames so the file name order will match the anming tag order
    LLM_output = sorted([k for k in os.listdir(mosq.config_dict['LLM_OUTPUT_DIR']) if '.csv' in k])
    disease_name = args.disease_name.lower()
    file_name_list = list()
    for i in range(0, len(LLM_output)):
        LLM_data = pd.read_csv(os.path.join(mosq.config_dict['LLM_OUTPUT_DIR'], LLM_output[i]))
        date_df = LLM_data[['start_date', 'end_date', 'year']]
        param_data = LLM_data[['r', 'r_s', 'K', 'K_s', 'Sv']]
        data_list = list()
        new_initial_states = dict()
        #set output path and add to file name list
        output_file = f'{disease_name}_model_output_{naming_tag[i]}.csv'
        file_name_list.append(output_file)
        for k in range(0,len(param_data.index)):
            params_k = param_data.iloc[k,].to_dict()
            dates_k = date_df.iloc[k,].to_dict()
            date_vec = pd.date_range(start=f"{dates_k['year']}-{dates_k['start_date']}", end = f"{dates_k['year'] + 1}-{dates_k['end_date']}")
            
            if len(new_initial_states) > 0:
                del params_k['Sv']
                params_k = {**params_k, **new_initial_states}
                
#             if disease_name == 'dengue':
#                  disease = DengueSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
            #elif disease_name == 'wnv':
            if disease_name == 'wnv':
                #disease = WNVSEIRModel.param_dict(config_file = 'config/run_LLM_human.yaml', param_dict = params_k)
                disease = WNVSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
            
            #reset the duration based on the dates + year being run
            disease.config_dict['DURATION'] = len(date_vec) - 1
            disease.t = (0, disease.config_dict['DURATION'])
            #need to add the +1 to get the correct step size
            disease.t_eval = np.linspace(0, disease.config_dict['DURATION'], disease.config_dict['DURATION']*disease.config_dict['RESOLUTION'] + 1)
            
            #disease.logger.info(disease)
            disease.run_model(disease_name, verbose = False)
            #add if statement so that if the solver breaks we just fill the missing bottom of the data with nans
            if len(disease.df) != len(date_vec):
                empty = np.empty((len(date_vec) - len(disease.df), len(disease.df.columns)))
                empty[:] = np.nan
                disease.df = pd.concat([disease.df, pd.DataFrame(empty, columns = disease.df.columns)], axis = 0).reset_index()
            #add this back in later
            #if disease_name == 'wnv':
                #pull_new_states_df = disease.df.drop('Ih', axis = 1)
                #new_initial_states = pull_new_states_df.iloc[-1].to_dict()
            disease.df['Date'] = date_vec
            data_list.append(disease.df)
            disease.logger.info(f'{naming_tag[i]}, {dates_k["year"]} success')
            
        disease.logger.info(f'{naming_tag[i]} SUCCESS')
            
        final_df = pd.concat(data_list).rename(columns = disease.state_names_order)
        output_path = os.path.join(disease.config_dict['OUTPUT_DIR'], output_file) 
        final_df.to_csv(output_path, index=False)
   
    #final success statement - log it and print into the console
    file_name_list_add = [os.path.join(disease.config_dict['OUTPUT_DIR'], k) for k in file_name_list]
    success_vec = [os.path.isfile(k) for k in file_name_list_add]
    if sum(success_vec) == len(mosq_data_paths):
        disease.logger.info(f'ALL HPU RUNS SUCCESSFUL')
        print(f'ALL HPU RUNS SUCCESSFUL')
    else:
        false_index = [i for i in range(0, len(success_vec)) if success_vec[i] == False]
        fail_hpus = [naming_tag[i] for i in false_index]
        log_fail = ' and '.join(fail_hpus)
        disease.logger.info(f'Run FAILED for the following HPUs: {log_fail}')
        print(f'Run FAILED for the following HPUs: {log_fail}')
        


if __name__ == "__main__":
    main()
