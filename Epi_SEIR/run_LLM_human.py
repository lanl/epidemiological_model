"""Main script to execute vector borne disease models simulation
Run from run_main.sh script. Instantiates VBDM (Vector Borne Disease Model)
given paramters from config.yaml configuration file.

    Typical usage example:

    python models_main.py -c <absolute path to config.yaml> -d <disease name
    [dengue][wnv]> -f <show model output figures> -sf <save model output figures>
"""
import sys
import os

import pandas as pd

from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from LLM import LogisticLinkModel
from utils import create_arg_parser


def main():
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
    
    #im thinking i'm going to want the mosquito directory as an argument
    mosq_data = sorted(os.listdir(disease.config_dict['MOSQUITOES_DIR']))
    #pulls the HPU tag
    naming_tag = [k[12:-4] for k in mosq_data]
    mosq_data_paths = [os.path.join(disease.config_dict['MOSQUITOES_DIR'], k) for k in mosq_data]
    
    for i in range(0, len(mosq_data_paths)):
        #mosq = LogisticLinkModel(args.config_file, mosq_data_paths[i])
        mosq = LogisticLinkModel('config/run_LLM_human.yaml', mosq_data_paths[i])
        mosq.run_model()
        mosq.save_output(naming_tag[i])
    
    LLM_output = sorted([k for k in os.listdir(disease.config_dict['LLM_OUTPUT_DIR']) if '.csv' in k])
    
    for i in range(0, len(LLM_output)):
        param_data = pd.read_csv(os.path.join(disease.config_dict['LLM_OUTPUT_DIR'], LLM_output[i]))
        data_list = list()
        new_initial_states = dict()
        for k in range(0,len(param_values.index)):
            params_k = param_values.iloc[k,].to_dict()
            
            if len(new_initial_states) > 0:
                del params_k['Sv']
                params_k = {**params_k, **new_initial_states}
                
            if disease_name == 'dengue':
                 disease = DengueSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
            elif disease_name == 'wnv':
                disease = WNVSEIRModel.param_dict(config_file = args.config_file, param_dict =  params_k)
            disease.logger.info(disease)
            disease.run_model(disease_name, verbose = False)
            #we need to see if we want to do with for Sv, I'm guessing we do
            new_initial_states = disease.df.iloc[-1].to_dict()
            data_list.append(self.df)
            disease.logger.info('SUCCESS')
            
        final_df = pd.concat(data_list).rename(columns = disease.self.state_names_order)
        output_path = os.path.join(disease.config_dict['OUTPUT_DIR'],
                                           f'{disease_name}_model_output_{naming_tag[i]}.csv')
        final_df.to_csv(output_path, index=False)
        


if __name__ == "__main__":
    main()
