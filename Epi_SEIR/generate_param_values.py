"""Script to generate parameter value combinations to sequence through for models_params.py
Generates parameter values based on values input in local_param_config.yaml
Outputs a csv of all possible parameter combinations in Epi_SEIR/parameters

    Typical usage example:
    
    param_data(<disease_name>)
    
    Output: 
    
    parameters/<disease_name>_param_values.csv
    
"""

import pandas as pd
import numpy as np
import yaml
import os


def param_data(disease_name):
    """Generates input data for sequencing through parameter values
    """
    with open('config/local_param_config.yaml', 'r') as in_file:
        config_dict = yaml.safe_load(in_file)
    
    keys = list(config_dict[disease_name.upper()].keys())
    out = []
    mult = 1
    for i in range(0,len(keys)):
        start = config_dict[disease_name.upper()][keys[i]]['start']
        stop = config_dict[disease_name.upper()][keys[i]]['stop']
        step = config_dict[disease_name.upper()][keys[i]]['step']
        out.append(np.arange(start, stop, step))
        mult = mult * len(out[i])
    
    #get df of all possible combinations of the parameters       
    df = pd.DataFrame(np.array(np.meshgrid(*out)).reshape(len(keys), mult).T, columns = keys) 
    output_path = os.path.join(config_dict['OUTPUT_DIR'], f'{disease_name}_param_values.csv')
    df.to_csv(output_path, index = False)
