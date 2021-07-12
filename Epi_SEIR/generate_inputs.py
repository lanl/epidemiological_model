import pandas as pd
#import random
import numpy as np
import yaml
#import pyarrow as pa
#import pyarrow.parquet as pq
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default='config/config.yaml')
args = parser.parse_args()

# generate mosquito population input
with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)

n = config_dict['DURATION']
output_path = config_dict['MOSQUITOES_FILE_PATH']

arr = np.random.randint(4000, 4001, size=n)
df = pd.DataFrame({'Sv': arr})
df.to_csv(output_path)
#pq.write_table(pa.Table.from_pandas(df), os.path.join(output_path,
#                                                          'mosq.parquet'))

# generate other initial states
#output_path = '/Users/jkeithley/Documents/CIMMID/human/dengue_model/epi_seir/initial_states_input'

#keys = ['Sh', 'Eh', 'Iha', 'Ihs', 'Rh', 'Ev', 'Iv']
#arr = np.random.randint(10, 100, size=len(keys))
#states = pd.DataFrame([dict(zip(keys, arr))])
#states.to_csv(os.path.join(output_path, 'initial_states.csv'), index=False)
#pq.write_table(pa.Table.from_pandas(states), os.path.join(output_path,
#                                                          'initial_states.parquet'))
