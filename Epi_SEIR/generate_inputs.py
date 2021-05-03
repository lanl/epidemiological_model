import pandas as pd
import random
import numpy as np
import yaml
import pyarrow as pa
import pyarrow.parquet as pq
import os

# generate mosquito population input
with open('/Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml', 'r') as in_file:
    n = yaml.safe_load(in_file)['DURATION']
output_path = '/Users/jkeithley/Documents/CIMMID/human/dengue_model/epi_seir/mosquitoes'

arr = np.random.randint(100, 1500, size=n)
df = pd.DataFrame({'Sv': arr})
df.to_csv(os.path.join(output_path, 'mosq.csv'))
pq.write_table(pa.Table.from_pandas(df), os.path.join(output_path,
                                                          'mosq.parquet'))

# generate other initial states
#output_path = '/Users/jkeithley/Documents/CIMMID/human/dengue_model/epi_seir/initial_states_input'

#keys = ['Sh', 'Eh', 'Iha', 'Ihs', 'Rh', 'Ev', 'Iv']
#arr = np.random.randint(10, 100, size=len(keys))
#states = pd.DataFrame([dict(zip(keys, arr))])
#states.to_csv(os.path.join(output_path, 'initial_states.csv'), index=False)
#pq.write_table(pa.Table.from_pandas(states), os.path.join(output_path,
#                                                          'initial_states.parquet'))
