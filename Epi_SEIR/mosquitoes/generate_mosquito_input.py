import pandas as pd
import random
import numpy as np
import yaml

with open('/Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml', 'r') as in_file:
    n = yaml.safe_load(in_file)['DURATION']

arr = np.random.randint(100, 1500, size=n)

df = pd.DataFrame({'Sv': arr})

df.to_csv('/Users/jkeithley/Documents/CIMMID/human/dengue_model/epi_seir/mosquitoes/mosq.csv')
