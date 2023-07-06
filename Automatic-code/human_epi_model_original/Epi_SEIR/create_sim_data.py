###OPEN THIS FILE IN THE REPO OR IT WONT WORK###
from dengue import DengueSEIRModel
from utils import create_arg_parser_fit
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml

#run model and turn output into df
disease_name = 'dengue'
self = DengueSEIRModel('config/local_test_config.yaml')
self.run_model(disease_name)
model_df = pd.DataFrame(dict(zip(list(self.initial_states.keys()), self.model_output.T)))

#get weekly cases and then simulate data based on a negative binomial distribution
week_out = model_df.iloc[::7, :].reset_index()
week_out['Dh'] = week_out['Ch'].diff().fillna(1)
week_out = week_out[week_out['Dh'].isna() == False]
#the 0.24 below is the dispersion value, so you can decide what to change it to
sigma_squared = week_out['Dh'] + 0.24*(week_out['Dh']**2)
p =  week_out['Dh']  / sigma_squared 
n =  week_out['Dh'] **2 / (sigma_squared -  week_out['Dh'])
#I believe the 4.0 below was hard coded, should be changed to np.unique(n) or whatever stable value n takes
sim_noise = [np.random.negative_binomial(n[0],k) for k in p]
        
#plot to see what the simulated data looks like
week_out['Time'] = self.t_eval[::7]
#week_out['Time'] = self.t_eval[::7][1:len(self.t_eval)]
plt.plot(week_out['Time'], sim_noise, 'ro', label=f'sim data')
plt.plot(week_out['Time'], week_out['Dh'], 'b-', label= f'model')
plt.legend(loc='best')

sim_noise = pd.DataFrame({'Dh': sim_noise})
#then you can save the dataset to the fit_data folder