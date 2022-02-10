from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser_fit
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml


disease_name = 'dengue'
self = DengueSEIRModel('config/local_test_config.yaml')
self.run_model(disease_name)

model_df = pd.DataFrame(dict(zip(list(self.initial_states.keys()), self.model_output.T)))

week_out = model_df.iloc[::7, :].reset_index()
week_out['Dh'] = week_out['Ch'].diff().fillna(1)
sim_noise = np.round(week_out['Dh'] + abs(np.random.normal(scale=100, size=len(week_out['Dh']))))
#sim_noise = [np.random.poisson(k) for k in week_out['Dh']]

week_out['Time'] = self.t_eval[::7]

plt.plot(week_out['Time'], sim_noise, 'ro', label=f'sim data')
plt.plot(week_out['Time'], week_out['Dh'], 'b-', label= f'model')
plt.legend(loc='best')

#sim_noise.to_csv('fit_data/sim_data.csv')


#made based on the following params
#     nu_h: 0.164 
#     nu_v: 0.164 
#     gamma_h: 0.167 
#     mu_v: .0714 
#     beta_h: 0.33
#     beta_v: 0.33
#     a_v: 0.75
#     r_v: .2286
#     #r_v: .5
#     #K_v: 200000
#     K_v: 4028296
    
#      Sh:
#     name: 'Susceptible Humans'
#     value: 6748000
#     position: 0
#   Eh:
#     name: 'Exposed Humans'
#     value: 0
#     position: 1
#   Ih:
#     name: 'Infected Humans'
#     value: 0
#     position: 2
#   Rh:
#     name: 'Recovered Humans'
#     value: 0
#     position: 3
#   Ch:
#     name: 'Cummulative Infected Humans' #added for fitting purposes
#     value: 0
#     position: 4
#   Sv:
#     name: 'Susceptible Vectors'
#     value: 2513848
#     position: 5
#   Ev:
#     name: 'Exposed Vectors'
#     value: 0
#     position: 6
#   Iv:
#     name: 'Infected Vectors'
#     #value: 500
#     value: 94
#     position: 7



