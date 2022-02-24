from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser_fit
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
import numdifftools

disease_name = 'dengue'
self = DengueSEIRModel('config/local_test_config.yaml')
self.fit_constants()


fit_disease = DengueSEIRModel.param_dict('config/local_test_config.yaml', param_dict =  self.fit_out.params)
fit_disease.run_model(disease_name)


model_df = pd.DataFrame(dict(zip(list(fit_disease.initial_states.keys()), fit_disease.model_output.T)))

week_out = model_df.iloc[::7, :].reset_index()
week_out['Dh'] = week_out['Ch'].diff().fillna(0)

#data = pd.read_csv('fit_data/rdj_2010_2011_season_peak.csv')
data = pd.read_csv('fit_data/sim_data2.csv')

week_out['Time'] = fit_disease.t_eval[::7]

plt.plot(week_out['Time'], data['Dh'], 'ro', label=f'data')
plt.plot(week_out['Time'], week_out['Dh'], 'b-', label= f'fit')
plt.legend(loc='best')

#This is incredibly variable for some of the parameters - specifically Sv is not behaving well
test = self.proflike()
#Sv parameter fit highly dependent on initial value...with a guess of around 2,000,000 but we had a higher NLL than 3,000,000
plt.plot(test[1]['a_v'], test[1]['nll'], 'ro', label = 'NLL')
plt.plot(test[1]['a_v'], [test[0]]*len(test[1]), 'b-', label = 'Threshold')
plt.legend(loc='best')
plt.title('Biting rate a_v')
#plt.savefig('plots/fit_test/a_v_0_001_simdata_no_round.png')
#plt.savefig('plots/fit_test/a_v_0_0003_riodata_no_round.png')
#plt.savefig('plots/fit_test/a_v_0_05_simdata2_nbinom.png')


plt.plot(test[2]['nu_h'], test[2]['nll'], 'ro', label = 'NLL')
plt.plot(test[2]['nu_h'], [test[0]]*len(test[2]), 'b-', label = 'Threshold')
plt.legend(loc='best')
plt.title('Biting rate a_v')
#plt.savefig('plots/fit_test/a_v_0_01.png')


# nll = list()
# for p in param_seq:
#                 if k in list(self.params.keys()):
#                     self.params[k] = p

#                     #just need to put something in the params category, doesn't really matter what 
#                 nll.append(fit_objective(self.params))
#                 print(self.params[k])
                
# plt.plot(test[3]['K_v'], test[3]['nll'], 'ro', label = 'NLL')
# plt.plot(test[3]['K_v'], [test[0]]*len(test[3]), 'b-', label = 'Threshold')
# plt.legend(loc='best')


#Results without fitting Sv

plt.plot(test[1]['Iv'], test[1]['nll'], 'ro', label = 'NLL')
plt.plot(test[1]['Iv'], [test[0]]*len(test[1]), 'b-', label = 'Threshold')
plt.legend(loc='best')


plt.plot(test[2]['K_v'], test[2]['nll'], 'ro', label = 'NLL')
plt.plot(test[2]['K_v'], [test[0]]*len(test[2]), 'b-', label = 'Threshold')
plt.legend(loc='best')