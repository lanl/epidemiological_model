from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser_fit
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
import numdifftools
#import loess
from scipy.interpolate import interp1d
import statsmodels.api as sm
import scipy

disease_name = 'dengue'
self = DengueSEIRModel('config/local_test_config.yaml')
self.fit_constants()
self.fit_out


fit_disease = DengueSEIRModel.param_dict('config/local_test_config.yaml', param_dict =  self.fit_out.params)
fit_disease = DengueSEIRModel.param_dict('config/local_test_config.yaml', param_dict =  {'a_v': 0.74932285})
fit_disease.run_model(disease_name)


model_df = pd.DataFrame(dict(zip(list(fit_disease.initial_states.keys()), fit_disease.model_output.T)))

week_out = model_df.iloc[::7, :].reset_index()
week_out['Dh'] = week_out['Ch'].diff().fillna(0)

#data = pd.read_csv('fit_data/rdj_2010_2011_season_peak.csv')
data = pd.read_csv('fit_data/sim_data3.csv')

week_out['Time'] = fit_disease.t_eval[::7]


plt.plot(week_out['Time'], data['Dh'], 'ro', label=f'data')
plt.plot(week_out['Time'], week_out['Dh'], 'b-', label= f'fit')
plt.title('Sim a_v = 0.75, fit a_v = 0.7493, Sim disp = 0.25, fit disp = 0.309')
plt.legend(loc='best')
#plt.savefig('plots/fit_test/nbinom/fit2_simdata3.png')

#This is incredibly variable for some of the parameters - specifically Sv is not behaving well
test = self.proflike()
#Sv parameter fit highly dependent on initial value...with a guess of around 2,000,000 but we had a higher NLL than 3,000,000
plt.plot(test[1]['a_v'], test[1]['nll'], 'ro', label = 'NLL')
plt.plot(test[1]['a_v'], [test[0]]*len(test[1]), 'b-', label = 'Threshold')
#plt.plot(test[1]['a_v'], lowess_y, 'p-', label = 'lowess fit')
plt.legend(loc='best')
plt.title('NLL of Biting rate a_v')
#plt.savefig('plots/fit_test/a_v_0_001_simdata_no_round.png')
#plt.savefig('plots/fit_test/a_v_0_0003_riodata_no_round.png')
#plt.savefig('plots/fit_test/a_v_0_05_simdata2_nbinom.png')
plt.savefig('plots/fit_test/nbinom/fit2_a_v_simdata3.png')

plt.plot(test[2]['dispersion'], test[2]['nll'], 'ro', label = 'NLL')
plt.plot(test[2]['dispersion'], [test[0]]*len(test[2]), 'b-', label = 'Threshold')
plt.legend(loc='best')
plt.title('NLL of Dispersion')
plt.savefig('plots/fit_test/nbinom/fit2_dispersion_simdata3.png')

#try a loess fit
#frac 1 makes it closer to the actual points, which is what I think we want
lowess = sm.nonparametric.lowess(test[1]['nll'], test[1]['a_v'], frac = .1)
lowess_x = list(zip(*lowess))[0]
lowess_y = list(zip(*lowess))[1]

f = interp1d(lowess_x, lowess_y, bounds_error=False)
f_thresh = lambda x: f(x) - test[0]

#five makes it a slightly arbitrary guess but not terrible
lb = scipy.optimize.newton(f_thresh, test[1]['a_v'][5])
ub = scipy.optimize.newton(f_thresh, test[1]['a_v'][len(test[1]) - 5])


plt.plot(test[2]['nu_h'], test[2]['nll'], 'ro', label = 'NLL')
plt.plot(test[2]['nu_h'], [test[0]]*len(test[2]), 'b-', label = 'Threshold')
plt.legend(loc='best')
plt.title('Biting rate a_v')
#plt.savefig('plots/fit_test/a_v_0_01.png')

###Test why simulated data + nbinom is being weird
#I can't tell why this is happening
self = DengueSEIRModel('config/local_test_config.yaml')
params_fit = self.init_fit_parameters()
test = obj()
resid = test[0] - test[1]

a7_mse = np.mean(resid**2)
a75_mse = np.mean(resid**2)
a9_mse = np.mean(resid**2)

time = np.arange(0, len(test[0]), 1)
plt.plot(time, test[0], 'ro', label=f'data')
plt.plot(time, test[1], 'b-', label= f'fit')
plt.legend(loc='best')

mod_out = test[1]
sigma_squared = mod_out + self.dispersion*(mod_out**2)
p = mod_out / sigma_squared #[k/sigma**2 for k in mod_out]
n = mod_out**2 / (sigma_squared - mod_out)

probs_a7 = nbinom.pmf(np.round(test[0]),n, np.array(p))
log_probs_a7 = np.log(probs_a7)
nll_a7 = -sum(log_probs_a7)
nll_fudge_a7 = -sum(np.log(probs_a7 + 1e-323))

probs_a75 = nbinom.pmf(np.round(test[0]),n, np.array(p))
log_probs_a75 = np.log(probs_a75)
nll_a75 = -sum(log_probs_a75)
nll_fudge_a75 = -sum(np.log(probs_a75 + 1e-323))

probs_a9 = nbinom.pmf(np.round(test[0]),n, np.array(p))
log_probs_a9 = np.log(probs_a9)
nll_a9 = -sum(log_probs_a9)
nll_fudge_a9 = -sum(np.log(probs_a9 + 1e-323))
