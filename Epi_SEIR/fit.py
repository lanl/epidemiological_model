"""Fitting Model Class.

Contains class for specific diseases that inherit from the
class FitModel. Inherits from the VectorBorneDiseaseModel class.

"""
import os
from utils import timer
from abc import abstractmethod
#need numdifftools installed to estimate covariance and then get std errors https://lmfit.github.io/lmfit-py/fitting.html
import numpy as np
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

from scipy.integrate import solve_ivp
from lmfit import Parameters, minimize, fit_report
from scipy.stats import poisson
from scipy.stats import norm
from scipy.stats import nbinom
from scipy.stats import chi2
from scipy.optimize import newton
from scipy.interpolate import interp1d

import vbdm



class FitModel(vbdm.VectorBorneDiseaseModel):
    """Defines a general fitting procedure.

    Inherits from VectorBorneDiseaseModel.

    Attributes:
        config_dict: All configurations read in from config file.\n
        params: ODE system parameters.\n
        initial_states: initial states for population sizes.\n
        fit_params: parameter names and range values to fit.\n
        fit_data: list of csvs to fit the model to.\n

    """

    def __init__(self, config_file, disease_name, param_dict = None): 
        
        super().__init__(config_file, disease_name)

    
        if self.config_dict[disease_name]['FIT'] == True:
            self.fit_params = list(self.config_dict[disease_name]['FIT_PARAMS'].keys())
            self.fit_method = self.config_dict[disease_name]['FIT_METHOD']
            self.calc_ci_bool = self.config_dict[disease_name]['CALC_CI']
            self.dispersion = self.config_dict[disease_name]['DISPERSION']
            self.fit_params_range = self.config_dict[disease_name]['FIT_PARAMS']
            self.fit_data_res = self.config_dict[disease_name]['FIT_DATA']['res']
            self.fit_data_compartment = self.config_dict[disease_name]['FIT_DATA']['compartment']
            self.fit_data = []
            for k in self.config_dict[disease_name]['FIT_DATA']['PATH']:
                    self.fit_data.append(pd.read_csv(k))
                    
            self.error_check_resolution()
            self.error_check_compartment_names()
            self.error_check_data_and_output_length()
            

        
    @classmethod
    def param_dict(cls, config_file, disease_name, param_dict):
        """Inherit vbdm param_dict class method"""
        return super(FitModel, cls).param_dict(config_file = config_file, disease_name = disease_name, param_dict =  param_dict)

    @abstractmethod
    def model_func(self, t, y):
        pass
    
    def init_fit_parameters(self):
        #create the Parameter objects
            params_obj = Parameters()
            #add selected parameters and guesses into the framework
            param_keys = [i for i in self.fit_params if i in list(self.params.keys())]
            init_keys = [i for i in self.fit_params if i in list(self.initial_states.keys())]
            for k in param_keys:
                params_obj.add(k, value = self.params[k], min = self.fit_params_range[k]['min'], max = self.fit_params_range[k]['max'])
                #params_obj.add(k, self.params[k])
            for j in init_keys:
                params_obj.add(j, value = self.initial_states[j], min = self.fit_params_range[j]['min'], max = self.fit_params_range[j]['max'])
                #params_obj.add(j, self.initial_states[j])
            if 'dispersion' in self.fit_params:
                params_obj.add('dispersion', value = self.dispersion, min = self.fit_params_range['dispersion']['min'], max = self.fit_params_range['dispersion']['max'])
            return(params_obj)
    
    #adding this for the proflikelihood calculation
    def _init_proflik_parameters(self, lik_fitparams):
        #create the Parameter objects
            params_obj = Parameters()
            #add selected parameters and guesses into the framework
            param_keys = [i for i in lik_fitparams if i in list(self.params.keys())]
            init_keys = [i for i in lik_fitparams if i in list(self.initial_states.keys())]
            for k in param_keys:
                params_obj.add(k, lik_fitparams[k].value, min = lik_fitparams[k].value*.1, max = lik_fitparams[k].value*1.9)
                #params_obj.add(k, lik_fitparams[k].value)
            for j in init_keys:
                params_obj.add(j, lik_fitparams[j].value, min = lik_fitparams[j].value*.1, max = lik_fitparams[j].value*1.9)
                #params_obj.add(j, lik_fitparams[j].value)
            if 'dispersion' in lik_fitparams:
                params_obj.add('dispersion', lik_fitparams['dispersion'].value, min = lik_fitparams['dispersion'].value*.1, max = lik_fitparams['dispersion'].value*1.9)
            return(params_obj)
         
       
    def fit_run_model(self, params_fit, disease_name):
        param_keys = [i for i in params_fit if i in list(self.params.keys())]
        init_keys = [i for i in params_fit if i in list(self.initial_states.keys())]
        
        for k in param_keys:
            self.params[k] = params_fit[k]
        for j in init_keys:
            self.initial_states[j] = params_fit[j]
        if 'dispersion' in params_fit:
            self.dispersion = params_fit['dispersion']
        
        if 'Dh' in self.fit_data_compartment:
            self.run_model(disease_name, False, True)
        else:
            self.run_model(disease_name, False)
            
    
    def fit_objective(self, params_fit, disease_name):
        if self.fit_method == 'pois':
            data = np.empty(shape = 0)
            mod_out = np.empty(shape = 0)
            self.fit_run_model(params_fit, disease_name)
            for i in range(0,len(self.fit_data)):
                res = self.fit_data_res[i]
                compartment = self.fit_data_compartment[i]
                df = self.fit_data[i][compartment]
                #drop nas in the compartment - mainly to drop the first obs if we are fitting daily cases
                prep_df = self.df.copy()
                prep_df = prep_df[prep_df[compartment].isna() == False].reset_index(drop = True)
                if res == "weekly":
                    if compartment == 'Dh':
                        week_agg = prep_df.groupby(prep_df.index // 7).sum()
                        out = week_agg[compartment]
                    else:
                        week_out = prep_df.iloc[::7, :].reset_index()
                        out = week_out[compartment]
                elif res == "daily":
                    out = prep_df[compartment]
                #think sometime about weighting this based on the amount of data each source has
                data = np.concatenate((data,df))
                mod_out = np.concatenate((mod_out,out))
                #add a if negative barrier so we can maybe get away without not rounding the model output
                mod_out = np.where(mod_out < 0, 0, mod_out)
            #adding a fudge factor for the log currently, because getting all zeros due to terrible fit
            return -sum(np.log(poisson.pmf(np.round(data),mod_out) + 1e-323))
                
        elif self.fit_method == 'nbinom':
            data = np.empty(shape = 0)
            mod_out = np.empty(shape = 0)
            self.fit_run_model(params_fit, disease_name)
            for i in range(0,len(self.fit_data)):
                res = self.fit_data_res[i]
                compartment = self.fit_data_compartment[i]
                df = self.fit_data[i][compartment]
                #drop nas in the compartment - mainly to drop the first obs if we are fitting daily cases
                prep_df = self.df.copy()
                prep_df = prep_df[prep_df[compartment].isna() == False].reset_index(drop = True)
                if res == "weekly":
                    if compartment == 'Dh':
                        week_agg = prep_df.groupby(prep_df.index // 7).sum()
                        out = week_agg[compartment]
                    else:
                        week_out = prep_df.iloc[::7, :].reset_index()
                        out = week_out[compartment]
                elif res == "daily":
                    out = prep_df[compartment]
                data = np.concatenate((data,df))
                mod_out = np.concatenate((mod_out,out))
                #add a if negative barrier so we can maybe get away without not rounding the model output
                #and make the filler not zero so we don't encounter NaN values
                mod_out = np.where(mod_out <= 0, 1e-323, mod_out)
            #parameterization as seen in https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html 
            sigma_squared = mod_out + self.dispersion*(mod_out**2)
            p = mod_out / sigma_squared #[k/sigma**2 for k in mod_out]
            #add below so the model won't error out when mod_out and sigma squared are 1e-323
            div_vec = sigma_squared - mod_out
            if 0 in div_vec:
                div_vec = [x if x !=0 else 1e-232 for x in div_vec]
            n = mod_out**2 / div_vec
            #can't have n be zero otherwise nbinom errors out
            if 0 in n:
                n = [x if x !=0 else 1e-232 for x in n]
            return -sum(np.log(nbinom.pmf(np.round(data),n, np.array(p)) + 1e-323))
        
        elif self.fit_method == 'norm':
            data = np.empty([0])
            mod_out = np.empty([0])
            self.fit_run_model(params_fit, disease_name)
            for i in range(0,len(self.fit_data)):
                res = self.fit_data_res[i]
                compartment = self.fit_data_compartment[i]
                df = self.fit_data[i][compartment]
                #drop nas in the compartment - mainly to drop the first obs if we are fitting daily cases
                prep_df = self.df.copy()
                prep_df = prep_df[prep_df[compartment].isna() == False].reset_index(drop = True)
                if res == "weekly":
                    if compartment == 'Dh':
                        week_agg = prep_df.groupby(prep_df.index // 7).sum()
                        out = week_agg[compartment]
                    else:
                        week_out = prep_df.iloc[::7, :].reset_index()
                        out = week_out[compartment]
                elif res == "daily":
                    out = prep_df[compartment]
                data = np.concatenate((data,df))
                mod_out = np.concatenate((mod_out,out))
                #check with ethan if below is correct for a normaly distribution as well
#                 sigma_squared = mod_out + self.dispersion*(mod_out**2)
#                 sigma_squared = np.where(sigma_squared <= 0, 1e-232, sigma_squared)
#                 sigma = np.sqrt(sigma_squared)
                #actually will just do sigma = self.dispersion*mod_out for now
                sigma = self.dispersion*mod_out
                sigma = np.where(sigma <= 0, 1e-232, sigma)
            return -sum(np.log(norm.pdf(data,mod_out,sigma)))
        
    
    @timer 
    def fit_constants(self, disease_name):
        try:
            self.fit_out = minimize(self.fit_objective, self.init_fit_parameters(), kws = {'disease_name': disease_name}, method = 'nelder')
            if self.fit_out.success != True:
                self.logger.exeption('Minimization algorithm for fitting failed')
        except Exception as e:
            self.logger.exception('Exception occurred when running minimization for model fitting')
            raise e
    
    #not going to add the specifics of the numbers calculated within the range to the config file at the moment 
    #need to make sure below is within the min/max range listed in the config file
    def _calc_param_range(self, param, perc = .01, num = 10):
        param_seq = np.linspace((1-perc)*param, (1 + perc)*param, 2*num)
        return param_seq
        
    def proflike(self, disease_name):
        threshold = self.fit_objective(self.fit_out.params, disease_name) + chi2.ppf(.95, len(self.fit_params))/2
        self.df_list = list()
        self.df_list.append(threshold)
        for k in self.fit_params:
            #reset self.params to the fit values for the start of every run
            param_keys = [i for i in self.fit_params if i in list(self.params.keys())]
            init_keys = [i for i in self.fit_params if i in list(self.initial_states.keys())]
            for j in param_keys:
                self.params[j] = self.fit_out.params[j].value
            for h in init_keys:
                self.initial_states[h] = self.fit_out.params[h].value
            if 'dispersion' in self.fit_params:
                self.dispersion = self.fit_out.params['dispersion'].value
            
            #only want to fit the other parameters
            lik_fitparams = self.fit_out.params.copy()
            del lik_fitparams[k] 
            
            #set up parameters to sequence through and lists for results
            
            if k in list(self.params.keys()):
                param_seq = self._calc_param_range(self.params[k], perc = self.fit_params_range[k]['proflike'])
            elif k in list(self.initial_states.keys()):
                param_seq = self._calc_param_range(self.initial_states[k], perc = self.fit_params_range[k]['proflike'])
            elif k == 'dispersion':
                param_seq = self._calc_param_range(self.dispersion, perc = self.fit_params_range['dispersion']['proflike'])
                
            nll = list()
            fit_success = list()
            
            #calculate nll through the sequenced parameters
            for p in param_seq:
                if k in list(self.params.keys()):
                    self.params[k] = p
                elif k in list(self.initial_states.keys()):
                    self.initial_states[k] = p
                elif k == 'dispersion':
                    self.dispersion = p
                
                if len(self.fit_params) > 1:
                    fit = minimize(self.fit_objective, self._init_proflik_parameters(lik_fitparams), kws = {'disease_name': disease_name}, method = 'nelder')
                    #update guesses for the next round
                    for i in lik_fitparams:
                        lik_fitparams[i].value = fit.params[i].value
                    fit_success.append(fit.success)
                    nll.append(self.fit_objective(fit.params, disease_name))
                else:
                    #just need to put something in the params category, doesn't really matter what 
                    fit_success.append(np.nan)
                    nll.append(self.fit_objective(self.params, disease_name))
            
            df = pd.DataFrame({k:param_seq, 'nll': nll, 'success': fit_success})
            self.df_list.append(df)
        #return self.df_list
    
    def plot_proflike(self):
        for i in range(1, len(self.df_list)):
            plt.plot(self.df_list[i][self.df_list[i].columns[0]], self.df_list[i]['nll'], 'ro', label = 'NLL')
            plt.plot(self.df_list[i][self.df_list[i].columns[0]], [self.df_list[0]]*len(self.df_list[i]), 'b-', label = 'Threshold')
            plt.legend(loc='best')
            plt.title(f'NLL of {self.df_list[i].columns[0]}')
            plt.show()
            plt.close()
        
    def calc_ci(self):
        self.CI_list = list()
        for i in range(0, len(self.df_list) - 1):
            #slightly weird indexing  because self.df_list[0] is the threshold value
            #also examine frac value later, smaller values yield fits closer to the actual data points
            lowess = sm.nonparametric.lowess(self.df_list[i + 1]['nll'], self.df_list[i + 1][self.fit_params[i]], frac = .1)
            
            lowess_x = list(zip(*lowess))[0]
            lowess_y = list(zip(*lowess))[1]

            f = interp1d(lowess_x, lowess_y, bounds_error=False)
            f_thresh = lambda x: f(x) - self.df_list[0]

            #five makes it a slightly arbitrary guess but not terrible since there are currently 20 points being ran
            lb = newton(f_thresh, self.df_list[i + 1][self.fit_params[i]][5])
            ub = newton(f_thresh, self.df_list[i + 1][self.fit_params[i]][len(self.df_list[i + 1]) - 5])
            
            self.CI_list.append({'lb':lb, 'ub':ub})
        
    def save_fit_output(self, disease_name):
        #write a csv file with the parameter outputs
        fitted_params = pd.DataFrame(columns=self.fit_params, index = [1])
        for k in fitted_params.columns:
            fitted_params[k] = self.fit_out.params[k].value
        path_param_values = os.path.join(self.config_dict['PARAMETER_FIT_DIR'],
                                       f'{disease_name}_parameter_values.csv')
        fitted_params.to_csv(path_param_values, index=False)
        #write a text file with the full output
        path_param_out = os.path.join(self.config_dict['PARAMETER_FIT_DIR'],
                                       f'{disease_name}_parameter_output.txt')
        with open(path_param_out, 'w') as fh:
            fh.write(fit_report(self.fit_out))
        if self.calc_ci_bool == True:
            fitted_params_CI = pd.DataFrame(columns=self.fit_params, index = ['Lower CI', 'Upper CI'])
            for i in range(0, len(self.CI_list)):
                CI = self.CI_list[i]
                fitted_params_CI.loc['Lower CI', self.fit_params[i]] = CI['lb']
                fitted_params_CI.loc['Upper CI', self.fit_params[i]] = CI['ub']
            path_param_values_CI = os.path.join(self.config_dict['PARAMETER_FIT_DIR'],
                                       f'{disease_name}_parameter_values_CI.csv')
            fitted_params_CI.to_csv(path_param_values_CI)
    

    def error_check_resolution(self):
        """check if resolution in configuration file is 1 for a daily run
        
        check if all fitting data resolutions are listed as either daily or weekly
        """
        try:
            if self.config_dict['RESOLUTION'] != 1:
                raise ValueError('Model run resolution must be one')
        except ValueError as e:
            self.logger.exception('Model run resolution must be one')
            raise e
        
        try:
            if len([x for x in self.fit_data_res if x != 'daily' and x != 'weekly']):
                raise ValueError('Fitting data resolution must be either "daily" or "weekly"')
        except ValueError as e:
            self.logger.exception('Fitting data resolution must be either "daily" or "weekly"')
            raise e
            
    def error_check_compartment_names(self):
        """check if fitting compartment names match an initial state name
        
        check if fitting compartment names matches a column name in corresponding fitting data
        """
        try:
            #added the two states that we sometimes add after running the model
            if sum([i in (list(self.initial_states.keys()) + ['Ih', 'Dh']) for i in self.fit_data_compartment]) != len(self.fit_data_compartment):
                raise ValueError('Fitting data compartment names must match initial state names')
        except ValueError as e:
            self.logger.exception('Fitting data compartment names must match initial state names')
            raise e
        
        data_colnames = [i.columns for i in self.fit_data]
        try:
            if len([x for x in self.fit_data_compartment if x in data_colnames[self.fit_data_compartment.index(x)]]) != len(self.fit_data_compartment):
                raise ValueError('Fitting data compartment name must match a data column name')
        except ValueError as e:
            self.logger.exception('Fitting data compartment name must match a data column name')
            raise e
    
    def error_check_data_and_output_length(self):
        """check if the number of rows in fitting data sets will match the number of rows in the corresponding model output 
           check if the fitting data has any NAs
        """
        data_rownum = [i.shape[0] for i in self.fit_data]
        model_rownum = []
        na_list = []
        for i in range(0,len(self.fit_data)):
            if self.fit_data_res[i] == 'daily':
                model_rownum.append(len(self.t_eval))
            elif self.fit_data_res[i] == 'weekly':
                model_rownum.append(len(self.t_eval[::7]))
            na_list.append(self.fit_data[i].isnull().values.any())
        rows_no_match = [i for i, item in enumerate(data_rownum) if item != model_rownum[i]]
        
        try:
            if len(rows_no_match) != 0:
                raise ValueError(f'Number of fitting data {", ".join(str(x) for x in rows_no_match)} rows will not match number of model output rows')
        except ValueError as e:
            self.logger.exception(f'Number of fitting data {", ".join(str(x) for x in rows_no_match)} rows will not match number of model output rows')
            raise e
        
        try:
            if sum(na_list) > 0:
                raise ValueError('At least one fitting data set contains NaN')
        except ValueError as e:
            self.logger.exception('At least one fitting data set contains NaN')
            raise e
        

