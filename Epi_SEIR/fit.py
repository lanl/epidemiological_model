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
import numdifftools
from scipy.integrate import solve_ivp
from lmfit import Parameters, minimize, fit_report
from scipy.stats import poisson
from scipy.stats import norm
from scipy.stats import nbinom
from scipy.stats import chi2

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
    
    def _init_fit_parameters(self):
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
            return(params_obj)
    
    #adding this for the proflikelihood calculation
    def _init_proflik_parameters(self, lik_fitparams):
        #create the Parameter objects
            params_obj = Parameters()
            #add selected parameters and guesses into the framework
            param_keys = [i for i in lik_fitparams if i in list(self.params.keys())]
            init_keys = [i for i in lik_fitparams if i in list(self.initial_states.keys())]
            for k in param_keys:
                params_obj.add(k, lik_fitparams[k].value, min = lik_fitparams[k].value*.5, max = lik_fitparams[k].value*1.5)
                #params_obj.add(k, self.params[k])
            for j in init_keys:
                params_obj.add(j, lik_fitparams[j].value, min = lik_fitparams[j].value*.5, max = lik_fitparams[j].value*1.5)
                #params_obj.add(j, self.initial_states[j])
            return(params_obj)
         
       
    def fit_run_model(self, params_fit):
        keys = list(self.initial_states.keys())
        self.model_output = np.empty([0, len(keys)])
        
        #param_keys = [i for i in self.fit_params if i in list(self.params.keys())]
        #init_keys = [i for i in self.fit_params if i in list(self.initial_states.keys())]
        
        param_keys = [i for i in params_fit if i in list(self.params.keys())]
        init_keys = [i for i in params_fit if i in list(self.initial_states.keys())]
        
        for k in param_keys:
            self.params[k] = params_fit[k]
        for j in init_keys:
            self.initial_states[j] = params_fit[j]
            
        try:
            #Note we are not not inputting the parameters as arguments, but it is working without that
            sol = solve_ivp(self.model_func, self.t, list(self.initial_states.values()), t_eval=self.t_eval)
            out = sol.y.T
        except Exception as e:
            self.logger.exception('Exception occurred running model for fitting')
            raise e
        #Note: columns are abbrevations (ex. Rh) instead of full name (ex. Recovered Humans) here
        self.model_df = pd.DataFrame(dict(zip(list(self.initial_states.keys()), out.T)))
    
    def fit_objective(self, params_fit):
        if self.fit_method == 'res':
            resid = np.empty([0])
            self.fit_run_model(params_fit)
            for i in range(0,len(self.fit_data)):
                res = self.fit_data_res[i]
                compartment = self.fit_data_compartment[i]
                df = self.fit_data[i][compartment]
                if res == "weekly":
                    week_out = self.model_df.iloc[::7, :].reset_index()
                    #added below for getting weekly cases
                    week_out['Dh'] = week_out['Ch'].diff().fillna(0)
                    out = week_out[compartment]
                elif res == "daily":
                     #added below for getting daily ccases
                    self.model['Dh'] = self.model_df['Ch'].diff().fillna(0)
                    out = self.model_df[compartment]
                    #think sometime about weighting this based on the amount of data each source has
                resid = np.concatenate((resid,out- df))
                #already flat here because of concatenating the residuals
            return resid
        elif self.fit_method == 'pois':
            data = np.empty([0])
            mod_out = np.empty([0])
            self.fit_run_model(params_fit)
            for i in range(0,len(self.fit_data)):
                res = self.fit_data_res[i]
                compartment = self.fit_data_compartment[i]
                df = self.fit_data[i][compartment]
                if res == "weekly":
                    week_out = self.model_df.iloc[::7, :].reset_index()
                    #added below for getting weekly cases
                    week_out['Dh'] = week_out['Ch'].diff().fillna(0)
                    out = week_out[compartment]
                elif res == "daily":
                     #added below for getting daily ccases
                    self.model['Dh'] = self.model_df['Ch'].diff().fillna(0)
                    out = self.model_df[compartment]
                    #think sometime about weighting this based on the amount of data each source has
                data = np.concatenate((data,df))
                mod_out = np.concatenate((mod_out,out))
            #adding a fudge factor for the log currently, because getting all zeros due to terrible fit
            return -sum(np.log(poisson.pmf(np.round(data),np.round(mod_out)) + 0.00001))
        elif self.fit_method == 'nbinom':
            data = np.empty([0])
            mod_out = np.empty([0])
            self.fit_run_model(params_fit)
            for i in range(0,len(self.fit_data)):
                res = self.fit_data_res[i]
                compartment = self.fit_data_compartment[i]
                df = self.fit_data[i][compartment]
                if res == "weekly":
                    week_out = self.model_df.iloc[::7, :].reset_index()
                    #added below for getting weekly cases
                    week_out['Dh'] = week_out['Ch'].diff().fillna(0)
                    out = week_out[compartment]
                elif res == "daily":
                     #added below for getting daily ccases
                    self.model['Dh'] = self.model_df['Ch'].diff().fillna(0)
                    out = self.model_df[compartment]
                    #think sometime about weighting this based on the amount of data each source has
                data = np.concatenate((data,df))
                mod_out = np.concatenate((mod_out,out))
            sigma = 0.1*np.mean(data)  # example WLS assuming sigma = 0.1*mean(data)
            #https://stackoverflow.com/questions/62454956/parameterization-of-the-negative-binomial-in-scipy-via-mean-and-std
            #I believe below is correct, but will want to double check
            p = [k/sigma**2 for k in mod_out]
            n = [k*p/(1.0-p) for k in mod_out]
            #n is same as here https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html 
            p = [k/sigma**2 for k in mod_out]
            return -sum(np.log(nbinom.pmf(np.round(data),np.round(n), p)))
        elif self.fit_method == 'norm':
            data = np.empty([0])
            mod_out = np.empty([0])
            self.fit_run_model(params_fit)
            for i in range(0,len(self.fit_data)):
                res = self.fit_data_res[i]
                compartment = self.fit_data_compartment[i]
                df = self.fit_data[i][compartment]
                if res == "weekly":
                    week_out = self.model_df.iloc[::7, :].reset_index()
                    #added below for getting weekly cases
                    week_out['Dh'] = week_out['Ch'].diff().fillna(0)
                    out = week_out[compartment]
                elif res == "daily":
                    #added below for getting daily ccases
                    self.model['Dh'] = self.model_df['Ch'].diff().fillna(0)
                    out = self.model_df[compartment]
                    #think sometime about weighting this based on the amount of data each source has
                data = np.concatenate((data,df))
                mod_out = np.concatenate((mod_out,out))
            return -sum(np.log(norm.pdf(data,mod_out,0.1*np.mean(data)))) # example WLS assuming sigma = 0.1*mean(data)
        
    
    @timer 
    def fit_constants(self):
        #should return self.fit_out.success here as well 
        try:
            if self.fit_method == 'res':
                self.fit_out = minimize(self.fit_objective, self._init_fit_parameters())
            else:
                self.fit_out = minimize(self.fit_objective, self._init_fit_parameters(), method = 'nelder')
        except Exception as e:
            self.logger.exception('Exception occurred when running minimization for model fitting')
            raise e
    
    #not going to add the specifics of the numbers calculated within the range to the config file at the moment 
    #need to make sure below is within the min/max range listed in the config file
    def _calc_param_range(self, param, perc = .01, num = 10):
        param_seq = np.linspace((1-perc)*param, (1 + perc)*param, 2*num)
        return param_seq
        
    def proflike(self):
        threshold = self.fit_objective(self.fit_out.params) + chi2.ppf(.95, len(self.fit_params))/2
        df_list = list()
        df_list.append(threshold)
        perc_dict = dict(zip(self.fit_params, [.1,.1, .01]))
        for k in self.fit_params:
            #reset self.params to the fit values for the start of every run
            param_keys = [i for i in self.fit_params if i in list(self.params.keys())]
            init_keys = [i for i in self.fit_params if i in list(self.initial_states.keys())]
            for j in param_keys:
                self.params[j] = self.fit_out.params[j].value
            for h in init_keys:
                self.initial_states[h] = self.fit_out.params[h].value
            
            #only want to fit the other parameters
            lik_fitparams = self.fit_out.params.copy()
            del lik_fitparams[k] 
            
            #set up parameters to sequence through and lists for results
            
            if k in list(self.params.keys()):
                param_seq = self._calc_param_range(self.params[k], perc = perc_dict[k])
            elif k in list(self.initial_states.keys()):
                param_seq = self._calc_param_range(self.initial_states[k], perc = perc_dict[k])
                
            nll = list()
            fit_success = list()
            
            #calculate nll through the sequenced parameters
            for p in param_seq:
                if k in list(self.params.keys()):
                    self.params[k] = p
                elif k in list(self.initial_states.keys()):
                    self.initial_states[k] = p
                fit = minimize(self.fit_objective, self._init_proflik_parameters(lik_fitparams), method = 'nelder')
                #update guesses for the next round
                for i in lik_fitparams:
                    lik_fitparams[i].value = fit.params[i].value
                fit_success.append(fit.success)
                nll.append(self.fit_objective(fit.params))
            
            df = pd.DataFrame({k:param_seq, 'nll': nll, 'success': fit_success})
            df_list.append(df)
        return df_list
           
        #next steps: 1) find some easy way to store the results for each parameter
        #            2) add chi square threshhold and find some easy way to find intersection between fit (maybe NLL~loess(param_values))?
        #.           3) find an easy way to save CIs for each parameter
        #            4) Test this out with dengue data!        
                                              
        
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
        #commenting this out, because this will not always be true with how we need to move to daily/weekly cases
        #try:
           # if sum([i in list(self.initial_states.keys()) for i in self.fit_data_compartment]) != len(self.fit_data_compartment):
               # raise ValueError('Fitting data compartment names must match initial state names')
        #except ValueError as e:
            #self.logger.exception('Fitting data compartment names must match initial state names')
            #raise e
        
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
        

