"""Fitting Model Class.

Contains class for specific diseases that inherit from the
class FitModel. Inherits from the VectorBorneDiseaseModel class.

"""
import os
from utils import timer
from abc import abstractmethod

import numpy as np
import yaml
import pandas as pd
from scipy.integrate import solve_ivp
from lmfit import Parameters, minimize, fit_report

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
                params_obj.add(k, self.params[k], min = self.fit_params_range[k]['min'], max = self.fit_params_range[k]['max'])
            for j in init_keys:
                params_obj.add(j, self.initial_states[j], min = self.fit_params_range[j]['min'], max = self.fit_params_range[j]['max'])
            return(params_obj)
        
    
       
    def fit_run_model(self, params_fit):
        keys = list(self.initial_states.keys())
        self.model_output = np.empty([0, len(keys)])
        
        param_keys = [i for i in self.fit_params if i in list(self.params.keys())]
        init_keys = [i for i in self.fit_params if i in list(self.initial_states.keys())]
        
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
        resid = np.empty([0])
        self.fit_run_model(params_fit)
        for i in range(0,len(self.fit_data)):
            res = self.fit_data_res[i]
            compartment = self.fit_data_compartment[i]
            df = self.fit_data[i][compartment]
            if res == "weekly":
                week_out = self.model_df.iloc[::7, :].reset_index()
                out = week_out[compartment]
            elif res == "daily":
                out = self.model_df[compartment]
            #think sometime about weighting this based on the amount of data each source has
            resid = np.concatenate((resid,out- df))
        #already flat here because of concatenating the residuals
        return resid
    
    @timer 
    def fit_constants(self):
        try:
            self.fit_out = minimize(self.fit_objective, self._init_fit_parameters())
        except Exception as e:
            self.logger.exception('Exception occurred when running minimization for model fitting')
            raise e
            
        
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
        try:
            if sum([i in list(self.initial_states.keys()) for i in self.fit_data_compartment]) != len(self.fit_data_compartment):
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
        

