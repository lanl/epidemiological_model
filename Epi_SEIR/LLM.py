"""
This script fits time-varying logistic growth parameters to the mosquito PBM output

"data" is obtained from the Mosquito PBM output

@author: Marina Mancuso (mmancuso@lanl.gov)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.integrate import odeint
from scipy.optimize import least_squares
import os
import yaml
from multiprocessing import Pool
from datetime import datetime
from utils import create_logger
from abc import ABC, abstractmethod

class LogisticLinkModel(ABC):
    """Defines the general logistic link model.

    Reads root config file (absolute path) and preps data to fit to.

    Attributes:
        config_dict: All configurations read in from config file.\n

    """
    def __init__(self, config_file, mosq_data):
        self.logger = create_logger(__name__, config_file)
        
        self._read_config(config_file)
        #self.error_fit_data_file_path()
        
        self.dates = self.config_dict['DATES']
        #self.location = self.config_dict['POPULATION']['location']
        self.pop_type = self.config_dict['POPULATION']['type']
        self.data_col_names = self.config_dict['DATA_COL_NAMES'] 
        #new code
        #self.mosq_data_path = self.config_dict['MOSQUITOES_DIR']
        
        #prep data
        #self.full_df = pd.read_csv(self.config_dict['DATA_PATH'])
        #self.full_df = pd.read_csv(mosq_data)
        self.df = pd.read_csv(mosq_data)
        
        #self.error_fit_data_columns()
        #self.error_pop_type()
        #self.error_date_type()
        
        #self.df = self.full_df[self.full_df[self.data_col_names['location']] == self.location].sort_values(by=self.data_col_names['date'])
        
        #self.error_data_date_nan()
        
        self.ts_length = len(self.df[self.data_col_names['date']])
        self.ts_start = pd.to_datetime(self.df[self.data_col_names['date']]).min()
        self.ts_range = pd.date_range(start=self.ts_start, end=pd.to_datetime(self.df[self.data_col_names['date']]).max())
        
        if self.pop_type == 'Total':
            self.PBM_type = self.data_col_names['total_mosq']
        elif self.pop_type == 'Active':
            self.PBM_type = self.data_col_names['active_mosq']
        
        date_mosq = self.df[[self.data_col_names['date'], self.PBM_type]]
        date_mosq[self.data_col_names['date']] = [datetime.strptime(k, '%Y-%m-%d') for k in date_mosq[self.data_col_names['date']]]
        #do this dataset construction and merge to make sure the date numbers all fit together even if the mosquito data has some missing dates
        date_num_trans = pd.DataFrame({'day_index': np.arange(0, len(self.ts_range)), self.data_col_names['date']:self.ts_range})
        self.date_mosq_final = date_mosq.merge(date_num_trans, how = 'left', left_on = self.data_col_names['date'], right_on = self.data_col_names['date'])
        
        #self.error_fit_data()
        
    
    def _prep_dates(self):
        """Create an array where columns are start date, end date, year and rows are all possible combinations of those three by inputs"""
        #Get array of all possible start dates
        start_range = pd.date_range(start=datetime.strptime(self.dates['start1'], '%m-%d'), end=datetime.strptime(self.dates['start2'], '%m-%d'))
        start_range_format = np.array([k.strftime('%m-%d') for k in start_range])

        years_format = np.arange(int(self.dates['start_year']), int(self.dates['end_year']) + 1, 1)
        
        if self.pop_type == 'Total':
                    #Get array of all possible end dates
            end_range = pd.date_range(start=datetime.strptime(self.dates['end1'], '%m-%d'), end=datetime.strptime(self.dates['end2'], '%m-%d'))
            end_range_format = np.array([k.strftime('%m-%d') for k in end_range])

            #create dataset of all possible start, end, and year combinations where each row contains a start date, end date, and year in that order
            all_vals = np.array(np.meshgrid(*[start_range_format, end_range_format, years_format])).reshape(3, len(start_range_format) * len(end_range_format) * len(years_format)).T
        elif self.pop_type == 'Active':
            end_dates = start_range - pd.offsets.Day()
            end_dates_format = np.array([k.strftime('%m-%d') for k in end_dates])
            
            start_end = pd.DataFrame({'start':start_range_format, 'end':end_dates_format})
            start_year = pd.DataFrame(np.array(np.meshgrid(*[start_range_format, years_format])).reshape(2, len(start_range_format) * len(years_format)).T, columns = ['start','year'])
            start_year['year'] = start_year['year'].apply(int)
            all_vals = start_year.merge(start_end, how = 'left')[['start', 'end', 'year']].values
            
        
        return all_vals
    
        
    def _compute_r_and_K(self, rb, rs, Kb, Ks, t):
        """Computes growth rate (r) and carrying capacity (K) values"""
        r = rb - rs * np.cos(2 * np.pi * t / 365)
        K = Kb - Ks * np.cos(2 * np.pi * t / 365)
        
        return (r,K)
    
    
    def _logistic_model(self, r, K, S):
        """The differential logistic equation for the mosquito population"""
        dS = r * S * (1 - S / K)
        
        return dS
    
    
    def _RK4(self, x0, delta_t, params, t):
        """4th order Runge Kutta solver for non-autonomous logistic model"""
        # total number of steps = 1/(step size)
        n = int((1 / delta_t) * len(t))
        # create a vector for S
        S = np.zeros(int(n))
        # extract parameter components from input    
        rb, rs, Kb, Ks = params
        # assign initial condition to first element of S
        S[0] = x0

        for i in np.arange(0,n-1):
            # compute r and K for time step t_i
            t_i = i * delta_t
            r1, K1 = self._compute_r_and_K(rb, rs, Kb, Ks, t_i)
            k1 = delta_t * self._logistic_model(r1, K1, S[i])
            # compute r and K for time step t_i + (delta_t)/2
            r23, K23 = self._compute_r_and_K(rb, rs, Kb, Ks, t_i + (delta_t)/2)
            k2 = delta_t * self._logistic_model(r23, K23, S[i] + k1/2)
            k3 = delta_t * self._logistic_model(r23, K23, S[i] + k2/2)
            # compute r and K for time step t_i + (delta_t)/2
            r4, K4 = self._compute_r_and_K(rb, rs, Kb, Ks, t_i + (delta_t))
            k4 = delta_t * self._logistic_model(r4, K4, S[i] + k3)
            # new computed numerical value
            comp_value = S[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            # if new computed value exceeds the max possible
            #   carrying capacity, set it equal to the previous value
            #   This avoids the risk of exceeding carrying capacity
            if comp_value <= Kb + Ks:
                S[i + 1] = comp_value
            else:
                   S[i + 1] = S[i] 
                    
        return S
    
    
    def _offset_val(self, ts_range):
        """The nature of the non-autonomous parameter functions make fitting sensitive to the start of the mosquito PBM time series 
           This function calculates the offset value: the difference between the initial PBM start date and the location's typical season beginning"""
       # if self.location == 'Toronto':
        #will probably want to add the offset 
        offset_start = pd.to_datetime(f'May 1, {self.ts_start.year}')
        offset = (self.ts_start - offset_start).days

        return offset
    
    
    def _cost(self, init_params, mosq, mosq_dates, ts_range):
        """Cost function for parameter optimization"""
        init_params = tuple(init_params)
        mosq_start, mosq_end = mosq_dates

        # initial condition
        #x0 = list(mosq)[0]
        x0 = np.max([0.01, list(mosq)[0]])

        # step size
        delta_t = 1

        # time frame for fitting
        offset = self._offset_val(ts_range)
            # NOTE: t(0) starts on offset day, not the first day of mosquito PBM time series
        t = np.arange(mosq_start, mosq_end) + offset

        # logistic link model simulation
        model = self._RK4(x0, delta_t, init_params, t)

        # error between mosquito PBM time series and model
        res = model[::int(1/delta_t)] - mosq

        return res.flatten()
    
    
    def _fit_data(self, init_params, mosq, mosq_start, mosq_end, ts_range):
        offset = self._offset_val(ts_range)
        mosq_dates = [mosq_start, mosq_end]

        # bound constraints for parameters
        if self.pop_type == 'Total':
            # bd = ([-0.1, -0.1, 70000, -2e5],[0.1, 0, 6e5, 6e5]) # baseline
            # bd = ([-0.2, -0.3, 1, 0],[0.2, 0, 6e5, 6e5]) # baseline
            bd = ([-0.2, -0.4, 1, -1e5],[0.2, 0, 1e5, 1e5]) # baseline
            # bd = ([-0.11, -0.11, 70000, -2e5],[0.11, 0, 6e5, 6e5]) # r10
            # bd = ([-0.125, -0.125, 70000, -2e5],[0.125, 0, 6e5, 6e5]) # r25
            # bd = ([-0.15, -0.15, 70000, -2e5],[0.15, 0, 6e5, 6e5]) # r50
            # bd = ([-0.2, -0.2, 70000, -2e5],[0.2, 0, 6e5, 6e5]) # r100
            # bd = ([-0.1, -0.1, 63000, -2.2e5],[0.1, 0, 6.6e5, 6.6e5]) # K10
            # bd = ([-0.1, -0.1, 52500, -2.5e5],[0.1, 0, 7.5e5, 7.5e5]) # K25
            # bd = ([-0.1, -0.1, 35000, -3e5],[0.1, 0, 9e5, 9e5]) # K50
            # bd = ([-0.1, -0.1, 0, -4e5],[0.1, 0, 10.2e5, 10.2e5]) # K100
            # bd = ([-0.11, -0.11, 63000, -2.2e5],[0.11, 0, 6.6e5, 6.6e5]) # rK10
            # bd = ([-0.125, -0.125, 52500, -2.5e5],[0.125, 0, 7.5e5, 7.5e5]) # rK25
            # bd = ([-0.15, -0.15, 35000, -3e5],[0.15, 0, 9e5, 9e5]) # rK50
            # bd = ([-0.2, -0.2, 0, -4e5],[0.2, 0, 10.2e5, 10.2e5]) # rK100
        elif self.pop_type == 'Active':
            # bd = ([-0.1, -0.1, 1, -2e5],[0.1, 0, 6e5, 6e5]) # attempt 1 (original)
            # bd = ([-0.1, -0.2, 1, -1e4],[0.1, 0, 1e4, 1e4]) # attempt 2
            # bd = ([-0.1, -0.2, 1, 0],[0.1, 0, 1e4, 1e4]) # attempt 3
            # bd = ([-0.1, -0.3, 1, 0],[0.1, 0, 5e4, 1e4]) # attempt 4
            # bd = ([-0.1, -0.25, 1, 0],[0.1, 0, 5e4, 1e4]) # attempt 5
            # bd = ([-0.1, -0.3, 1, 0],[0.1, 0, 1e4, 1e4]) # attempt 6
            bd = ([-0.1, -0.25, 1, 0],[0.1, 0, 1e4, 1e4]) # attempt 7

        
#         # optimization run
#         fitted_params = least_squares(cost, init_params, bounds = bd, method = 'trf', args=(mosq, mosq_dates, ts_range)).x

#         # simulate model with optimized parameters
#         x0 = list(mosq)[0]
#         t = np.arange(mosq_start, mosq_end) + offset
#         delta_t = 1
#         fit_x = RK4(x0, delta_t, fitted_params, t)

        optim = least_squares(self._cost, init_params, bounds = bd, method = 'trf', args=(mosq, mosq_dates, ts_range))
        fitted_params = optim.x
        resid = optim.fun

        #return (fit_x[::int(1/delta_t)], fitted_params, t)
        return (resid, fitted_params)
    
    def _all_func(self, times):
        #print(times)
        start = times[0]
        end = times[1]
        year = str(times[2])
        
        if self.pop_type == 'Active':
            year_next = str(times[2] + 1)

        start_full = datetime.strptime(year + '-' + start, '%Y-%m-%d')
        
        if self.pop_type == 'Total':
            end_full = datetime.strptime(year + '-' + end, '%Y-%m-%d')
        elif self.pop_type == 'Active':
            end_full = datetime.strptime(year_next + '-' + end, '%Y-%m-%d')

        mosq_start = self.date_mosq_final[self.date_mosq_final[self.data_col_names['date']] == start_full]['day_index'].reset_index()['day_index'][0]
        mosq_end = self.date_mosq_final[self.date_mosq_final[self.data_col_names['date']] == end_full]['day_index'].reset_index()['day_index'][0]

        # find subset of mosq PBM output
        #end is not less than or equal too - confirm with marina but this matches the previous indexing
        mosq = np.array(self.date_mosq_final[(self.date_mosq_final['date'] >= start_full) & (self.date_mosq_final['date'] < end_full)][self.PBM_type])

        # initialize parameter values
        r_b_init = 0
        r_s_init = -0.07

        if self.pop_type == 'Total':
            K_b_init = 8e4
            K_s_init = 100
        elif self.pop_type == 'Active':
            K_b_init = np.mean(mosq)
            K_s_init = 1

        init_params = np.array([r_b_init, r_s_init, K_b_init, K_s_init])

        # run data fitting
        mosq_fit = self._fit_data(init_params, mosq, mosq_start, mosq_end, self.ts_range)

        #pull fit params
        rb, rs, Kb, Ks = mosq_fit[1]

#         norm_vec = ts_sub - mosq_fit[0].T
#         days_fit = len(mosq_fit[0])    # normalize w.r.t amount of days fit

#         # mean squared error
#         MSE = np.linalg.norm(norm_vec)/days_fit
#         RMSE = (MSE)**(1/2)
        RMSE = np.sqrt(np.mean(mosq_fit[0]**2))
        return [RMSE, rb, rs, Kb, Ks, times[0], times[1], times[2]]
    
    def run_model(self):
        self.all_date_vals = self._prep_dates()
        list_date_vals = self.all_date_vals.tolist()
        #fit_out = np.apply_along_axis(self._all_func, 1, self.all_date_vals)
        
        with Pool() as pool:
            res = pool.map(self._all_func, list_date_vals)
            #print(res)
            out = list(res)
    
        fit_out = np.array(out).reshape((len(out), 8))
        
        df_full = pd.DataFrame(fit_out, columns = ['RMSE', 'r', 'r_s', 'K', 'K_s', 'start_date', 'end_date', 'year'])

        
        if self.pop_type == 'Total':
            df_agg = df_full.groupby('year').agg({'RMSE':'min'}).reset_index()
            df_best = df_full[df_full['RMSE'].isin(df_agg['RMSE'])].sort_values('year').reset_index(drop = True)
        elif self.pop_type == 'Active':
            df_agg = df_full.groupby('start_date').agg({'RMSE':'sum'}).reset_index()
            df_min = df_agg[df_agg['RMSE'] == df_agg['RMSE'].min()].reset_index()
            df_best = df_full[df_full['start_date'] == df_min['start_date'][0]].sort_values('year').reset_index(drop = True)
        
        #join up to get the mosquito population data for each start date with the output data
        df_best['full_start_date'] = [datetime.strptime(str(df_best['year'][i]) + '-' + df_best['start_date'][i], '%Y-%m-%d') for i in range(0, len(df_best))]
        df_join = df_best.merge(self.date_mosq_final, how = 'left', left_on = 'full_start_date', right_on ='date')
        df_join = df_join.drop(['full_start_date', self.data_col_names['date'], 'day_index'], axis = 1).rename(columns = {self.PBM_type: f'{self.pop_type}_mosq_start'})
        #join up to get the mosquito population data for May 1st each year with the output data
        df_join['May_First'] = [datetime.strptime(str(k) + '-05-01', '%Y-%m-%d') for k in df_join['year']]
        df_join_final = df_join.merge(self.date_mosq_final, how = 'left', left_on = 'May_First', right_on ='date')
        self.full_out_df = df_join_final.drop(['May_First', self.data_col_names['date'], 'day_index'], axis = 1).rename(columns = {self.PBM_type: f'{self.pop_type}_mosq_May1'})
        #self.integration_df = self.full_out_df[['r', 'r_s', 'K', 'K_s', f'{self.pop_type}_mosq_May1']].rename(columns = {f'{self.pop_type}_mosq_May1': 'Sv'})
        #switching to the chosen start date - ask Marina about this
        self.integration_df = self.full_out_df[['r', 'r_s', 'K', 'K_s', f'{self.pop_type}_mosq_start']].rename(columns = {f'{self.pop_type}_mosq_start': 'Sv'})

    
    def save_output(self, name_tag):
        output_path = os.path.join(self.config_dict['LLM_OUTPUT_DIR'], f'full_output/{self.pop_type}_output_{name_tag}.csv')
        output_path2 = os.path.join(self.config_dict['LLM_OUTPUT_DIR'], f'{self.pop_type}_integration_output_{name_tag}.csv')
        
        self.full_out_df.to_csv(output_path, index=False)
        self.integration_df.to_csv(output_path2, index=False)
        
#         output_path = os.path.join(self.config_dict['OUTPUT_DIR'],
#                                            #f'{self.location}_{self.pop_type}_model_output_2005_2008.csv')
#                                             f'{self.location}_{self.pop_type}_output.csv')
#         output_path2 = os.path.join(self.config_dict['OUTPUT_DIR'],
#                                            #f'{self.location}_{self.pop_type}_model_output_2005_2008.csv')
#                                             f'{self.location}_{self.pop_type}_integration_output.csv')

        
    def _read_config(self, config_file):
        """Reads root configuration file"""
        try:
            with open(config_file, 'r') as in_file:
                self.config_dict = yaml.safe_load(in_file)
        except OSError as e:
            self.logger.exception('Exception occurred opening configuration file')
            raise e
        else:
            self.logger.info('Configuration file successfully opened')
            
    
#     def error_fit_data_file_path(self):
#         """check if fit data file is a valid file
        
#         check if fit data file is a csv
#         """
#         try:
#             if not os.path.isfile(self.config_dict['DATA_PATH']):
#                 raise ValueError('Fit data file must be a valid file')
#         except ValueError as e:
#             self.logger.exception('Fit data file must be a valid file')
#             raise e
    
#         try:
#             if 'csv' not in self.config_dict['DATA_PATH']:
#                 raise ValueError('Fit data file must be a csv')
#         except ValueError as e:
#             self.logger.exception('Fit data file must be a csv')
#             raise e
            
    
#     def error_fit_data_columns(self):
#         """check if config listed data columns actually match columns in the fit data
        
#         check if mosquito population location is in the location column of the fit data
#         """
#         try:
#             if len([i for i in list(self.data_col_names.values()) if i in self.full_df.columns]) != len(list(self.data_col_names.values())):
#                 raise ValueError('Data column names listed in config file do not match actual column names in the fitting data')
#         except ValueError as e:
#             self.logger.exception('Data column names listed in config file do not match actual column names in the fitting data')
#             raise e
        
#         try:
#             if self.location not in np.unique(self.full_df[self.data_col_names['location']]):
#                 raise ValueError('Population location not found in location column of fitting data')
#         except ValueError as e:
#             self.logger.exception('Population location not found in location column of fitting data')
#             raise e
            
            
#     def error_pop_type(self):
#         """check if mosquito population is listed as 'Total' or 'Active'
#         """
#         try:
#             if self.pop_type != "Total" and self.pop_type != "Active":
#                 raise ValueError('Mosquito population type must be "Total" or "Active"')
#         except ValueError as e:
#             self.logger.exception('Mosquito population type must be "Total" or "Active"')
#             raise e
    
    
#     def error_date_type(self):
#         """check if dates are inputted as strings
        
#         check if dates are correctly formatted as monthd day
#         """
#         try:
#             if len([i for i in list(self.dates.values()) if isinstance(i, str)]) != len(list(self.dates.values())):
#                 raise ValueError('Start and end dates must be string')
#         except ValueError as e:
#             self.logger.exception('Start and end dates must be string')
#             raise e
        
#         #fix below with current format of start, end, year config inputs
# #         date_options = ['start1', 'start2', 'end1', 'end2']
# #         for i in date_options:
# #             k = self.dates[i]
# #             try:
# #                 datetime.strptime(k, '%b %d')
# #             except ValueError as e:
# #                 self.logger.exception('Start and end dates must be in format "First three letters of month name day number" (ex. "Mar 1").')
# #                 raise e
    
#     def error_data_date_nan(self):
#         """check if there are missing mosquito population values
        
#         check if there are missing date values
#         """ 
#         try:
#             if self.df[self.data_col_names['date']].isnull().values.any() == True:
#                 raise ValueError('Date missing')
#         except ValueError as e:
#             self.logger.exception('Date missing')
#             raise e
            
    
#     def error_fit_data(self):
#         """check if there are missing mosquito population values
        
#         check if length of time series matches date ranges
#         """
#         try:
#             if np.isnan(np.sum(self.date_mosq_final[self.PBM_type])) == True:
#                 raise ValueError('Population value missing')
#         except ValueError as e:
#             self.logger.exception('Population value missing')
#             raise e
#         try:
#             if self.ts_length != len(self.ts_range):
#                 raise ValueError('Length of time series does not match range of dates')
#         except ValueError as e:
#             self.logger.exception('Length of time series does not match range of dates')
#             raise e
            
    

        
    

