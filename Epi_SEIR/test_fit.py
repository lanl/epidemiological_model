from utils import create_logger
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
#will need to add below to the condo environment
from lmfit import Parameters, minimize, fit_report
import vbdm

class FitSEIRModel(vbdm.VectorBorneDiseaseModel):
    def __init__(self, config_file, disease_name):
            self.logger = create_logger(__name__, config_file)

            super().__init__(config_file, disease_name)

            #not sure if this correct, check later (used in init_parameters and model_func_fit)
            #in config can just put parameters_fit as under dengue PARAMETERS_FIT: [ "gamma_h", "nu_h"]
            self.params_to_fit = [k: self.params[k] for k in self.config_dict[disease_name]['PARAMETERS_FIT']]
            #OLD: need to add this to config file
            #OLD: self.fit_data = pd.read_csv(self.config_dict['FIT_DATA_FILE_PATH'])
            #structure 'FIT_DATA_FILE_PATH' as a list [file_path1, file_path2, etc.]
            self.fit_data = []
            #length is correct, will only go 1 less than the length
            for i in range(0, len(self.config_dict['FIT_DATA_FILE_PATH'])):
                self.fit_data[i] = pd.read_csv(self.config_dict['FIT_DATA_FILE_PATH'][i])
            

@abstractmethod
    def model_func_fit(self, params_fit):
        pass

#some way to specify what compartments we are pulling, ideally outputs list of those comparments through the inital states
#would ideally return dictionary with values {'Im':self.initial_states['Im'], 'Ih':self.initial_states['Ih']}, could also be function of mult compartments
@abstractmethod
    def fit_compartments(self):
        pass

    #i think I want to do the first underscore in the function this way
    def _init_parameters(self):
        #create the Parameter objects
            params_obj = Parameters()
            #add selected parameters and guesses into the framework
            for k in self.params_to_fit:
                params_obj.add(k, self.params_to_fit[k])
            return(params_obj)
    
    def fit_model_out(self, params_fit):
        keys = list(self.initial_states.keys())
        initial_values = np.array(list(self.fit_compartments().values()))
        self.model_output = np.empty([0, len(initial_values)])
        #set initial conditions as first obs
        self.model_output = np.concatenate((self.model_output, initial_values))
        
        
        #t_eval can be the same as t since we will want to pull max 1 model output per day (will pull 1 evaluation for each day)
        t = (0,1)
        t_eval = (0,1)
        
        #would need to add 'FIT TIME STOP' to config
        #starting with one for indexing model output correctly
        for i in range(1, self.config_dict['FIT_TIME_STOP']):
            #changed indexing to i-1 because changed range start to 1
            self.initial_states['Sv'] = self.mosq[i-1]

            try:
                #will need params_fit as an argument in self.model_func
                sol = solve_ivp(self.model_func_fit, t, list(self.initial_states.values()), t_eval=t_eval, args = params_fit)
                out = sol.y.T
                last_out = out[-1]
            except Exception as e:
                self.logger.exception('Exception occured running model')
                raise e
    
            self.initial_states = dict(zip(keys, last_out))
            
            #OLD: since inital states are reset above, can still use fit_compartments function
            #OLD: week_out[i % self.config_dict['FIT DATA RESOULTION']] = list(self.fit_compartments().values())
            #below, if you have a data resolulation length model output you pull the last entry
            #OLD: if i % self.config_dict['FIT DATA RESOULTION'] == self.config_dict['FIT DATA RESOULTION'] - 1:
                #Currently taking the last row of each weeks/time resolution output
                #so probably don't necessarily have to actually store it every week, but I will for now since things may get more complicated later
                #OLD: self.model_output = np.concatenate((self.model_output, week_out[-1]))
            
            #since inital states are reset above, can still use fit_compartments function (note: appending values needs to be in double brackets)
            self.model_output = np.concatenate((self.model_output, [list(self.fit_compartments().values())]), axis = 0)
    
    def fit_objective(self, params_fit):
        #OLD: resid = np.empty((len(self.fit_data.columns), len(self.fit_data.index)))
        #OLD: calculating residual between each data column and corresponding model output
        #OLD: fit_keys = list(self.fit_compartments().keys())
        #OLD: self.fitting_output = pd.DataFrame(self.fit_model_out(params_fit))
        #OLD: self.fitting_output.columns = fit_keys
        #OLD: for k in fit_keys:
            #OLD: resid[fit_keys.index(k)] = self.fit_data[k] - self.fitting_output[k]
        # now flatten this to a 1D array, as minimize() needs
        #return resid.flatten()
        resid = np.empty([0])
        #prep names, and read in model output
        fit_keys = list(self.fit_compartments().keys())
        fitting_output = pd.DataFrame(self.fit_model_out(params_fit))
        fitting_output.columns = fit_keys
        #create time column for joining: NOTE then each fit data needs to have a column that counts from 0 to end date
        fitting_output['time'] = np.arange(0, self.config_dict['FIT_TIME_STOP'])
        for k in fit_keys:
            #make this config dict a dictionary where the keys are the fitting compartments 
            df = self.config_dict['FIT DATA RESOULTION'][k]
            out = fitting_output[['time', k]]
            #pull correct times by selecting times
            out_select = out[out['time'].isin(df['time'])
            #calculate the residual
            resid = np.concatenate((resid,out_select[k] - df[k]))
        #already flat here
        return resid
            
        
        
    
    initial_params = self._init_parameters()
    def fit_parameters(self, params_fit = initial_params):
        self.fit_out = minimize(self.fit_objective, params_fit)
        
    def save_fit_output(self, disease_name)
    #Need to create fit_output folder
        #write a csv file with the parameter outptus
        fitted_params = pd.DataFrame(columns=self.params_to_fit.keys(), index = [1])
        for k in fitted_params.columns:
            fitted_params[k] = self.fit_out.params[k].value
        path_param_values = os.path.join(self.config_dict['PARAMETER_FIT_DIR'],
                                       f'{disease_name}_parameter_values.csv')
        fitted_params.to_csv(output_path, index=False)
        #write a text file with the full output
        path_param_out = os.path.join(self.config_dict['PARAMETER_FIT_DIR'],
                                       f'{disease_name}_parameter_output.txt')
        with open(path_param_out, 'w') as fh:
            fh.write(fit_report(self.fit_out))
    
    #add into dengue.py
    #going to need to add another super_init thing, I want the fitting file to be in between the two others
    def model_func_fit(self, t, y, params_fit):
        for k in in self.params_to_fit:
            self.params[k] = params_fit[k]
        return self.model_func(t,y)

    #this is for dengue
    def fit_compartments(self):
        return {'Ch':self.initial_states['Ch']}  
                             
                             
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser
import sys

def fit():
    """Fits vector born disease model to data.

    Instaniates each disease from configuration file
    """
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser()
        args, unknown = parser.parse_known_args()

    disease_name = args.disease_name.lower()

    if disease_name == 'dengue':
        disease = DengueSEIRModel(args.config_file)
    elif disease_name == 'wnv':
        disease = WNVSEIRModel(args.config_file)

    disease.logger.info(disease)
    disease.fit_parameters()
    disease.save_fit_output(disease_name)
    disease.logger.info('SUCCESS')

#ask jeff about this: he says it is helpful for unit testing (more details in mattermost)
if __name__ == "__main__":
    fit()