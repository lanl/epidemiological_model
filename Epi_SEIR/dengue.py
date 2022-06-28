"""Dengue SEIR Model Class

Contains class for dengue disease model. Inherits from
FitModel and Vector Borne Disease Model class.

    Typical usage example:

    den = DengueSEIRModel(<config_file_path>, <command_line_arguments>)
    
Class method param_dict allows you to change parameters values from configuration file inputs 
by inputting a dictionary of different parameter values
    
    Typical usage example:
    
    den_param_dict = DengueSEIRModel.param_dict(<config_file_path>, <parameter dictionary>)
"""

from utils import create_logger
import sys
import math
import fit
import pickle
from scipy.interpolate import splev

class DengueSEIRModel(fit.FitModel):

    """Models the spread of dengue.

    Inherits from the FitModel class. Specifies ODE system
    of equations.

    """

    def __init__(self, config_file, param_dict = None):
        self.logger = create_logger(__name__, config_file)
        if __name__ == "dengue_test":
            self.logger.disabled = True
        
        super().__init__(config_file, 'DENGUE')

        # Read in time dependent parameter splines
        try:
            #with open(self.config_dict['DENGUE']['PARAMETERS']['biting_rate'], "rb") as file:
            with open('parameters/test_spline.pkl', "rb") as file:
                self.biting_rate = pickle.load(file)
        except FileNotFoundError as e:
            self.logger.exception('Biting rate parameter spline file not found.')
            raise e
        else:
            self.logger.info('Biting rate parameter spline file successfully opened.')
        
        self.error_zero_constants()
        self.error_zero_to_one_params()
    
    @classmethod
    def param_dict(cls, config_file, param_dict):
        """Inherit vbdm param_dict class method"""
        return super(DengueSEIRModel, cls).param_dict(config_file = config_file, disease_name = 'DENGUE', param_dict =  param_dict)
    
    def _population_sizes(self):
        """Calculates population sizes of human and vector compartments"""
        self.Nh = sum([self.states['Sh'], self.states['Eh'], self.states['Ih'], self.states['Rh']])
        self.Nv = sum([self.states['Sv'], self.states['Ev'], self.states['Iv']])
    
    def _biting_rate(self):
        """Calculates biting rate DEPRECATED"""
        # evaluate a spline in main function. make spline a class attribute read in from a
        # .pkl

        # OLD TWO LINES
        #self.b = self.params['sigma_h'] * self.params['sigma_v'] / \
        #    (self.params['sigma_h'] * self.Nh + self.params['sigma_v'] * self.Nv)
    
    def _force_of_infection(self):
        """Calculates force of infection"""
        self.lambda_h = self.b * self.params['beta_h'] * self.states['Iv'] 
        self.lambda_v = self.b * self.params['beta_v'] * self.states['Ih']
        
    #def _birth_rate(self,t):
        #"""Caclualtes vector natural birth rate"""
        #if t < 180:
            #self.psi_v = self.params['r_v'] + self.params['mu_v']
        #else:
            #self.psi_v = 0
            #self.r_v = self.psi_v - self.params['mu_v']
    
    
    def _mosq_population_values(self, t):
        self.K_v = self.params['K'] - self.params['K_s'] * math.cos((2 * math.pi / 365))
        self.r_v = self.params['r'] - self.params['r_s'] * math.cos((2 * math.pi / 365))
    
    def _initial_human_pop(self):
        self.H0 = sum([self.initial_states['Sh'], self.initial_states['Eh'], self.initial_states['Ih'], self.initial_states['Rh']])
    
    def _birth_rate(self):
        """Caclualtes vector natural birth rate"""
        self.psi_v = self.r_v + self.params['mu_v']
    
    def _calc_hNv(self):
        self.hNv = self.psi_v - self.r_v * self.Nv / self.K_v

    def model_func(self, t, y):
        """Defines system of ODEs for dengue model.

        Initial State Names:
            Sh: Susceptible human population.\n
            Eh: Exposed human population.\n
            Ih: Infectious human population.\n
            Rh: Recovered human population.\n
            Nh: Human population.\n
            Sv: Susceptible vector population.\n
            Ev: Exposed vector population.\n
            Iv: Infectious vector population.\n
            Nv: Vector population.\n

        Parameters:
            H0: Initial human population.\n
            psi_h: Per capita birth rate for humans.\n
            psi_v: Per capita birth rate for mosquitoes.\n
            beta_h: Probability of vector to host transmission.\n
            beta_v: Probability of host to vector transmission.\n
            sigma_h: Maximum number of times one mosquito bite a human per day.\n
            sigma_v: Maximum number of mosquito bites a human can sustain per day.\n
            lambda_h: Human force of infection.\n
            lambda_v: Vector force of infection.\n
            alpha_h: Rate of lost dengue immunity.\n
            nu_h: Human latent period.\n
            nu_v: Vector latent period.\n
            gamma_h: Human infectious period.\n
            mu_v: Vector natural death rate.\n
            mu_h: Human natural death rate.\n
            r: Baseline vector growth rate factor.\n
            r_s: Scaling vector growth rate factor.\n
            r_v: Time-varying vector growth rate.\n
            K: Baseline vector carrying capacity factor.\n
            K_s: Scaling vector carrying capacity factor.\n
            K_v: Time-varying vector carrying capacity.\n

        """
        ddt = self.initial_states.copy()
        self.states = dict(zip(self.initial_states.keys(), y))

        # Find population size
        self._population_sizes()

        # Find biting rate
        #self._biting_rate()
        self.b = splev([t], self.biting_rate)[0]

        # Find force of infection
        self._force_of_infection()

        # Find mosquito carrying capcity and growth rate from LLM fit
        self._mosq_population_values(t)
        
        # Find initial human population
        self._initial_human_pop()
        
        # Find vector natural birth rate
        self._birth_rate()
        
        # Find hNv value
        self._calc_hNv()
        
        # System of equations
        ddt['Sh'] = self.params['psi_h'] * self.H0 - \
            self.lambda_h * self.states['Sh'] + \
            self.params['alpha_h'] * self.states['Rh'] - \
            self.params['mu_h'] * self.states['Sh']
        ddt['Eh'] = self.lambda_h * self.states['Sh'] - \
            self.params['nu_h'] * self.states['Eh'] - \
            self.params['mu_h'] * self.states['Eh']
        ddt['Ih'] = self.params['nu_h'] * self.states['Eh'] - \
            self.params['gamma_h'] * self.states['Ih'] - \
            self.params['mu_h'] * self.states['Ih']
        ddt['Rh'] = self.params['gamma_h'] * self.states['Ih'] - \
            self.params['alpha_h'] * self.states['Rh'] - \
            self.params['mu_h'] * self.states['Rh']
        ddt['Ch'] = self.params['nu_h'] * self.states['Eh']
        ddt['Sv'] = self.hNv * self.Nv - \
            self.lambda_v * self.states['Sv'] - \
            self.params['mu_v'] * self.states['Sv']
        ddt['Ev'] = self.lambda_v * self.states['Sv'] - \
            self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Ev']
        ddt['Iv'] = self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Iv']

        return tuple(ddt.values())
    
    
    def error_zero_constants(self):
        """check if human population is non-zero.

        check if vector carrying capcity is non-zero.
        """

        # EXTRACT list of position
        try:
            if sum([self.initial_states['Sh'], self.initial_states['Eh'], self.initial_states['Ih'], self.initial_states['Rh']]) == 0:
                raise ValueError('Human population size cannot be zero')
        except ValueError as e:
            self.logger.exception('Human population size cannot be zero')
            raise e
    
    def error_zero_to_one_params(self):
        """check if parameters that should be in [0,1] are
        """
        try:
            if self.params['nu_v'] < 0 or self.params['nu_v'] > 1:
                raise ValueError('nu_v must be in [0,1]')
        except ValueError as e:
            self.logger.exception('nu_v must be in [0,1]')
            raise e
            
        try:
            if self.params['nu_h'] < 0 or self.params['nu_h'] > 1:
                raise ValueError('nu_h must be in [0,1]')
        except ValueError as e:
            self.logger.exception('nu_h must be in [0,1]')
            raise e
        
        try:
            if self.params['mu_v'] < 0 or self.params['mu_v'] > 1:
                raise ValueError('mu_v must be in [0,1]')
        except ValueError as e:
            self.logger.exception('mu_v must be in [0,1]')
            raise e
        
        try:
            if self.params['gamma_h'] < 0 or self.params['gamma_h'] > 1:
                raise ValueError('gamma_h must be in [0,1]')
        except ValueError as e:
            self.logger.exception('gamma_h must be in [0,1]')
            raise e
        
        try:
            if self.params['beta_h'] < 0 or self.params['beta_h'] > 1:
                raise ValueError('beta_h must be in [0,1]')
        except ValueError as e:
            self.logger.exception('beta_h must be in [0,1]')
            raise e
            
        try:
            if self.params['beta_v'] < 0 or self.params['beta_v'] > 1:
                raise ValueError('beta_v must be in [0,1]')
        except ValueError as e:
            self.logger.exception('beta_v must be in [0,1]')
            raise e
