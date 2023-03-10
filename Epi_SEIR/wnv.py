"""West Nile Virus SEIR Model Class

Contains class for West Nile Virus disease model. Inherits from
FitModel and Vector Borne Disease Model class.

    Typical usage example:

    wnv = WNVSEIRModel(<config_file_path>)
    
Class method param_dict allows you to change parameters values from configuration file inputs 
by inputting a dictionary of different parameter values
    
    Typical usage example:
    
    wnv_param_dict = WNVSEIRModel.param_dict(<config_file_path>, <parameter dictionary>)
"""

from utils import create_logger
import fit
import numpy as np
import math


class WNVSEIRModel(fit.FitModel):

    """Models the spread of WNV.

    Inherits from the FitModel class. Specifies ODE system
    of equations.

    """

    def __init__(self, config_file, param_dict = None):
        self.logger = create_logger(__name__, config_file)
        if __name__ == "wnv_test":
            self.logger.disabled = True

        super().__init__(config_file, 'WNV')
        
        self.error_zero_constants()
        self.error_zero_to_one_params()
        
    @classmethod
    def param_dict(cls, config_file, param_dict):
        """Inherit vbdm param_dict class method"""
        return super(WNVSEIRModel, cls).param_dict(config_file = config_file, disease_name = 'WNV', param_dict =  param_dict)

    def _population_sizes(self):
        """Calculates population sizes of human and vector compartments"""
        self.Nv = sum([self.states['Sv'], self.states['Ev'], self.states['Iv']])
        self.Nb = sum([self.states['Sb'], self.states['Eb'], self.states['Ib'], self.params['m'] * self.states['Rb']])

    def _force_of_infection(self):
        """Calculates force of infection"""
        self.lambda_v = self.params['beta_v']*self.params['alpha']/self.Nb
        self.lambda_b = self.params['beta_b']*self.params['alpha']/self.Nb
    
    def _mosq_population_values(self, t):
        self.K_v = self.params['K'] - self.params['K_s'] * math.cos((2 * math.pi / 365) * t)
        self.r_v = self.params['r'] - self.params['r_s'] * math.cos((2 * math.pi / 365) * t)
    
    def _bird_population_values(self, t):
        self.psi_b = self.params['rho'] + self.params['rho_s'] * math.sin(self.params['theta'] * t)

    def _birth_rate(self):
        """Caclualtes vector natural birth rate"""
        self.psi_v = self.r_v + self.params['mu_v']   

    def model_func(self, t, y):
        """Defines system of ODEs for dengue model

        Initial State Names:
            Sv: Susceptible vector population.\n
            Ev: Exposed vector population. \n
            Iv: Infected vector population.\n
            Nv: Vector population.\n
            Sb: Susceptible bird population.\n
            Eb: Exposed bird population.\n
            Ib: Infected bird population.\n
            Rb: Recovered bird population. \n
            Nb: Bird population.\n
            Ih: Infectious human population.\n

        Parameters:
            mu_v: Vector death rate.\n
            mu_b: Bird natural death rate. \n
            alpha: Biting rate under frequency dependence. \n
            beta_v: Probability of virus transmission to vector, per infectious bite.\n
            beta_b: Probability of virus transmission to bird, per infectious bit.\n
            psi_v: Vector birth rate. \n
            psi_b: bird recruitment rate. \n
            K_v: Vector carrying capacity.\n
            r_v: Vector growth rate.\n
            lambda_v: vector force of infection. \n
            lambda_b: Bird force of infection. \n
            nu_v: Vector latent period. \n
            nu_b: Bird latent period. \n
            gamma_b: Bird recovery rate. \n
            eta: Contact rate * probability of transmission to humans. \n
            m: Proportion of recovered birds that survived the WNV infection.


        """
        ddt = self.initial_states.copy()

        self.states = dict(zip(self.initial_states.keys(), y))

        # Find population size
        self._population_sizes()

        # Find forces of infection
        self._force_of_infection()

        # Find vector carrying capacity and growth rate
        self._mosq_population_values(t)
        
        # Find vector natural birth rate
        self._birth_rate()
        
        # Find vector carrying capacity and growth rate
        self._bird_population_values(t)

        ddt['Sv'] = (self.psi_v - self.r_v * self.Nv / self.K_v) * self.Nv - \
             self.lambda_v * self.states['Sv'] * self.states['Ib'] - \
             self.params['mu_v'] * self.states['Sv']
        ddt['Ev'] = self.lambda_v * self.states['Sv'] * self.states['Ib'] - \
            self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Ev']
        ddt['Iv'] = self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Iv']
        ddt['Sb'] = self.psi_b * self.Nb - \
            self.lambda_b * self.states['Iv'] * self.states['Sb'] - \
            self.params['mu_b'] * self.states['Sb']
        ddt['Eb'] = self.lambda_b * self.states['Iv'] * self.states['Sb'] - \
            self.params['nu_b'] * self.states['Eb'] - \
            self.params['mu_b'] * self.states['Eb']
        ddt['Ib'] = self.params['nu_b'] * self.states['Eb'] - \
            self.params['gamma_b'] * self.states['Ib'] - \
            self.params['mu_b'] * self.states['Ib'] 
        ddt['Rb'] = self.params['gamma_b'] * self.states['Ib'] - \
            self.params['mu_b'] * self.states['Rb']

        return tuple(ddt.values())
    
    
    def error_zero_constants(self):
        """check if human population is non-zero.

        check if vector carrying capcity is non-zero.
        """

        # EXTRACT list of position
        try:
            if sum([self.initial_states['Sb'], self.initial_states['Eb'], self.initial_states['Ib'], self.params['m'] * self.initial_states['Rb']]) == 0:
                raise ValueError('Bird population size cannot be zero')
        except ValueError as e:
            self.logger.exception('Bird population size cannot be zero')
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
            if self.params['nu_b'] < 0 or self.params['nu_b'] > 1:
                raise ValueError('nu_b must be in [0,1]')
        except ValueError as e:
            self.logger.exception('nu_b must be in [0,1]')
            raise e
        
        try:
            if self.params['mu_v'] < 0 or self.params['mu_v'] > 1:
                raise ValueError('mu_v must be in [0,1]')
        except ValueError as e:
            self.logger.exception('mu_v must be in [0,1]')
            raise e
        
        try:
            if self.params['mu_b'] < 0 or self.params['mu_b'] > 1:
                raise ValueError('mu_b must be in [0,1]')
        except ValueError as e:
            self.logger.exception('mu_b must be in [0,1]')
            raise e
        
        try:
            if self.params['alpha'] < 0 or self.params['alpha'] > 1:
                raise ValueError('alpha must be in [0,1]')
        except ValueError as e:
            self.logger.exception('alpha must be in [0,1]')
            raise e
            
        try:
            if self.params['eta'] < 0 or self.params['eta'] > 1:
                raise ValueError('eta must be in [0,1]')
        except ValueError as e:
            self.logger.exception('eta must be in [0,1]')
            raise e

