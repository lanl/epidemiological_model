"""West Nile Virus SEIR Model Class

Contains class for West Nile Virus disease model. Inherits from
Vector Borne Disease Model class.

    Typical usage example:

    wnv = WNVSEIRModel(<config_file_path>)
    
Class method param_dict allows you to change parameters values from configuration file inputs 
by inputting a dictionary of different parameter values
    
    Typical usage example:
    
    wnv_param_dict = WNVSEIRModel.param_dict(<config_file_path>, <parameter dictionary>)
"""

from utils import create_logger
import vbdm
import numpy as np
import math


class WNVSEIRModel(vbdm.VectorBorneDiseaseModel):

    """Models the spread of WNV.

    Inherits from the VectorBorneDiseaseModel class. Solves ODE system
    of equations and plots the resulting curves.

    Attributes:
        logger: python logging object.\n
        long_state_names: more compartment values for output.\n

    """

    def __init__(self, config_file, param_dict = None):
        self.logger = create_logger(__name__, config_file)

        super().__init__(config_file, 'WNV')
        
    @classmethod
    def param_dict(cls, config_file, param_dict):
        """Inherit vbdm param_dict class method"""
        return super(WNVSEIRModel, cls).param_dict(config_file = config_file, disease_name = 'WNV', param_dict =  param_dict)

    def _population_sizes(self):
        """Calculates population sizes of human and vector compartments"""
        self.Nv = sum([self.states['Sv'], self.states['Ev'], self.states['Iv']])
        self.Nb = sum([self.states['Sb'], self.states['Ib']])

    def _force_of_infection(self):
        """Calculates force of infection"""
        self.lambda_v = self.params['alpha_v']*self.params['beta_b']/self.Nb
        self.lambda_b = self.params['alpha_b']*self.params['beta_b']/self.Nb
    
    def _mosq_population_values(self, t):
        self.K_v = self.params['K_b'] + self.params['K_s'] * math.sin((2 * math.pi / 365) * t - math.pi / 2)
        self.r_v = self.params['r_b'] + self.params['r_s'] * math.sin((2 * math.pi / 365) * t - math.pi / 2)

    def model_func(self, t, y):
        """Defines system of ODEs for dengue model

        Initial State Names:
            Sv: Susceptible mosquito population.\n
            Ev: Exposed mosquito population. \n
            Iv: Infected mosquito population.\n
            Nv: Mosquito population.\n
            Sb: Susceptible bird population.\n
            Eb: Susceptible bird population.\n
            Ib: Infected bird population.\n
            Nb: Bird population.\n
            Ih: Infectious human population.\n

        Parameters:
            mu_v: Mosquito death rate.\n
            mu_b: Bird death/recovery rate. \n
            beta_b: Biting rate under frequency dependence. \n
            alpha_v: Probability of virus transmission to mosquito per infectious bite.\n
            alpha_b: Probability of virus transmission to bird,per infectious bit.\n
            K_b: For K_v function, baseline mosquito carrying capacity.\n
            K_s: For K_v function, scaling factor for the mosquito carrying capacity.\n
            K_v: Time varying carrying capacity for mosquitoes.\n
            r_b: For r_v function, baseline mosquito growth rate.\n
            r_s: For r_v function, scaling factor for mosquito growth rate.\n
            r_v: Time varying growth rate for mosquitoes.\n
            lambda_v: Mosquito force of infection. \n
            lambda_b: Bird force of infection. \n
            nu_v: Mosquito latent period. \n
            nu_b: Bird latent period. \n
            eta: Contact rate * probability of transmission to humans


        """
        ddt = self.initial_states.copy()

        self.states = dict(zip(self.initial_states.keys(), y))

        # Find population size
        self._population_sizes()

        # Find forces of infection
        self._force_of_infection()

        #Find mosquito carrying capacity and growth rate
        self._mosq_population_values(t)

        ddt['Sv'] = self.r_v * (1 - (self.Nv / self.K_v)) * self.Nv - \
            self.lambda_v * self.states['Sv'] * self.states['Ib'] - \
            self.params['mu_v'] * self.states['Sv']
        ddt['Ev'] = self.lambda_v * self.states['Sv'] * self.states['Ib'] - \
            self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Ev']
        ddt['Iv'] = self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Iv']
        ddt['Sb'] = -self.lambda_b * self.states['Iv'] * self.states['Sb']
        ddt['Eb'] = self.lambda_b * self.states['Iv'] * self.states['Sb'] - \
            self.params['nu_b'] * self.states['Eb']
        ddt['Ib'] = self.params['nu_b'] * self.states['Eb'] - \
            self.params['mu_b'] * self.states['Ib']


        rng = np.random.default_rng()
        try:
            ddt['Ih'] = rng.poisson(lam=self.params['eta'] * self.states['Iv'])
        except ValueError:
            self.logger.exception(f"Used Normal distribution, lam = {math.trunc(self.params['eta']*self.states['Iv'])}")
            ddt['Ih'] = math.trunc(rng.normal(loc = self.params['eta']*self.states['Iv'], scale = math.sqrt(self.params['eta']*self.states['Iv'])))

        return tuple(ddt.values())
