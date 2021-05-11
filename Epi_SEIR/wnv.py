"""West Nile Virus SEIR Model Class

Contains class for West Nile Virus disease model. Inherits from
Vector Borne Disease Model class.

    Typical usage example:

    wnv = WNVSEIRModel(<config_file_path>)
"""

from utils import create_logger
import vbdm
#import numpy as np
import math


class WNVSEIRModel(vbdm.VectorBorneDiseaseModel):

    """Models the spread of WNV.

    Inherits from the VectorBorneDiseaseModel class. Solves ODE system
    of equations and plots the resulting curves.

    Attributes:
        logger: python logging object.
        long_state_names: more compartment values for output.

    """

    def __init__(self, config_file):
        self.logger = create_logger(__name__, vbdm.args.config_file)
        self.long_state_names = ['Susceptible Vectors', 'Infected Vectors',
                                 'Vector Population Size', 'Susceptible Birds',
                                 'Infected Birds', 'Bird Population Size',
                                 'Infected Humans']

        super().__init__(config_file, 'WNV')

    def set_y0(self):
        """Sets initial states to be passed into model_func"""
        y0 = self.initial_states['Sv'], self.initial_states['Iv'], \
            self.initial_states['Nv'], self.initial_states['Sb'], \
            self.initial_states['Ib'], self.initial_states['Nb'], \
            self.initial_states['Ih']

        return y0

    def model_func(self, y, t):
        """Defines system of ODEs for dengue model

        Initial State Names:
            Sv: Susceptible mosquito population.\n
            Iv: Infected mosquito population.\n
            Nv: Mosquito population.\n
            Sb: Susceptible bird population.\n
            Ib: Infected bird population.\n
            Nb: Bird population.\n
            Ih: Infectious human population.\n

        Parameters:
            mu_m: Mosquito birth/death rate.\n
            beta: Contact rate, probability of transmission between birds and
                  mosquitoes at time t. \n
            alpha: Rate of WNV seeding into the local model domain before day
                   200.\n
            delta_b: Recovery rate of birds.\n
            eta: Contact rate, probability of transmission from mosquitoes to
                 humans.\n

            A: Lower asymptote for beta.\n
            K: Upper asymptote for beta.\n
            r: Growth rate for beta.\n
            t0: Inflection point for beta.\n

        """
        # States and population
        Sv, Iv, Nv, Sb, Ib, Nb, Ih = y

        beta = self.params['A'] + (self.params['K'] - self.params['A']) / \
            (1 + math.exp(-self.params['r'] * (t - self.params['t0'])))

        dSv = self.params['mu_m'] * Nv - beta * Sv * Ib / Nb - \
            (self.params['mu_m'] + self.params['alpha']) * Sv
        dIv = beta * Sv * Ib / Nb - self.params['mu_m'] * Iv + \
            self.params['alpha'] * Sv
        dSb = -beta * Iv * Sb / Nb
        dIb = beta * Iv * Sb / Nb - Ib / self.params['delta_b']
        dIh = self.params['eta'] * Iv  # TODO wrap in poisson random number generator

        return dSv, dIv, Nv, dSb, dIb, Nb, dIh
