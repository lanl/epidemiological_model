"""West Nile Virus SEIR Model Class

Contains class for West Nile Virus disease model. Inherits from
Vector Borne Disease Model class.

    Typical usage example:

    wnv = WNVSEIRModel(<config_file_path>)
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
        logger: python logging object.
        long_state_names: more compartment values for output.
        day_counter: keeps track of time passed to determine alpha.

    """

    def __init__(self, config_file):
        self.logger = create_logger(__name__, vbdm.args.config_file)
        self.day_counter = 0

        super().__init__(config_file, 'WNV')

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
            mu_v: Mosquito birth/death rate.\n
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
        ddt = self.initial_states.copy()

        states = dict(zip(self.initial_states.keys(), y))

        self.day_counter += 1
        if self.day_counter >= 200:
            self.params['alpha'] = 0

        beta = self.params['A'] + (self.params['K'] - self.params['A']) / \
            (1 + math.exp(-self.params['r'] * (t - self.params['t0'])))

        ddt['Sv'] = self.params['mu_v'] * states['Nv'] - \
            beta * states['Sv'] * states['Ib'] / states['Nb'] - \
            (self.params['mu_v'] + self.params['alpha']) * states['Sv']
        ddt['Iv'] = beta * states['Sv'] * states['Ib'] / states['Nb'] - \
            self.params['mu_v'] * states['Iv'] + \
            self.params['alpha'] * states['Sv']
        ddt['Sb'] = -beta * states['Iv'] * states['Sb'] / states['Nb']
        ddt['Ib'] = beta * states['Iv'] * states['Sb'] / states['Nb'] - \
            states['Ib'] / self.params['delta_b']

        # rng = np.random.default_rng()
        # dIh = rng.poisson(lam=self.params['eta'] * Iv)
        # print('FLAG -------- dIh, Iv:', dIh, Iv)

        ddt['Ih'] = self.params['eta'] * states['Iv']

        return tuple(ddt.values())
