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
        self.logger = create_logger(__name__, config_file)
        self.day_counter = 0

        super().__init__(config_file, 'WNV')

    def _population_sizes(self):
        self.Nv = sum([self.states['Sv'], self.states['Iv']])
        self.Nb = sum([self.states['Sb'], self.states['Ib']])

    def _force_of_infection(self, t):
        self.beta = self.params['A'] + (self.params['K'] - self.params['A']) / \
            (1 + math.exp(-self.params['r'] * (t - self.params['t0'])))

    def model_func(self, t, y):
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

        self.states = dict(zip(self.initial_states.keys(), y))

        # Find population size
        self._population_sizes()

        # Find force of infection
        self._force_of_infection(t)

        self.day_counter += 1
        if self.day_counter >= 200:
            self.params['alpha'] = 0

        ddt['Sv'] = self.params['mu_v'] * self.Nv - \
            self.beta * self.states['Sv'] * self.states['Ib'] / self.Nb - \
            (self.params['mu_v'] + self.params['alpha']) * self.states['Sv']
        ddt['Iv'] = self.beta * self.states['Sv'] * self.states['Ib'] / self.Nb - \
            self.params['mu_v'] * self.states['Iv'] + \
            self.params['alpha'] * self.states['Sv']
        ddt['Sb'] = -self.beta * self.states['Iv'] * self.states['Sb'] / self.Nb
        ddt['Ib'] = self.beta * self.states['Iv'] * self.states['Sb'] / self.Nb - \
            self.states['Ib'] / self.params['delta_b']

        rng = np.random.default_rng()
        try:
            ddt['Ih'] = rng.poisson(lam=self.params['eta'] * self.states['Iv'])
        except ValueError:
            self.logger.exception(f"Used Normal distribution, lam = {math.trunc(self.params['eta']*self.states['Iv'])}")
            ddt['Ih'] = math.trunc(rng.normal(loc = self.params['eta']*self.states['Iv'], scale = math.sqrt(self.params['eta']*self.states['Iv'])))

        return tuple(ddt.values())
