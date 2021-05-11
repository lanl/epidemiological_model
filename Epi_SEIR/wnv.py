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

    """

    def __init__(self, config_file):
        self.logger = create_logger(__name__, vbdm.args.config_file)

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
        # NOTE: wnv system not modeled yet in this function

        # TODO Implement actual WNV system of equations

        # TODO Update parameter comment block

        # TODO update config file with dummy values

        # TODO update population values input with dummy values

        # TODO uncomment import wnv in models_main module

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


        # ----- OLD EQUATIONS
        # N_h = sum([Sh, Eh, Iha, Ihs, Rh])
        # N_v = sum([Sv, Ev, Iv])

        # Biting rate
        #b = self.params['sigma_h'] * self.params['sigma_v'] / \
        #    (self.params['sigma_h'] * N_h + self.params['sigma_v'] * N_v)
        #b_h = b * N_v
        #b_v = b * N_h

        # Force of infecton
        #lambda_h = b_h * self.params['beta_h'] * Iv / N_v
        #lambda_v = b_v * self.params['beta_v'] * (Iha + Ihs) / N_h

        # System of equations
        #dSh = -lambda_h * Sh
        #dEh = lambda_h * Sh - self.params['nu_h'] * Eh
        #dIha = self.params['psi'] * self.params['nu_h'] * Eh - \
        #    self.params['gamma_h'] * Iha
        #dIhs = (1 - self.params['psi']) * self.params['nu_h'] * \
        #    Eh - self.params['gamma_h'] * Ihs
        #dRh = self.params['gamma_h'] * (Iha + Ihs)
        #dSv = -lambda_v * Sv - self.params['mu_v'] * Sv
        #dEv = lambda_v * Sv - self.params['nu_v'] * Ev - self.params['mu_v'] * Ev
        #dIv = self.params['nu_v'] * Ev - self.params['mu_v'] * Iv

        #return dSh, dEh, dIha, dIhs, dRh, dSv, dEv, dIv
        # -----
