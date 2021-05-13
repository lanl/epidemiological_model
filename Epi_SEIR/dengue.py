"""Dengue SEIR Model Class

Contains class for dengue disease model. Inherits from
Vector Borne Disease Model class.

    Typical usage example:

    den = DengueSEIRModel(<config_file_path>)
"""

from utils import create_logger
import vbdm


class DengueSEIRModel(vbdm.VectorBorneDiseaseModel):

    """Models the spread of dengue.

    Inherits from the VectorBorneDiseaseModel class. Solves ODE system
    of equations and plots the resulting curves.

    Attributes:
        logger: python logging object.
        long_state_names: more compartment values for output.

    """

    def __init__(self, config_file):
        self.logger = create_logger(__name__, vbdm.args.config_file)
        #self.long_state_names = ['Susceptible Humans', 'Exposed Humans',
        #                         'Asymptomatic Infected Humans', 'Symptomatic Infected Humans',
        #                         'Recovered Humans', 'Susceptible Vectors',
        #                         'Exposed Vectors', 'Infected Vectors']

        super().__init__(config_file, 'DENGUE')

    def set_y0(self):
        """Sets initial states to be passed into model_func"""
        y0 = self.initial_states['Sh'], self.initial_states['Eh'], \
            self.initial_states['Iha'], self.initial_states['Ihs'], \
            self.initial_states['Rh'], self.initial_states['Sv'], \
            self.initial_states['Ev'], self.initial_states['Iv']

        return y0

    def model_func(self, y, t):
        """Defines system of ODEs for dengue model.

        Initial State Names:
            Sh: Susceptible human population.\n
            Eh: Exposed human population.\n
            Iha: Asymptomatic infectious human population.\n
            Ihs: Symptomatic infectious human population\n
            Rh: Recovered human population.\n
            Sv: Susceptible vector population.\n
            Ev: Exposed vector population.\n
            Iv: Infectious vector population.\n

        Parameters:
            lambda_h: Human force of infection.\n
            lambda_v: Vector force of infection.\n
            nu_h: Human latent period.\n
            nu_v: Vector latent period.\n
            psi: Proportion of infections reported.\n
            gamma_h: Human infectious period.\n
            mu_v: Vector life expectancy.\n
            sigma_h: Maximum number of bites a human can support per unit time.\n
            sigma_v: Maximum vector biting rate.\n
            beta_h: Probability of vector to host transmission.\n
            beta_v: Probability of host tp vector transmission.\n
            b_h: Biting rate (1 / (day * human))\n
            b_v: Biting rate (1 / (day * mosquito))\n

        """
        # States and population
        Sh, Eh, Iha, Ihs, Rh, Sv, Ev, Iv = y
        N_h = sum([Sh, Eh, Iha, Ihs, Rh])
        N_v = sum([Sv, Ev, Iv])

        # Biting rate
        b = self.params['sigma_h'] * self.params['sigma_v'] / \
            (self.params['sigma_h'] * N_h + self.params['sigma_v'] * N_v)
        b_h = b * N_v
        b_v = b * N_h

        # Force of infecton
        lambda_h = b_h * self.params['beta_h'] * Iv / N_v
        lambda_v = b_v * self.params['beta_v'] * (Iha + Ihs) / N_h

        # System of equations
        dSh = -lambda_h * Sh
        dEh = lambda_h * Sh - self.params['nu_h'] * Eh
        dIha = self.params['psi'] * self.params['nu_h'] * Eh - \
            self.params['gamma_h'] * Iha
        dIhs = (1 - self.params['psi']) * self.params['nu_h'] * \
            Eh - self.params['gamma_h'] * Ihs
        dRh = self.params['gamma_h'] * (Iha + Ihs)
        dSv = -lambda_v * Sv - self.params['mu_v'] * Sv
        dEv = lambda_v * Sv - self.params['nu_v'] * Ev - self.params['mu_v'] * Ev
        dIv = self.params['nu_v'] * Ev - self.params['mu_v'] * Iv

        return dSh, dEh, dIha, dIhs, dRh, dSv, dEv, dIv
