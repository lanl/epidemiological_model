"""Dengue SEIR Model Class

Contains class for dengue disease model. Inherits from
Vector Borne Disease Model class.

    Typical usage example:

    den = DengueSEIRModel(<config_file_path>, <command_line_arguments>)
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
        Nh: human population size
        Nh: vector population size

    """

    def __init__(self, config_file):
        self.logger = create_logger(__name__, config_file)

        super().__init__(config_file, 'DENGUE')

    def _population_sizes(self):
        """Calculates population sizes of human and vector compartments"""
        self.Nh = sum([self.states['Sh'], self.states['Eh'], self.states['Iha'], self.states['Ihs'], self.states['Rh']])
        self.Nv = sum([self.states['Sv'], self.states['Ev'], self.states['Iv']])

    def _biting_rate(self):
        """Calculates biting rate"""
        b = self.params['sigma_h'] * self.params['sigma_v'] / \
            (self.params['sigma_h'] * self.Nh + self.params['sigma_v'] * self.Nv)
        self.b_h = b * self.Nv
        self.b_v = b * self.Nh

    def _force_of_infection(self):
        """Calculates force of infection"""
        self.lambda_h = self.b_h * self.params['beta_h'] * self.states['Iv'] / self.Nv
        self.lambda_v = self.b_v * self.params['beta_v'] * (self.states['Iha'] + self.states['Ihs']) / self.Nh

    def model_func(self, t, y):
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
            beta_v: Probability of host to vector transmission.\n
            b_h: Biting rate (1 / (day * human))\n
            b_v: Biting rate (1 / (day * mosquito))\n

        """
        ddt = self.initial_states.copy()
        self.states = dict(zip(self.initial_states.keys(), y))

        # Find population size
        self._population_sizes()

        # Find biting rate
        self._biting_rate()

        # Find force of infection
        self._force_of_infection()

        # System of equations
        ddt['Sh'] = -self.lambda_h * self.states['Sh']
        ddt['Eh'] = self.lambda_h * self.states['Sh'] - \
            self.params['nu_h'] * self.states['Eh']
        ddt['Iha'] = self.params['psi'] * self.params['nu_h'] * self.states['Eh'] - \
            self.params['gamma_h'] * self.states['Iha']
        ddt['Ihs'] = (1 - self.params['psi']) * self.params['nu_h'] * self.states['Eh'] - \
            self.params['gamma_h'] * self.states['Ihs']
        ddt['Rh'] = self.params['gamma_h'] * (self.states['Iha'] + self.states['Ihs'])
        ddt['Sv'] = -self.lambda_v * self.states['Sv'] - \
            self.params['mu_v'] * self.states['Sv']
        ddt['Ev'] = self.lambda_v * self.states['Sv'] - \
            self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Ev']
        ddt['Iv'] = self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Iv']

        return tuple(ddt.values())
