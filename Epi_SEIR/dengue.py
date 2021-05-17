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

        super().__init__(config_file, 'DENGUE')

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
        ddt = self.initial_states.copy()
        states = dict(zip(self.initial_states.keys(), y))

        Nh = sum([states['Sh'], states['Eh'], states['Iha'], states['Ihs'], states['Rh']])
        Nv = sum([states['Sv'], states['Ev'], states['Iv']])

        # Biting rate
        b = self.params['sigma_h'] * self.params['sigma_v'] / \
            (self.params['sigma_h'] * Nh + self.params['sigma_v'] * Nv)
        b_h = b * Nv
        b_v = b * Nh

        # Force of infecton
        lambda_h = b_h * self.params['beta_h'] * states['Iv'] / Nv
        lambda_v = b_v * self.params['beta_v'] * (states['Iha'] + states['Ihs']) / Nh

        # System of equations
        ddt['Sh'] = -lambda_h * states['Sh']
        ddt['Eh'] = lambda_h * states['Sh'] - \
            self.params['nu_h'] * states['Eh']
        ddt['Iha'] = self.params['psi'] * self.params['nu_h'] * states['Eh'] - \
            self.params['gamma_h'] * states['Iha']
        ddt['Ihs'] = (1 - self.params['psi']) * self.params['nu_h'] * states['Eh'] - \
            self.params['gamma_h'] * states['Ihs']
        ddt['Rh'] = self.params['gamma_h'] * (states['Iha'] + states['Ihs'])
        ddt['Sv'] = -lambda_v * states['Sv'] - \
            self.params['mu_v'] * states['Sv']
        ddt['Ev'] = lambda_v * states['Sv'] - \
            self.params['nu_v'] * states['Ev'] - \
            self.params['mu_v'] * states['Ev']
        ddt['Iv'] = self.params['nu_v'] * states['Ev'] - \
            self.params['mu_v'] * states['Iv']

        return tuple(ddt.values())
