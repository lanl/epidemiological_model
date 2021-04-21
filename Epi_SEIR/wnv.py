
"""West Nile Virus SEIR Model Class

Contains class for West Nile Virus disease model. Inherits from
Vector Borne Disease Model class.

    Typical usage example:

    wnv = WNVSEIRModel(<config_file_path>, <duration_in_days>)
"""

import sys
from scipy.integrate import odeint
from utils import create_logger
import VBDM


class WNVSEIRModel(VBDM.VectorBorneDiseaseModel):

    """Models the spread of WNV.

    Inherits from the VectorBorneDiseaseModel class. Solves ODE system
    of equations and plots the resulting curves.

    """

    def __init__(self, config_file, days):
        super().__init__(config_file, 'WNV', days)
        self.initial_states['Sv'] = self.mosq[0]
        print("WNV initial", self.initial_states)

    def run_model(self):
        """Runs ODE solver to generate model output"""
        y0 = self.initial_states['Sh'], self.initial_states['Eh'], \
            self.initial_states['Iha'], self.initial_states['Ihs'], \
            self.initial_states['Rh'], self.initial_states['Sv'], \
            self.initial_states['Ev'], self.initial_states['Iv']

        try:
            self.model_output = odeint(self._model_WNV, y0,
                                       self.t, args=(self,))
        except Exception:
            logger.exception('Exception occured running WNV model')
            sys.exit(1)
            # self.success = False
        else:
            logger.info('WNV model run complete')

    def _model_WNV(self, y, t, p):
        """Defines system of ODEs for WNV model"""
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
        dSv = -lambda_v * Sh
        dEv = lambda_v * Sh - self.params['nu_v'] * Ev
        dIv = self.params['nu_v'] * Ev - self.params['mu_v'] * Iv

        return dSh, dEh, dIha, dIhs, dRh, dSv, dEv, dIv


if sys.argv[0].endswith('sphinx-build'):
    logger = create_logger(__name__, '_default_config.yaml')
else:
    logger = create_logger(__name__, VBDM.args.config_file)