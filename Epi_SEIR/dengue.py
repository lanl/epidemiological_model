"""Dengue SEIR Model Class

Contains class for dengue disease model. Inherits from
Vector Borne Disease Model class.

    Typical usage example:

    den = DengueSEIRModel(<config_file_path>)
"""

import sys
from scipy.integrate import odeint
from utils import create_logger
import VBDM
import numpy as np


class DengueSEIRModel(VBDM.VectorBorneDiseaseModel):

    """Models the spread of dengue.

    Inherits from the VectorBorneDiseaseModel class. Solves ODE system
    of equations and plots the resulting curves.

    """

    def __init__(self, config_file):
        self.logger = create_logger(__name__, VBDM.args.config_file)

        super().__init__(config_file, 'DENGUE')
        self.initial_states['Sv'] = self.mosq[0]
        # print("DENGUE INITIAL", self.initial_states)

    def run_model(self):
        """Runs ODE solver to generate model output"""
        y0 = self.initial_states['Sh'], self.initial_states['Eh'], \
            self.initial_states['Iha'], self.initial_states['Ihs'], \
            self.initial_states['Rh'], self.initial_states['Sv'], \
            self.initial_states['Ev'], self.initial_states['Iv']

        t = np.linspace(0, 1, self.config_dict['RESOLUTION'] + 1)

        try:
            self.model_output = odeint(self._model_dengue, y0,
                                       t, args=(self,))
        except Exception:
            self.logger.exception('Exception occured running dengue model')
            sys.exit(1)

        for i in range(1, self.config_dict['DURATION']):
            self.initial_states['Sv'] = self.mosq[i]
            y0 = tuple(self.model_output[-1])
            self.model_output = self.model_output[:-1]

            try:
                out = odeint(self._model_dengue, y0,
                             t, args=(self,))
            except Exception:
                self.logger.exception('Exception occured running dengue model')
                sys.exit(1)

            self.model_output = np.concatenate((self.model_output, out))

    def _model_dengue(self, y, t, p):
        """Defines system of ODEs for dengue model"""
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
