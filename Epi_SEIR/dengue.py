"""Dengue SEIR Model Class

Contains class for dengue disease model. Inherits from
Vector Borne Disease Model class.

    Typical usage example:

    den = DengueSEIRModel(<config_file_path>, <command_line_arguments>)
    
Class method param_dict allows you to change parameters values from configuration file inputs 
by inputting a dictionary of different parameter values
    
    Typical usage example:
    
    den_param_dict = DengueSEIRModel.param_dict(<config_file_path>, <parameter dictionary>)
"""

from utils import create_logger
import vbdm


class DengueSEIRModel(vbdm.VectorBorneDiseaseModel):

    """Models the spread of dengue.

    Inherits from the VectorBorneDiseaseModel class. Solves ODE system
    of equations and plots the resulting curves.

    Attributes:
        logger: python logging object.\n
        long_state_names: more compartment values for output.\n
        Nh: human population size.\n
        Nh: vector population size.\n

    """

    def __init__(self, config_file, param_dict = None): 
        self.logger = create_logger(__name__, config_file)
        
        super().__init__(config_file, 'DENGUE')
    
    @classmethod
    def param_dict(cls, config_file, param_dict):
        """Inherit vbdm param_dict class method"""
        return super(DengueSEIRModel, cls).param_dict(config_file = config_file, disease_name = 'DENGUE', param_dict =  param_dict)
    
    def _population_sizes(self):
        """Calculates population sizes of human and vector compartments"""
        self.Nh = sum([self.states['Sh'], self.states['Eh'], self.states['Ih'], self.states['Rh']])
        self.Nv = sum([self.states['Sv'], self.states['Ev'], self.states['Iv']])
    
    #for complex force of infection: will ultimately use this later
    #def _biting_rate(self):
        #"""Calculates biting rate"""
        #b = self.params['sigma_h'] * self.params['sigma_v'] / \
            #(self.params['sigma_h'] * self.Nh + self.params['sigma_v'] * self.Nv)
        #self.b_h = b * self.Nv
        #self.b_v = b * self.Nh
    
    #complex force of infection: will ultimately use this later
    #def _force_of_infection(self):
        #"""Calculates force of infection"""
        #self.lambda_h = self.b_h * self.params['beta_h'] * self.states['Iv'] / self.Nv
        #self.lambda_v = self.b_v * self.params['beta_v'] * (self.states['Ih']) / self.Nh
    
    def _force_of_infection(self):
        """Calculates force of infection"""
        self.lambda_h = self.params['a_v'] * self.params['beta_h'] * self.states['Iv'] / self.Nh
        self.lambda_v = self.params['a_v'] * self.params['beta_v'] * (self.states['Ih']) / self.Nh
        
    def _birth_rate(self):
        """Caclualtes vector natural birth rate"""
        self.psi_v = self.params['r_v'] + self.params['mu_v']

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
            gamma_h: Human infectious period.\n
            psi_v: Vector natural birth rate .\n
            mu_v: Vector natural death rate.\n
            r_v: Vector instrinsic growth rate.\n
            K_v: Vector carrying capacity.\n
            a_v: Average number of bites to a human for each mosquito per day.\n
            beta_h: Probability of vector to host transmission.\n
            beta_v: Probability of host to vector transmission.\n

        """
        ddt = self.initial_states.copy()
        self.states = dict(zip(self.initial_states.keys(), y))

        # Find population size
        self._population_sizes()

        # Find biting rate
        #self._biting_rate()

        # Find force of infection
        self._force_of_infection()
        
        # Find vector natural birth rate
        self._birth_rate()

        # System of equations
        ddt['Sh'] = -self.lambda_h * self.states['Sh']
        ddt['Eh'] = self.lambda_h * self.states['Sh'] - \
            self.params['nu_h'] * self.states['Eh']
        ddt['Ih'] = self.params['nu_h'] * self.states['Eh'] - \
            self.params['gamma_h'] * self.states['Ih']
        ddt['Rh'] = self.params['gamma_h'] * self.states['Ih']
        ddt['Sv'] = (self.psi_v - self.params['r_v'] * self.Nv / self.params['K_v']) * self.Nv - \
            self.lambda_v * self.states['Sv'] - \
            self.params['mu_v'] * self.states['Sv']
        ddt['Ev'] = self.lambda_v * self.states['Sv'] - \
            self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Ev']
        ddt['Iv'] = self.params['nu_v'] * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Iv']

        return tuple(ddt.values())
