"""Dengue SEIR Model Class

Contains class for dengue disease model. Inherits from
FitModel and Vector Borne Disease Model class.

    Typical usage example:

    den = DengueSEIRModel(<config_file_path>, <command_line_arguments>)
    
Class method param_dict allows you to change parameters values from configuration file inputs 
by inputting a dictionary of different parameter values
    
    Typical usage example:
    
    den_param_dict = DengueSEIRModel.param_dict(<config_file_path>, <parameter dictionary>)
"""

from utils import create_logger
import sys
import math
import fit
import os
import pandas as pd
import numpy as np
from scipy.interpolate import BSpline
from scipy import interpolate

class DengueSEIRModel(fit.FitModel):

    """Models the spread of dengue.

    Inherits from the FitModel class. Specifies ODE system
    of equations.

    """

    def __init__(self, config_file, param_dict = None):
        self.logger = create_logger(__name__, config_file)
        if __name__ == "dengue_test":
            self.logger.disabled = True
        
        super().__init__(config_file, 'DENGUE')
        
        self.error_zero_constants()
        self.error_zero_to_one_params()
        #------------------------------------------------------------
        # Temperature paths reading
 
        #self.Temp = np.array(pd.read_csv(self.config_dict['TEMPERATURE_FILE_PATH'])['t2m']) #added temperature
        
        
        
        
        #------------------------------------------
        #automatic code
        
        self.K= np.array(pd.read_csv(self.config_dict['CARRYING_CAPACITY'])['K']) # carrying capacity
        self.K_s= np.array(pd.read_csv(self.config_dict['CARRYING_CAPACITY'])['K_s']) #initial carrying capac
        self.r= np.array(pd.read_csv(self.config_dict['CARRYING_CAPACITY'])['r']) #mosquito growth rate
        self.r_s= np.array(pd.read_csv(self.config_dict['CARRYING_CAPACITY'])['r_s']) #initial growth rate 
        self.year= np.array(pd.read_csv(self.config_dict['CARRYING_CAPACITY'])['year']) #data year
        self.population=np.array(pd.read_csv(self.config_dict['CARRYING_CAPACITY'])['Population']) #population 
        self.year_number=self.config_dict['YEAR']
        self.municipality_name=self.config_dict['MUNICIPALITY']
        #----------------------------------
        #SPLINE SECTION
        
        #Initial biting rates
        #self.initial_biting_rates=np.random.uniform(low=0, high=0.5, size=(self.params['degrees_freedom']+1,)).tolist()
        # Define knots
        #print(self.config_dict['DURATION'])
        self.psis_knots=[(((i+1)-1)/(self.params['degrees_freedom']))*self.config_dict['DURATION'] for i in range (0,self.params['degrees_freedom']+1)]
        #print (self.initial_biting_rates)
        #print (self.psis_knots)
        
        #----------------------------------
    
    @classmethod
    def param_dict(cls, config_file, param_dict):
        """Inherit vbdm param_dict class method"""
        return super(DengueSEIRModel, cls).param_dict(config_file = config_file, disease_name = 'DENGUE', param_dict =  param_dict)
    
    def _population_sizes(self):
        """Calculates population sizes of human and vector compartments"""
        self.Nh = sum([self.states['Sh'], self.states['Eh'], self.states['Ih'], self.states['Rh']])
        self.Nv = sum([self.states['Sv'], self.states['Ev'], self.states['Iv']])
        
        
        #-----------------------------------------------------------------------------
    #=============================================================================
    # SPLINE SECTION MODIFICATION
    
    def Spline_1storder(self, index, t):
        if index+1==1:
            return 1
        elif t<self.psis_knots[index-1]:
            return 0
        elif self.psis_knots[index-1]<=t and t<self.psis_knots[index]:
            return t-self.psis_knots[index-1]
        elif self.psis_knots[index] <=t:
            return self.psis_knots[index]-self.psis_knots[index-1]
        
        
    def Spline_2ndorder(self, index, t):
        if index+1==1:
            return 1
        elif index+1==2:
            return t
        elif t<self.psis_knots[index]:
            return 0
        elif self.psis_knots[index]<=t and t<self.psis_knots[index+1]:
            return (t-self.psis_knots[index-2])**2
        elif self.psis_knots[index+1] <=t:
            return 2*(self.psis_knots[index-1]-self.psis_knots[index-2])*(t-self.psis_knots[index-1])+(self.psis_knots[index-1]-self.psis_knots[index-2])**2
                  
    
    
    def Splines_piecewise(self, t):
        #self.sp=0;
        #-----------------------------------------
        # Spline first order
        #self.sp=self.params['beta1']*self.Spline_1storder(0,t)+self.params['beta2']*self.Spline_1storder(1,t)+\
        #self.params['beta3']*self.Spline_1storder(2,t)
        #-------------------------------------------
        # Spline second order
        #self.sp=self.params['beta1']*self.Spline_2ndorder(0,t+self.params['intercept'])+\
        #self.params['beta2']*self.Spline_2ndorder(1,t+self.params['intercept'])+\
        #self.params['beta3']*self.Spline_2ndorder(2,t+self.params['intercept'])
        #--------------------------------------------
        # BSpline Python
        #self.coeff=[self.params['beta1'], self.params['beta2'],self.params['beta3'], self.params['beta4']]
        #self.sp= BSpline(self.psis_knots,self.coeff,3,extrapolate='periodic')
        #==================================================
        # BSpline Python matching ends
        self.ys=[self.params['y1'], self.params['y2'], self.params['y3'], self.params['y4'],self.params['y5'], self.params['y6'],self.params['y7'],self.params['y8'], self.params['y9'], self.params['y10'],self.params['y11'],
                self.params['y12']]
        self.xs=np.linspace(0,self.config_dict['DURATION'],len(self.ys))
        self.tck= interpolate.splrep(self.xs, self.ys, s=0, k=3)
        self.sp=BSpline(self.tck[0],self.tck[1],self.tck[2])
        
        
        return self.sp(t)
#         for j in range(0,self.params['degrees_freedom']):
#             self.sp=self.sp+coeff_beta[j]*self.Spline_1storder(j,t)
        #return self.sp
        
    def _biting_rate(self,t):
        """Calculates biting rate"""
        #self.b = self.params['sigma_h'] * self.Splines_piecewise(t) / \
        #    (self.params['sigma_h'] * self.Nh + self.Splines_piecewise(t) * self.Nv)
        #-------------------------------------------------------------------------------
        # BSPline Python
        #if t==0:
        #    self.b = self.params['sigma_h']*(-0.00766667) / \
        #    (self.params['sigma_h']*self.Nh + (-0.00766667)*self.Nv)
        #else:
        self.b = self.params['sigma_h']*self.Splines_piecewise(t+self.params['intercept']) / \
        (self.params['sigma_h']*self.Nh + self.Splines_piecewise(t+self.params['intercept'])*self.Nv)
                
            
#    
#             else
#                 self.b = self.params['sigma_h']*self.Splines_piecewise(t+self.params['intercept']) / \
#                 (self.params['sigma_h']*self.Nh + self.Splines_piecewise(t+self.params['intercept'])*self.Nv)  
          
#         self.b = self.params['sigma_h']*self.Splines_piecewise(t+self.params['intercept']) / \
#            (self.params['sigma_h']*self.Nh + self.Splines_piecewise(t+self.params['intercept'])*self.Nv)
        
        
        
    
    def _force_of_infection(self):
       """Calculates force of infection"""
       self.lambda_h = self.b * self.params['beta_h'] * self.states['Iv'] 
       self.lambda_v = self.b * self.params['beta_v'] * self.states['Ih']
    
   
    #-----------------------------------------------------------------------------
    #===============================================================================
    
#     def _biting_rate(self):
#         """Calculates biting rate"""
#         self.b = self.params['sigma_h'] * self.params['sigma_v'] / \
#             (self.params['sigma_h'] * self.Nh + self.params['sigma_v'] * self.Nv)
    
#     def _force_of_infection(self):
#         """Calculates force of infection"""
#         self.lambda_h = self.b * self.params['beta_h'] * self.states['Iv'] 
#         self.lambda_v = self.b * self.params['beta_v'] * self.states['Ih']
    #===============================================================================
    #def _birth_rate(self,t):
        #"""Caclualtes vector natural birth rate"""
        #if t < 180:
            #self.psi_v = self.params['r_v'] + self.params['mu_v']
        #else:
            #self.psi_v = 0
            #self.r_v = self.psi_v - self.params['mu_v']
    
    
    def _mosq_population_values(self, t):
        
        #-------------------------------------------------------------------------
#         cont=-1
#         if t%365==0:
#             cont+=1
#             string+=
#         else:
#              self.K_v = self.params['K'] - self.params['K_s'] * math.cos((2 * math.pi*t / 365))
#              self.r_v = self.params['r'] - self.params['r_s'] * math.cos((2 * math.pi*t / 365))
            
        #t=list(range(20))
        
#         year_dict={0: 2000, 1: 2001, 2: 2002, 3: 2003,4: 2004 ,5: 2005, 6: 2006,7: 2007,\
#                    8: 2008, 9: 2009, 10: 2010 ,11: 2011, 12: 2012 , 13: 2013, 14: 2014,\
#                    15: 2015, 16: 2016, 17: 2017, 18: 2018}
          
#         #for i in range(0,len(t)):
        
#         if (t%365==0):
#             print("New year",self.year_number)
#             print(t)
#             self.K_v = float(self.K[self.year_number]) - float(self.K_s[self.year_number]) * math.cos((2 * math.pi*t / 365))
#             self.r_v = float(self.r[self.year_number]) - float(self.r_s[self.year_number]) * math.cos((2 * math.pi*t / 365))
#             self.year_number+=1
#         else:
#             self.K_v = float(self.K[self.year_number-1]) - float(self.K_s[self.year_number-1]) * math.cos((2 * math.pi*t / 365))
#             self.r_v = float(self.r[self.year_number-1]) - float(self.r_s[self.year_number-1]) * math.cos((2 * math.pi*t / 365))
#             print("Still in year",self.year_number-1)
#             print(t)
        
        
           
        if (t>=365 and t<730):
            self.year_number=1
          
        if (t>=730 and t<1095):
            self.year_number=2
          
        if (t>=1095 and t<1460):
            self.year_number=3
     
        if (t>=1460 and t<1825):
            self.year_number=4
     
        if (t>=1825 and t<2190):
            self.year_number=5
     
        if (t>=2190 and t<2555):
            self.year_number=6
       
        if (t>=2555 and t<2920):
            self.year_number=7
   
        if (t>=2920 and t<3285):
            self.year_number=8
      
        if (t>=3285 and t<3650):
            self.year_number=9
  
        if (t>=3650 and t<4015):
            self.year_number=10
 
        if (t>=3650 and t<4015):
            self.year_number=11
 
        if (t>=4015 and t<4380):
            self.year_number=12
      
        if (t>=4380 and t<4745):
            self.year_number=13
    
        if (t>=4745 and t<5110):
            self.year_number=14
      
        if (t>=5110 and t<5475):
            self.year_number=15
       
        if (t>=5475 and t<5840):
            self.year_number=16
  
        if (t>=5840 and t<6205):
            self.year_number=17
    
        if (t>=6205 and t<6570):
            self.year_number=18
            
        self.K_v = float(self.K[self.year_number])*float(self.human_pop[self.year_number]) - float(self.K_s[self.year_number])*float(self.human_pop[self.year_number]) * math.cos((2 * math.pi*t / 365))
        self.r_v = float(self.r[self.year_number]) - float(self.r_s[self.year_number]) * math.cos((2 * math.pi*t / 365))
        #print(self.K_v)
        #print(self.r_v)
        #-------------------------------------------------------------------------
        #general code
        #self.K_v = self.params['K'] - self.params['K_s'] * math.cos((2 * math.pi*t / 365))
        #self.r_v = self.params['r'] - self.params['r_s'] * math.cos((2 * math.pi*t / 365))
     
#============================================================ 
#AUTOMATIC CODE YEAR
    def _mosq_population_values_year(self, t, year):
           
            
        self.K_v = float(self.K[year])*float(self.human_pop[year]) - float(self.K_s[year])*float(self.human_pop[year]) * math.cos((2 * math.pi*t / 365))
        self.r_v = float(self.r[year]) - float(self.r_s[year]) * math.cos((2 * math.pi*t / 365))
        #print (self.K[year])
        
    #============================================================     
    
    def _initial_human_pop(self):
        self.H0 = sum([self.initial_states['Sh'], self.initial_states['Eh'], self.initial_states['Ih'], self.initial_states['Rh']])
        print (self.initial_states['Sh'])
        print (self.initial_states['Eh'])
        print (self.initial_states['Ih'])
        print (self.initial_states['Rh'])
    
    def _birth_rate(self):
        """Calculates vector natural birth rate"""
        self.psi_v = self.r_v + self.params['mu_v']
    
    def _calc_hNv(self):
        self.hNv = self.psi_v - self.r_v * self.Nv / self.K_v
        
    #========================================================
    #AUTOMATIC CODE 
    def calculate_current_temp_year(self,t,year):
        #self.currentT= self.Temp[int(t)]
        temperature_path = os.path.join(self.config_dict['TEMPERATURE_FILE_PATH_YEAR'],
                                               f'{year+2000}.csv')
        temperature_df=pd.read_csv(temperature_path,sep=',')
        temperature_df=temperature_df[(temperature_df['NM_MUN']==self.municipality_name)]
        self.Temp = np.array(temperature_df['t2m'])
        
        self.EIP_Temp= np.exp(0.8-0.20*self.Temp[int(t)])
    #========================================================
    def calculate_current_temp(self,t):
        #self.currentT= self.Temp[int(t)]
        self.EIP_Temp= np.exp(0.8-0.20*self.Temp[int(t)])
        
    #def model_func(self, t, y):
    
    
    #=======================================================
    #AUTOMATIC CODE
    def model_func_year(self, t, y, year):
        """Defines system of ODEs for dengue model.

        Initial State Names:
            Sh: Susceptible human population.\n
            Eh: Exposed human population.\n
            Ih: Infectious human population.\n
            Rh: Recovered human population.\n
            Nh: Human population.\n
            Sv: Susceptible vector population.\n
            Ev: Exposed vector population.\n
            Iv: Infectious vector population.\n
            Nv: Vector population.\n

        Parameters:
            H0: Initial human population.\n
            psi_h: Per capita birth rate for humans.\n
            psi_v: Per capita birth rate for mosquitoes.\n
            beta_h: Probability of vector to host transmission.\n
            beta_v: Probability of host to vector transmission.\n
            sigma_h: Maximum number of times one mosquito bite a human per day.\n
            sigma_v: Maximum number of mosquito bites a human can sustain per day.\n
            lambda_h: Human force of infection.\n
            lambda_v: Vector force of infection.\n
            alpha_h: Rate of lost dengue immunity.\n
            nu_h: Human latent period.\n
            nu_v: Vector latent period.\n
            gamma_h: Human infectious period.\n
            mu_v: Vector natural death rate.\n
            mu_h: Human natural death rate.\n
            r: Baseline vector growth rate factor.\n
            r_s: Scaling vector growth rate factor.\n
            r_v: Time-varying vector growth rate.\n
            K: Baseline vector carrying capacity factor.\n
            K_s: Scaling vector carrying capacity factor.\n
            K_v: Time-varying vector carrying capacity.\n

        """
        ddt = self.initial_states.copy()
        self.states = dict(zip(self.initial_states.keys(), y))

        # Find population size
        self._population_sizes()

        # Find biting rate
        #self._biting_rate()
        self._biting_rate(t)

        # Find force of infection
        self._force_of_infection()

        # Find mosquito carrying capcity and growth rate from LLM fit
        self._mosq_population_values_year(t, year)
        
        # Find initial human population
        self._initial_human_pop()
      
        
        
        # Find vector natural birth rate
        self._birth_rate()
        
        # Find hNv value
        self._calc_hNv()
        
        self.calculate_current_temp_year(t, year) # ADDED TEMPERATURE ARRAY

        # System of equations
        ddt['Sh'] = self.params['psi_h'] * self.H0 - \
            self.lambda_h * self.states['Sh'] + \
            self.params['alpha_h'] * self.states['Rh'] - \
            self.params['mu_h'] * self.states['Sh']
        ddt['Eh'] = self.lambda_h * self.states['Sh'] - \
            self.params['nu_h'] * self.states['Eh'] - \
            self.params['mu_h'] * self.states['Eh']
        ddt['Ih'] = self.params['nu_h'] * self.states['Eh'] - \
            self.params['gamma_h'] * self.states['Ih'] - \
            self.params['mu_h'] * self.states['Ih']
        ddt['Rh'] = self.params['gamma_h'] * self.states['Ih'] - \
            self.params['alpha_h'] * self.states['Rh'] - \
            self.params['mu_h'] * self.states['Rh']
        ddt['Ch'] = self.params['nu_h'] * self.states['Eh']
        ddt['Sv'] = self.hNv * self.Nv - \
            self.lambda_v * self.states['Sv'] - \
            self.params['mu_v'] * self.states['Sv']
#         ddt['Ev'] = self.lambda_v * self.states['Sv'] - \
#            self.params['nu_v'] * self.states['Ev'] - \
#            self.params['mu_v'] * self.states['Ev']
#         ddt['Iv'] = self.params['nu_v'] * self.states['Ev'] - \
#            self.params['mu_v'] * self.states['Iv']
        
        
        
        ddt['Ev'] = self.lambda_v * self.states['Sv'] - \
            self.EIP_Temp* self.states['Ev'] - \
            self.params['mu_v'] * self.states['Ev']
        ddt['Iv'] = self.EIP_Temp * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Iv']
        
       
        return tuple(ddt.values())
    #==================================================================
    
    
    
    
    
    
    def model_func(self, t, y):
        """Defines system of ODEs for dengue model.

        Initial State Names:
            Sh: Susceptible human population.\n
            Eh: Exposed human population.\n
            Ih: Infectious human population.\n
            Rh: Recovered human population.\n
            Nh: Human population.\n
            Sv: Susceptible vector population.\n
            Ev: Exposed vector population.\n
            Iv: Infectious vector population.\n
            Nv: Vector population.\n

        Parameters:
            H0: Initial human population.\n
            psi_h: Per capita birth rate for humans.\n
            psi_v: Per capita birth rate for mosquitoes.\n
            beta_h: Probability of vector to host transmission.\n
            beta_v: Probability of host to vector transmission.\n
            sigma_h: Maximum number of times one mosquito bite a human per day.\n
            sigma_v: Maximum number of mosquito bites a human can sustain per day.\n
            lambda_h: Human force of infection.\n
            lambda_v: Vector force of infection.\n
            alpha_h: Rate of lost dengue immunity.\n
            nu_h: Human latent period.\n
            nu_v: Vector latent period.\n
            gamma_h: Human infectious period.\n
            mu_v: Vector natural death rate.\n
            mu_h: Human natural death rate.\n
            r: Baseline vector growth rate factor.\n
            r_s: Scaling vector growth rate factor.\n
            r_v: Time-varying vector growth rate.\n
            K: Baseline vector carrying capacity factor.\n
            K_s: Scaling vector carrying capacity factor.\n
            K_v: Time-varying vector carrying capacity.\n

        """
        ddt = self.initial_states.copy()
        self.states = dict(zip(self.initial_states.keys(), y))

        # Find population size
        self._population_sizes()

        # Find biting rate
        #self._biting_rate()
        self._biting_rate(t)

        # Find force of infection
        self._force_of_infection()

        # Find mosquito carrying capcity and growth rate from LLM fit
        self._mosq_population_values(t)
        
        # Find initial human population
        self._initial_human_pop()
        
        # Find vector natural birth rate
        self._birth_rate()
        
        # Find hNv value
        self._calc_hNv()
        
        self.calculate_current_temp(t) # ADDED TEMPERATURE ARRAY

        # System of equations
        ddt['Sh'] = self.params['psi_h'] * self.H0 - \
            self.lambda_h * self.states['Sh'] + \
            self.params['alpha_h'] * self.states['Rh'] - \
            self.params['mu_h'] * self.states['Sh']
        ddt['Eh'] = self.lambda_h * self.states['Sh'] - \
            self.params['nu_h'] * self.states['Eh'] - \
            self.params['mu_h'] * self.states['Eh']
        ddt['Ih'] = self.params['nu_h'] * self.states['Eh'] - \
            self.params['gamma_h'] * self.states['Ih'] - \
            self.params['mu_h'] * self.states['Ih']
        ddt['Rh'] = self.params['gamma_h'] * self.states['Ih'] - \
            self.params['alpha_h'] * self.states['Rh'] - \
            self.params['mu_h'] * self.states['Rh']
        ddt['Ch'] = self.params['nu_h'] * self.states['Eh']
        ddt['Sv'] = self.hNv * self.Nv - \
            self.lambda_v * self.states['Sv'] - \
            self.params['mu_v'] * self.states['Sv']
#         ddt['Ev'] = self.lambda_v * self.states['Sv'] - \
#            self.params['nu_v'] * self.states['Ev'] - \
#            self.params['mu_v'] * self.states['Ev']
#         ddt['Iv'] = self.params['nu_v'] * self.states['Ev'] - \
#            self.params['mu_v'] * self.states['Iv']
        
        
        
        ddt['Ev'] = self.lambda_v * self.states['Sv'] - \
            self.EIP_Temp* self.states['Ev'] - \
            self.params['mu_v'] * self.states['Ev']
        ddt['Iv'] = self.EIP_Temp * self.states['Ev'] - \
            self.params['mu_v'] * self.states['Iv']
        
       
        return tuple(ddt.values())
    
    
    def error_zero_constants(self):
        """check if human population is non-zero.

        check if vector carrying capcity is non-zero.
        """

        # EXTRACT list of position
        try:
            if sum([self.initial_states['Sh'], self.initial_states['Eh'], self.initial_states['Ih'], self.initial_states['Rh']]) == 0:
                raise ValueError('Human population size cannot be zero')
        except ValueError as e:
            self.logger.exception('Human population size cannot be zero')
            raise e
    
    def error_zero_to_one_params(self):
        """check if parameters that should be in [0,1] are
        """
        try:
            if self.params['nu_v'] < 0 or self.params['nu_v'] > 1:
                raise ValueError('nu_v must be in [0,1]')
        except ValueError as e:
            self.logger.exception('nu_v must be in [0,1]')
            raise e
            
        try:
            if self.params['nu_h'] < 0 or self.params['nu_h'] > 1:
                raise ValueError('nu_h must be in [0,1]')
        except ValueError as e:
            self.logger.exception('nu_h must be in [0,1]')
            raise e
        
        try:
            if self.params['mu_v'] < 0 or self.params['mu_v'] > 1:
                raise ValueError('mu_v must be in [0,1]')
        except ValueError as e:
            self.logger.exception('mu_v must be in [0,1]')
            raise e
        
        try:
            if self.params['gamma_h'] < 0 or self.params['gamma_h'] > 1:
                raise ValueError('gamma_h must be in [0,1]')
        except ValueError as e:
            self.logger.exception('gamma_h must be in [0,1]')
            raise e
        
        try:
            if self.params['beta_h'] < 0 or self.params['beta_h'] > 1:
                raise ValueError('beta_h must be in [0,1]')
        except ValueError as e:
            self.logger.exception('beta_h must be in [0,1]')
            raise e
            
        try:
            if self.params['beta_v'] < 0 or self.params['beta_v'] > 1:
                raise ValueError('beta_v must be in [0,1]')
        except ValueError as e:
            self.logger.exception('beta_v must be in [0,1]')
            raise e
