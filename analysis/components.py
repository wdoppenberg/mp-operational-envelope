import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dataclasses import dataclass
from .functions import vdkh

ex = lambda x, p: x in p

class OperationalEnvelope:
    def __init__(self, params):

        # Thruster parameters
        self.A_t = params['A_t'] #Nozzle throat area
        self.l_tube = params['l_tube'] # Propellant tubing length
        self.d_tube = params['d_tube'] # Propellant tubing diameter
        self.V_tube = (self.l_tube*np.pi*(self.d_tube**2))/4 # Propellant tubing volume
        self.eff_Q = params['eff_Q'] # Heating efficiency 
        
        # Propellant properties
        self.p = params['propellant']

        # Input parameters
        self.p_0 = params['p_0']
        self.V_fraction = params['V_fraction']
        
        if 'Tc_0' in params:
            self.Tc_0 = params['Tc_0']
        else:
            self.T_c0 = self.p.h_vap*self.p.T_vap0\
                    /(self.p.T_vap0*self.p.R_vap\
                    *np.log(self.p.p_vap0/self.p_0)\
                    +self.p.h_vap)

        # Quality factors
        self.C_d = params['C_d'] # Discharge coefficient
        self.xi_s = params['xi_s'] # I_sp quality        

    @property
    def Gamma(self):
        return vdkh(self.p.gamma)

    def __repr__(self):
        s = ''
        for k, v in vars(self).items():
            s+=f"{k}\t\t\t{v}\n"
        return s
