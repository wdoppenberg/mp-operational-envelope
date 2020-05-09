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
        self.T_0 = params['T_0']
        
        # Propellant properties
        self.prop = params['propellant']

        # Quality factors
        self.C_d = params['C_d'] # Discharge coefficient
        self.xi_s = params['xi_s'] # I_sp quality  
        self.I_sp = params['I_sp']*self.xi_s # Specific impulse

        # Input parameters
        self.p_0 = params['p_0']
        self.V_fraction = params['V_fraction']
        self.V_0 = self.V_tube*self.V_fraction
        
        if 'Tc_0' in params:
            self.Tc_0 = params['Tc_0']
        else:
            self.T_c0 = self.prop.h_vap*self.prop.T_vap0\
                    /(self.prop.T_vap0*self.prop.R_vap\
                    *np.log(self.prop.p_vap0/self.p_0)\
                    +self.prop.h_vap)

        self.mdot_0 = self.p_0*self.A_t*self.prop.Gamma\
            /(np.sqrt(self.prop.R_constant*self.T_c0))*self.C_d

    def simulate(self, dt, t_end, dim=None):
        self.t = np.arange(0, t_end, dt, dtype='f')
        if dim is None:
            dim = len(self.t)

        self.mdot = np.zeros(dim)
        self.p = np.zeros(dim)
        self.T_c = np.zeros(dim)

        self.mdot[0] = self.mdot_0
        self.p[0] = self.p_0
        self.T_c[0] = self.T_c0
        temp = 0

        for ii, _ in enumerate(self.t[1:], 1):
            self.p[ii] = self.V_0*self.p_0/(self.V_0+temp)
            self.mdot[ii] = ((self.p[ii-1]*self.A_t*self.prop.Gamma)\
                /np.sqrt(self.prop.R_constant*self.T_c[ii-1]))*self.C_d
            self.T_c[ii] = self.prop.h_vap*self.prop.T_vap0\
                /(self.prop.T_vap0*self.prop.R_vap*np.log(self.prop.p_vap0\
                /self.p[ii-1])+self.prop.h_vap)
            temp = temp+self.mdot[ii-1]*dt/self.prop.rho

        self.T_vap = self.prop.h_vap*self.prop.T_vap0\
            /(self.prop.T_vap0*self.prop.R_vap*np.log(self.prop.p_vap0/self.p)+self.prop.h_vap)
        self.Q = self.mdot * (self.T_c-self.T_0)*self.prop.c_l + self.prop.h*self.mdot
        self.V_t = self.V_0 * (self.p_0 / self.p)
        self.m = (self.V_tube - self.V_t) * self.prop.rho

        self.F_T = self.mdot*self.I_sp*9.81


    def __repr__(self):
        s = 'OperationalEnvelope(\n'
        for k, v in vars(self).items():
            s+=f"{k}\t\t\t{v:.4e}\n"
        s += ')'
        return s
