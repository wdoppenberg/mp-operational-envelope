import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dataclasses import dataclass
from .functions import vdkh

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

    def describe(self):
        print(f'''
        Nitrogen pressure \t[{min(self.p):.2e}, {max(self.p):.2e}] Pa
        Mass flow \t\t[{min(self.mdot):.2e}, {max(self.mdot):.2e}] kg/s
        Thrust \t\t\t[{min(self.F_T):.2e}, {max(self.F_T):.2e}] N
        Chamber temperature \t[{min(self.T_c):.2e}, {max(self.T_c):.2e}] K
        Propellant mass \t[{min(self.m):.2e}, {max(self.m):.2e}] kg
        Required power \t\t[{min(self.Q):.2e}, {max(self.Q):.2e}] W
        Input power \t\t[{min(self.Q/self.eff_Q):.2e}, {max(self.Q/self.eff_Q):.2e}] W
        ''')

    def plot(self, return_fig=False):
        plt.style.use('ggplot')

        fig, ax = plt.subplots(6, 1, figsize=(10,10), sharex=True)
        ax[0].plot(self.t, self.p)
        ax[0].set_ylabel('Nitrogen\npressure [Pa]')
        ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        ax[1].plot(self.t, self.mdot)
        ax[1].set_ylabel('Mass\nflow [kg/s]')
        ax[1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        ax[2].plot(self.t, self.F_T)
        ax[2].set_ylabel('Thrust [N]')
        ax[2].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        ax[3].plot(self.t, self.T_c)
        ax[3].set_ylabel('Chamber\ntemperature [K]')
        ax[3].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        ax[4].plot(self.t, self.m)
        ax[4].set_ylabel('Propellant\nmass [kg]')
        ax[4].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        ax[5].plot(self.t, self.Q)
        ax[5].set_ylabel('Required\npower [W]')
        ax[5].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        ax[5].set_xlabel('Time [s]')

        if return_fig:
            return fig
        

    def __repr__(self):
        s = 'OperationalEnvelope(\n'
        for k, v in vars(self).items():
            s+=f"{k}\t\t\t{v}\n"
        s += ')'
        return s
