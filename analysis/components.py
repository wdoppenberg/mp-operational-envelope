import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dataclasses import dataclass

class OperationalEnvelope:
    def __init__(self, kwargs):
        self.A_t = kwargs['A_t'] #Nozzle throat area
        self.l = kwargs['l'] # Propellant tubing length
        self.R_constant = kwargs['R_constant'] # Gas constant water vapour
        self.rho = kwargs['rho'] # Density propellant
        self.gamma = kwargs['gamma']
        self.d = kwargs['d'] # Inner diameter propellant tubing
        self.T_0 = kwargs['T_0'] # Ambient temperature in satellite
        self.h_vap = kwargs['h_vap'] # Heat of vapourization water vapour
        self.c_pl = kwargs['c_pl'] # Specific heat of liquid
        self.c_pv = kwargs['c_pv'] # Specific heat of vapour

        self.df = pd.DataFrame(columns=['Symbol', 'Value'])
        for i in kwargs.items():
            self.df = self.df.append({'Symbol': i[0], 'Value':i[1]}, ignore_index=True)


    @property
    def Gamma(self):
        g = self.gamma
        return np.sqrt(g) * ((2/(g+1)) ** ((g+1)/(2*(g-1))))

    def __repr__(self):
        return repr(self.df)