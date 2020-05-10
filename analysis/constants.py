import numpy as np
from dataclasses import dataclass

from .functions import vdkh

# Propellant properties - water
@dataclass
class Water:
    h = 2256e3      #J/kg
    c_l = 4187      #J/K/kg liquid water
    c_v = 1996      #J/K/kg water vapour
    gamma = 1.33
    R_constant = 461.67      #J/K/kg
    rho = 997       #kg/m^3
    p_vap0 = 1e5
    T_vap0 = 373
    h_vap = 40e3
    R_vap = 8.341
    Gamma = vdkh(1.33)
