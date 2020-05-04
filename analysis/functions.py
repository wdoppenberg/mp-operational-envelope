import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def mdot(p, A_t, Gamma, R, T_c):
    """
    Mass flow function.

    Accepts type `np.array` for `p` and `T_c`
    """
    
    return (p*A_t*Gamma) / np.sqrt(R*T_c)


def p_t(V_0, p_0, m_exit_t, rho):
    """
    Pressure vs time function.
    """
    return (V_0 * p_0) / (V_0 + (m_exit_t/rho))


def vdkh(g):
    """
    Vanderkerckhoven function
    """
    return np.sqrt(g) * ((2/(g+1)) ** ((g+1)/(2*(g-1))))