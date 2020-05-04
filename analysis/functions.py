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

def CA_2D(theta_nd):
    """
    2D conical nozzle
    """
    return np.sin(theta_nd)/theta_nd


def CA_3D(theta_nd):
    """
    3D conical nozzle
    """
    return (1+np.cos(theta_nd))/2


def C_F_divloss(CA, C_F0, p_e, p_a, p_c, AR):
    """
    Divergence loss
    """
    return CA*C_F0 + (p_e/p_a - p_a/p_c)*AR


def C_F_visc(AR, Re):
    """
    Viscous loss
    """
    return (17.6 * np.exp(0.0032*AR))/np.sqrt(0.773*Re)


def R_mod(Re, R_curv, R_t):
    """
    Modified Reynolds number
    """
    return Re*np.sqrt(R_curv/R_t)


def C_d_bdloss(g, Re, R_curv, R_t):
    """
    Throat boundary layer loss
    """
    A = ((g+1)/2)**(3/4) * ( (72-32*np.sqrt(6))/(3*(g+1)) + (4*np.sqrt(6))/3 )\
            * (1/np.sqrt(R_mod(Re, R_curv, R_t)))
    B = ((2*np.sqrt(2)*(g-1)*(g+2))/(3*np.sqrt(g+1)))*(1/R_mod(Re, R_curv, R_t))

    return 1 - A + B


def dp_pipe(f, L, D, rho, v):
    """
    Pressure loss pipe
    """
    return f*(L/(2*D))*rho*v**2