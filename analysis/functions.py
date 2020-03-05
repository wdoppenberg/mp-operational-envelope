import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def mdot(p, A_t, Gamma, R, T_c):
    """
    Mass flow function.

    Accepts type `np.array` for `p` and `T_c`
    """
    return np.divide(p*A_t*Gamma, np.sqrt(R*T_c))

