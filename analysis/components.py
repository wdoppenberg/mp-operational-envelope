import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dataclasses import dataclass

@dataclass
class Thruster:
    A_t: float

@dataclass
class Propellant:
    h: float
    c_l: float
    c_v: float
    R: float
    rho: float
    Gamma: float

@dataclass
class Environment:
    T_0: float
    l: float
    d: float

class OperationalEnvelope:
    def __init__(
            self, 
            thruster:       Thruster,
            propellant:     Propellant,
            environment:    Environment
            ):
            
        self.th = thruster
        self.prop = propellant
        self.env = environment
