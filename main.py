from analysis import constants as c, functions as f, components as cmp

input_vars = dict(
    A_t=4.5e-9,
    l_tube=0.3,
    d_tube=1.57e-3,
    T_0 = 283,
    I_sp = 74,
    eff_Q=0.6,
    propellant=c.Water,
    p_0=1.1e5,
    V_fraction=0.12,
    C_d=0.66,
    xi_s=0.2
)

oe = cmp.OperationalEnvelope(input_vars)

