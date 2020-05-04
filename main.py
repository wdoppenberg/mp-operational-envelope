from analysis import constants as c, functions as f, components as cmp

input_vars = dict(
    A_t=4.5e-9,
    l=0.3,
    R_constant=461.67,
    rho=997,
    gamma=1.3,
    d=1.57e-3,
    T_0=283,
    h_vap=2256,
    c_pl=4187,
    c_pv=1996
)

oe = cmp.OperationalEnvelope(input_vars)

