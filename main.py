from analysis import constants as c, functions as f, components as comp

prop = comp.Propellant(
    c.h,
    c.c_l,
    c.c_v,
    c.R,
    c.rho,
    c.Gamma
)

th = comp.Thruster(
    c.A_t
)

env = comp.Environment(
    c.T_0,
    c.l,
    c.d
)

op = comp.OperationalEnvelope(
    th,
    prop,
    env
)

print(op)