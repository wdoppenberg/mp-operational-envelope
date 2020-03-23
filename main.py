from analysis import constants as c, functions as f, components as cmp

prop = cmp.Propellant(
    c.h,
    c.c_l,
    c.c_v,
    c.R,
    c.rho,
    c.Gamma
)

th = cmp.Thruster(
    c.A_t
)

env = cmp.Environment(
    c.T_0,
    c.l,
    c.d
)

op = cmp.OperationalEnvelope(
    th,
    prop,
    env
)

if __name__ == "__main__":
    pass