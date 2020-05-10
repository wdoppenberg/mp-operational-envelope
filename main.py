from analysis import constants as c, OperationalEnvelope, Experiment
import matplotlib.pyplot as plt

input_vars1 = dict(
    l_tube=0.4,
    d_tube=1.57e-3,
    propellant=c.Water,
    p_0=1.1e5,
    V_fraction=0.12,
    C_d=0.66,
    xi_s=0.3
)

input_vars2 = dict(
    l_tube=0.3,
    d_tube=1.57e-3,
    propellant=c.Water,
    p_0=2e5,
    V_fraction=0.12,
    C_d=0.66,
    xi_s=0.3
)

input_vars3 = dict(
    l_tube=0.3,
    d_tube=1.57e-3,
    propellant=c.Water,
    p_0=1.1e5,
    V_fraction=0.2,
    C_d=0.66,
    xi_s=0.6
)

oe1 = OperationalEnvelope(input_vars1, name='oe1')
oe2 = OperationalEnvelope(input_vars2, name='oe2')
oe3 = OperationalEnvelope(input_vars3, name='oe3')

dt, t_end = 0.1, 1200
ex = Experiment(oe1, oe2, oe3)

ex.comparison(dt, t_end)

plt.show()