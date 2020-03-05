import numpy as np
import analysis.functions as f

def test_mdot():
    tol = 1e-2

    mdot_clc = 9.48686
    mdot_fcn = f.mdot(1e5, 0.1, 0.6, 400, 1000)

    assert abs(mdot_fcn - mdot_clc) < tol
    
    a = np.array([1,2,3,4,5])
    mdot_clc_vec = a * mdot_clc
    p = a * 1e5

    mdot_fcn_vec = f.mdot(p, 0.1, 0.6, 400, 1000)

    assert np.allclose(mdot_fcn_vec, mdot_clc_vec, tol)

