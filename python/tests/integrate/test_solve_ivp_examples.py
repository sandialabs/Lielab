"""
Tests solve_ivp can be used with "real" problems.
"""

def test_solve_ivp_Lorenz_RNxRN_RN():
    """
    Tests solve_ivp with the Lorenz equations.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPOptions
    import numpy as np

    def eoms(t, y):
        ybar = y.space[0].serialize()

        dybar = np.zeros((3,))
        sigma = 10.0
        rho = 28.0
        b1 = 8.0
        b2 = 3.0
        beta = b1/b2

        dybar[0] = -beta*ybar[0] + ybar[1]*ybar[2]
        dybar[1] = -sigma*ybar[1] + sigma*ybar[2]
        dybar[2] = -ybar[0]*ybar[1] + rho*ybar[1] - ybar[2]
        return CompositeAlgebra([rn.from_vector(dybar)])

    y0bar = np.zeros((3,))
    y0bar[0] = 25.0
    y0bar[1] = 0.0
    y0bar[2] = -20.0

    tspan = [0.0, 5.0]

    dynamics = HomogeneousIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, CompositeManifold([RN.from_vector(y0bar)]))

    assert np.abs(curve.t[0] - 0.0) < 1e-15
    assert curve.ybar[0, 0] == y0bar[0]
    assert curve.ybar[0, 1] == y0bar[1]
    assert curve.ybar[0, 2] == y0bar[2]

    assert np.abs(curve.t[-1] - 5.0) < 1e-15
    assert np.abs(curve.ybar[-1, 0] - 15.230) < 1e-1
    assert np.abs(curve.ybar[-1, 1] + 0.797) < 1e-1
    assert np.abs(curve.ybar[-1, 2] + 1.473) < 1e-1

def test_solve_ivp_Composite_Coadjoint():
    """
    Tests solve_ivp against a function with composite custom action (SE x GLR) x (SE x RN) -> (SE x RN).
    Also checks the integration of coadjoint actions.

    TODO
    ----
        - Rebasing every step will significantly improve accuracy
    """

    from lielab.domain import se, RN, CompositeAlgebra, CompositeManifold
    from lielab.functions import coad, exp
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPOptions
    import numpy as np

    def eoms(t, y):
        V = 1.0
        lambdabar = y.space[1].serialize()
        u = -lambdabar[2]

        dx = se.from_vector([V, 0.0, u])
        dlambdabar = -coad(dx)
        
        return CompositeAlgebra([dx, dlambdabar])

    def action(g, y):
        g0 = g.space[0]
        y0 = y.space[0]
        coAdyhat = g.space[1].get_matrix()
        lambdabar = y.space[1].serialize()
        lambdanext = RN.from_vector(np.dot(coAdyhat,lambdabar))
        return CompositeManifold([y0*g0, lambdanext])

    x0bar = np.array([0.0, 0.0, np.pi/2.0])
    x0 = se.from_vector(x0bar)

    lambda0bar = np.zeros((3,))
    lambda0bar[0] = 1.15407533e-03
    lambda0bar[1] = -3.17495766e+01
    lambda0bar[2] = -4.41935411e+00

    tspan = np.linspace(0.0, 1.0, 100) # TODO:

    y0 = CompositeManifold([exp(x0), RN.from_vector(lambda0bar)])

    dynamics = HomogeneousIVPSystem(eoms, action=action)

    curve = solve_ivp(dynamics, tspan, y0)

    assert np.abs(curve.t[0] - 0.0) < 1e-15
    assert curve.ybar[0, 1] == -1.0
    assert curve.ybar[0, 3] == 1.0
    assert curve.ybar[0, 2] == x0bar[0]
    assert curve.ybar[0, 5] == x0bar[1]
    assert curve.ybar[0, 9] == lambda0bar[0]
    assert curve.ybar[0, 10] == lambda0bar[1]
    assert curve.ybar[0, 11] == lambda0bar[2]

    assert np.abs(curve.t[-1] - 1.0) < 1e-15
    # assert np.abs(curve.ybar[-1, 2] - 0.12732395447351627) < 1e-1 # TODO:
    # assert np.abs(curve.ybar[-1, 5] + 0.0) < 1e-1
    # assert np.abs(curve.ybar[-1, 9] - lambda0bar[0]) < 1e-1
    # assert np.abs(curve.ybar[-1, 10] + lambda0bar[1]) < 1e-1
    # assert np.abs(curve.ybar[-1, 11] - lambda0bar[2]) < 1e-1


