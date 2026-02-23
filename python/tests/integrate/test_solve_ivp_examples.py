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
        ybar = y[0].serialize()

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

def test_solve_ivp_Pleiades():
    """
    N-body problem with n=7 "Pleiades Problem".

    This could probably get higher accuracy with GLR actions but needs rebasing every step in MK.

    References
    ----------
        [1] Hairer, Ernst, Gerhard Wanner, and Syvert P. Nørsett. Solving ordinary differential equations I:
            Nonstiff problems. Berlin, Heidelberg: Springer Berlin Heidelberg, 1993.
        [2] http://archimede.dm.uniba.it/~testset/report/plei.pdf
    """
    
    from lielab.integrate import solve_ivp, EuclideanIVPSystem, IVPOptions
    import numpy as np

    m = np.array([1,2,3,4,5,6,7])

    def eoms(t, y, m):
        n = int(y.size/4)

        px = y[0:n]
        py = y[n:2*n]
        vx = y[2*n:3*n]
        vy = y[3*n:4*n]

        pxd = vx
        pyd = vy
        vxd = np.zeros((n,))
        vyd = np.zeros((n,))

        for ii in range(n):
            for jj in range(n):
                if ii != jj:
                    rad = np.pow((px[ii] - px[jj])**2 + (py[ii] - py[jj])**2, 1.5)
                    vxd[ii] += m[jj]*(px[jj] - px[ii])/rad
                    vyd[ii] += m[jj]*(py[jj] - py[ii])/rad
        
        return np.concatenate([pxd, pyd, vxd, vyd])


    dynamics = EuclideanIVPSystem(lambda _t, _y: eoms(_t, _y, m))
    options = IVPOptions()
    options.abstol = 1e-10
    options.reltol = 1e-10

    y0 = np.array([3.0, 3.0, -1.0, -3.0, 2.0, -2.0, 2.0, 3.0, -3.0, 2.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.75, -1.5, 0.0, 0.0, 0.0, -1.25, 1.0, 0.0, 0.0])
    path = solve_ivp(dynamics, [0.0, 3.0], y0, options)

    tol = 1e-7

    assert np.abs(path.ybar[-1, 0] - 0.3706139143970502) < tol
    assert np.abs(path.ybar[-1, 1] - 3.237284092057233) < tol
    assert np.abs(path.ybar[-1, 2] + 3.222559032418324) < tol
    assert np.abs(path.ybar[-1, 3] - 0.6597091455775310) < tol
    assert np.abs(path.ybar[-1, 4] - 0.3425581707156584) < tol
    assert np.abs(path.ybar[-1, 5] - 1.562172101400631) < tol
    assert np.abs(path.ybar[-1, 6] + 0.7003092922212495) < tol

    assert np.abs(path.ybar[-1, 7] + 3.943437585517392) < tol
    assert np.abs(path.ybar[-1, 8] + 3.271380973972550) < tol
    assert np.abs(path.ybar[-1, 9] - 5.225081843456543) < tol
    assert np.abs(path.ybar[-1, 10] + 2.590612434977470) < tol
    assert np.abs(path.ybar[-1, 11] - 1.198213693392275) < tol
    assert np.abs(path.ybar[-1, 12] + 0.2429682344935824) < tol
    assert np.abs(path.ybar[-1, 13] - 1.091449240428980) < tol

    assert np.abs(path.ybar[-1, 14] - 3.417003806314313) < tol
    assert np.abs(path.ybar[-1, 15] - 1.354584501625501) < tol
    assert np.abs(path.ybar[-1, 16] + 2.590065597810775) < tol
    assert np.abs(path.ybar[-1, 17] - 2.025053734714242) < tol
    assert np.abs(path.ybar[-1, 18] + 1.155815100160448) < tol
    assert np.abs(path.ybar[-1, 19] + 0.8072988170223021) < tol
    assert np.abs(path.ybar[-1, 20] - 0.5952396354208710) < tol

    assert np.abs(path.ybar[-1, 21] + 3.741244961234010) < tol
    assert np.abs(path.ybar[-1, 22] - 0.3773459685750630) < tol
    assert np.abs(path.ybar[-1, 23] - 0.9386858869551073) < tol
    assert np.abs(path.ybar[-1, 24] - 0.3667922227200571) < tol
    assert np.abs(path.ybar[-1, 25] + 0.3474046353808490) < tol
    assert np.abs(path.ybar[-1, 26] - 2.344915448180937) < tol
    assert np.abs(path.ybar[-1, 27] + 1.947020434263292) < tol


def test_solve_ivp_Particle_Magnetic():
    """
    Solves the motion of a charged particle in a magnetic dipole field. Ex 11.3 from Ref 1. Note that
    the authors construct the action as (GL(6), R^6) -> R^6 and do exponential coordinates on
    the entire gl/GL(6). Numerically this isn't very accurate, especially as theta drifts away from the
    0 of the algebra. So I rewrote the action as (R^3 x SO(3), R^3 x R^3) -> (R^3 x R^3) since SO(3)
    is much more well behaved. Resulting path is a particle trapped in the van Allen belt.

    References
    ----------
        [1] Iserles, Arieh, Hans Z. Munthe-Kaas, Syvert P. Nørsett, and Antonella Zanna. "Lie-group methods."
            Acta numerica 9 (2000): 215-365.
    """
    
    from lielab.domain import so, rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import solve_ivp, HomogeneousIVPSystem, IVPOptions
    import numpy as np

    tspan = [0.0, 500.0]
    y0 = CompositeManifold([RN.from_vector([0, -2.5, 0]), RN.from_vector([0, 0, 0.012])])

    def vectorfield(t, y):
        p = y[0].serialize()
        v = y[1].serialize()
        azimuth = np.arctan2(p[1], p[0])
        elevation = np.arctan2(p[2], np.sqrt(p[0]**2 + p[1]**2))
        rad = np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)

        # m = [0, 0, 1] # Magnetic dipole vector
        br = 2*np.sin(elevation)/rad**3
        bx = br*np.cos(elevation)*np.cos(azimuth)
        by = br*np.cos(elevation)*np.sin(azimuth)
        bz = br*np.sin(elevation)
        btheta = np.cos(elevation)/rad**3
        bthx = btheta*np.cos(elevation - np.pi/2)*np.cos(azimuth)
        bthy = btheta*np.cos(elevation - np.pi/2)*np.sin(azimuth)
        bthz = btheta*np.sin(elevation - np.pi/2)
        b1 = bx + bthx
        b2 = by + bthy
        b3 = bz + bthz

        return CompositeAlgebra([rn.from_vector(v), so.from_vector([b1,b2,b3])])
    
    def action(g, y):
        pnext = g[0].serialize() + y[0].serialize()
        vnext = np.dot(g[1].get_matrix(), y[1].serialize())
        return CompositeManifold([RN.from_vector(pnext), RN.from_vector(vnext)])
    
    dynamics = HomogeneousIVPSystem(vectorfield, action=action)

    options = IVPOptions()
    options.dt_max = 1.0

    path = solve_ivp(dynamics, tspan, y0, options)

    assert (path.success == True)

    assert (np.abs(path.ybar[-1, 0] - -0.527505080797981) < 1e-3)
    assert (np.abs(path.ybar[-1, 1] - -2.434436925360713) < 1e-3)
    assert (np.abs(path.ybar[-1, 2] - -0.145842164780228) < 1e-3)
    assert (np.abs(path.ybar[-1, 3] - 0.000006752678154) < 1e-3)
    assert (np.abs(path.ybar[-1, 4] - 0.001167463126786) < 1e-3)
    assert (np.abs(path.ybar[-1, 5] - -0.011943072646892) < 1e-3)

    # MK integration perfectly preserves Hamiltonian = |v|^2 since acceleration is orthogonal to velocity.
    assert (np.abs(np.linalg.norm(path.ybar[0, 3:6]) - np.linalg.norm(path.ybar[-1, 3:6])) < 1e-16)

def test_solve_ivp_Composite_Coadjoint():
    """
    Tests solve_ivp against a function with composite custom action (SE x GLR) x (SE x RN) -> (SE x RN).
    Also checks the integration of coadjoint actions and custom connections.

    """

    from lielab.domain import se, RN, CompositeAlgebra, CompositeManifold
    from lielab.functions import coad, exp, dexpinv
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPOptions
    import numpy as np

    def eoms(t, y):
        V = 1.0
        lambdabar = y[1].serialize()
        u = -lambdabar[2]

        dx = se.from_vector([V, 0.0, u])
        dlambdabar = -coad(dx)
        
        return CompositeAlgebra([dx, dlambdabar])

    def action(g, y):
        g0 = g[0]
        y0 = y[0]
        coAdyhat = g[1].get_matrix()
        lambdabar = y[1].serialize()
        lambdanext = RN.from_vector(np.dot(coAdyhat,lambdabar))
        return CompositeManifold([y0*g0, lambdanext])

    def connection(theta, xi):
        return CompositeAlgebra([dexpinv(-theta[0], xi[0]), dexpinv(theta[1], xi[1])])

    x0bar = np.array([0.0, 0.0, np.pi/2.0])
    x0 = se.from_vector(x0bar)

    lambda0bar = np.zeros((3,))
    lambda0bar[0] = 1.15407533e-03
    lambda0bar[1] = -3.17495766e+01
    lambda0bar[2] = -4.41935411e+00

    tspan = np.linspace(0.0, 1.0, 2) # TODO:

    y0 = CompositeManifold([exp(x0), RN.from_vector(lambda0bar)])

    dynamics = HomogeneousIVPSystem(eoms, action=action, connection=connection)

    options = IVPOptions()
    options.dt_max = 0.1

    curve = solve_ivp(dynamics, tspan, y0)

    assert np.abs(curve.t[0] - 0.0) < 1e-15
    assert curve.ybar[0, 0] == x0bar[0]
    assert curve.ybar[0, 1] == x0bar[1]
    assert curve.ybar[0, 3] == -1.0
    assert curve.ybar[0, 4] == 1.0
    assert curve.ybar[0, 6] == lambda0bar[0]
    assert curve.ybar[0, 7] == lambda0bar[1]
    assert curve.ybar[0, 8] == lambda0bar[2]

    assert np.abs(curve.t[-1] - 1.0) < 1e-15
    assert np.abs(curve.ybar[-1, 0] - 0.12732395447351627) < 1e-1
    assert np.abs(curve.ybar[-1, 1] + 0.0) < 1e-1
    assert np.abs(curve.ybar[-1, 6] - lambda0bar[0]) < 1e-1
    assert np.abs(curve.ybar[-1, 7] + lambda0bar[1]) < 1e-1
    assert np.abs(curve.ybar[-1, 8] - lambda0bar[2]) < 1e-1


def test_solve_ivp_AEIF():
    """
    Tests solve_ivp with the Adaptive Exponential Integrate and Fire model. Primarily checks
    robustness of event handling.

    References
    ----------
        [1] https://en.wikipedia.org/wiki/Exponential_integrate-and-fire
        [2] https://nest-simulator.readthedocs.io/en/stable/model_details/aeif_models_implementation.html
    """

    from lielab.integrate import solve_ivp, EuclideanIVPSystem, IVPOptions, IVPStatus
    import numpy as np
    import copy

    V_reset = -58.0
    V_peak = 0.0
    V_th = -50.0
    I_e = 420.0
    g_L = 11.0
    tau_w = 300.0
    E_L = -70.0
    Delta_T = 2.0
    a = 3.0
    C_m = 200.0

    def eoms(t, y):
        v = min(y[0], V_peak)
        w = y[1]

        Ispike = 0.0

        if Delta_T != 0.0:
            Ispike = g_L * Delta_T * np.exp((v - V_th) / Delta_T)

        dv = (-g_L * (v - E_L) + Ispike - w + I_e) / C_m
        dw = (a * (v - E_L) - w) / tau_w

        return np.array([dv, dw])

    def event(t, y):
        return V_peak - y[0]


    dynamics = EuclideanIVPSystem(eoms, event=event)

    options = IVPOptions()
    options.dt_max = 0.01

    tf = 100.0
    first_run = True
    tf_reached = False
    n_events = 0

    while not tf_reached:
        if first_run:
            tspan = np.array([0.0, tf])
            y0 = np.array([E_L, 5.0])
        else:
            tspan = np.array([sol.t[-1], tf])
            y0 = np.array([V_reset, sol.ybar[-1, 1]])
        
        seg = solve_ivp(dynamics, tspan, y0, options)

        if (seg.status == IVPStatus.SUCCESS_EVENT or seg.status == IVPStatus.SUCCESS_EVENT_BUT_TOL):
            assert (np.abs(seg.ybar[-1, 0] - V_peak) < 0.1)
            n_events += 1

        if first_run:
            sol = copy.copy(seg)
        else:
            sol.t = np.concatenate((sol.t, seg.t))
            sol.ybar = np.vstack((sol.ybar, seg.ybar))
            sol.thetabar = np.vstack((sol.thetabar, seg.thetabar))
            sol.y = sol.y + seg.y
            sol.theta = sol.theta + seg.theta

        if sol.t[-1] >= tf:
            tf_reached = True

        first_run = False

    assert (n_events == 7)
    assert (np.abs(sol.t[-1] - tf) < 1e-15)
