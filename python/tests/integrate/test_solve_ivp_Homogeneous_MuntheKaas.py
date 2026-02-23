"""
Unit-like tests for solve_ivp with Homogeneous systems.
"""

import pytest

"""
Error handling
"""

def test_solve_ivp_Homogeneous_MuntheKaas_nan():
    """
    Tests solve_ivp for error handling when nans appear in the vectorfield.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPStatus, IVPOptions, IVPMethod
    import numpy as np

    m = -0.5
    tf = 10.0

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]*np.nan
        return CompositeAlgebra([rn.from_vector(dybar)])

    tspan = [0.0, tf]

    dynamics = HomogeneousIVPSystem(eoms)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert (curve.status == IVPStatus.ERROR_NANS_IN_VF)
    assert (curve.success == False)

def test_solve_ivp_Homogeneous_MuntheKaas_inf():
    """
    Tests solve_ivp for error handling when infs appear in the vectorfield.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPStatus, IVPOptions, IVPMethod
    import numpy as np

    m = -0.5
    tf = 10.0

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]*np.inf
        return CompositeAlgebra([rn.from_vector(dybar)])

    tspan = [0.0, tf]

    dynamics = HomogeneousIVPSystem(eoms)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert (curve.status == IVPStatus.ERROR_INFS_IN_VF)
    assert (curve.success == False)

def test_solve_ivp_Homogeneous_MuntheKaas_tspan_short():
    """
    Tests solve_ivp for error handling when tspan is too short.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPStatus, IVPOptions, IVPMethod
    import numpy as np

    m = -0.5
    tf = 10.0

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]
        return CompositeAlgebra([rn.from_vector(dybar)])

    tspan = [0.0]

    dynamics = HomogeneousIVPSystem(eoms)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    with pytest.raises(RuntimeError):
        curve = solve_ivp(dynamics, tspan, y0, options)

def test_solve_ivp_Homogeneous_MuntheKaas_tspan_repeat():
    """
    Tests solve_ivp for error handling when tspan has repeated values.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPStatus, IVPOptions, IVPMethod
    import numpy as np

    m = -0.5
    tf = 10.0

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]
        return CompositeAlgebra([rn.from_vector(dybar)])

    tspan = [0.0, 0.0]

    dynamics = HomogeneousIVPSystem(eoms)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    with pytest.raises(RuntimeError):
        curve = solve_ivp(dynamics, tspan, y0, options)

def test_solve_ivp_Homogeneous_MuntheKaas_tspan_descending():
    """
    Tests solve_ivp for error handling when tspan is in descending order.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPStatus, IVPOptions, IVPMethod
    import numpy as np

    m = -0.5
    tf = 10.0

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]
        return CompositeAlgebra([rn.from_vector(dybar)])

    tspan = [0.0, -tf]

    dynamics = HomogeneousIVPSystem(eoms)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    with pytest.raises(RuntimeError):
        curve = solve_ivp(dynamics, tspan, y0, options)

def test_solve_ivp_Homogeneous_MuntheKaas_topos():
    """
    Tests solve_ivp for error handling when user provided topology comes out incorrect.
    """

    from lielab.domain import so, SO, SU, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPOptions, IVPMethod, IVPStatus
    import numpy as np

    def eoms(t, y):
        return CompositeAlgebra([so.from_vector([1.0, 2.0, 3.0])])

    def action(g, y):
        return CompositeManifold([g[0]*y[0], SU(2)])

    tspan = np.linspace(0.0, 5.0, 2)

    dynamics = HomogeneousIVPSystem(eoms, action=action)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    with pytest.raises(RuntimeError):
        curve = solve_ivp(dynamics, tspan, CompositeManifold([SO(3)]), options)

"""
Basic features
"""

def test_solve_ivp_Homogeneous_MuntheKaas():
    """
    Tests solve_ivp against a classical problem with known solution.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPStatus, IVPOptions, IVPMethod
    import numpy as np

    m = -0.5
    tf = 10.0

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]
        return CompositeAlgebra([rn.from_vector(dybar)])

    tspan = [0.0, tf]

    dynamics = HomogeneousIVPSystem(eoms)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert (curve.status == IVPStatus.SUCCESS)
    assert (curve.success == True)

    assert (np.abs(curve.t[0] - 0.0) < 1e-15)
    assert (curve.ybar[0, 0] == y0bar[0])
    assert (curve.ybar[0, 1] == y0bar[1])
    assert (curve.ybar[0, 2] == y0bar[2])
    assert (curve.ybar[0, 3] == y0bar[3])


    assert (np.abs(curve.t[-1] - tf) < 1e-15)
    assert (np.abs(curve.ybar[-1, 0] - np.exp(m*tf)*y0bar[0]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 1] - np.exp(m*tf)*y0bar[1]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 2] - np.exp(m*tf)*y0bar[2]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 3] - np.exp(m*tf)*y0bar[3]) < 1e-5)

def test_solve_ivp_Homogeneous_MuntheKaas_tol():
    """
    Tests solve_ivp for reporting tolerance issues when the step cannot adapt properly.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPOptions, IVPStatus, IVPOptions, IVPMethod
    import numpy as np

    m = -0.5
    tf = 10.0

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]
        return CompositeAlgebra([rn.from_vector(dybar)])

    tspan = [0.0, tf]

    dynamics = HomogeneousIVPSystem(eoms)

    options = IVPOptions()
    options.dt_min = 5.0
    options.method = IVPMethod.MuntheKaas

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert (curve.status == IVPStatus.SUCCESS_BUT_TOL)
    assert (curve.success == True)

    assert (curve.t.size == 3)
    # assert (len(curve.y) == 3)
    assert (curve.ybar.shape[0] == 3)
    assert (curve.ybar.shape[1] == 4)
    # assert (len(curve.theta) == 3)
    assert (curve.thetabar.shape[0] == 3)
    assert (curve.thetabar.shape[1] == 4)

def test_solve_ivp_Homogeneous_MuntheKaas_event():
    """
    Tests solve_ivp handling of events with a classical problem with known solution.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPOptions, IVPStatus, IVPOptions, IVPMethod
    import numpy as np

    m = -0.5
    tf = 10.0
    y2cross = 3.5

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]
        return CompositeAlgebra([rn.from_vector(dybar)])

    def event(t, y, y2cross=y2cross):
        ybar = y.serialize()
        return ybar[2] - y2cross

    tspan = [0.0, tf]

    dynamics = HomogeneousIVPSystem(eoms, event=event)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert (curve.status == IVPStatus.SUCCESS_EVENT)
    assert (curve.success == True)

    assert (np.abs(curve.t[0] - 0.0) < 1e-15)
    assert (curve.ybar[0, 0] == y0bar[0])
    assert (curve.ybar[0, 1] == y0bar[1])
    assert (curve.ybar[0, 2] == y0bar[2])
    assert (curve.ybar[0, 3] == y0bar[3])

    tcross = np.log(y2cross/y0bar[2])/m

    assert (np.abs(curve.t[-1] - tcross) < options.reltol*10.0)
    assert (np.abs(curve.ybar[-1, 0] - np.exp(m*tcross)*y0bar[0]) < options.reltol*10.0)
    assert (np.abs(curve.ybar[-1, 1] - np.exp(m*tcross)*y0bar[1]) < options.reltol*10.0)
    assert (np.abs(curve.ybar[-1, 2] - np.exp(m*tcross)*y0bar[2]) < options.reltol*10.0)
    assert (np.abs(curve.ybar[-1, 3] - np.exp(m*tcross)*y0bar[3]) < options.reltol*10.0)

def test_solve_ivp_Homogeneous_MuntheKaas_event_tol():
    """
    Tests solve_ivp handling of events and tolerance issues with a classical problem with known solution.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPOptions, IVPStatus, IVPMethod
    import numpy as np

    m = -0.5
    tf = 10.0
    y2cross = 0.25

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]
        return CompositeAlgebra([rn.from_vector(dybar)])
    
    def event(t, y, y2cross=y2cross):
        ybar = y.serialize()
        return ybar[2] - y2cross

    tspan = [0.0, tf]

    dynamics = HomogeneousIVPSystem(eoms, event=event)

    options = IVPOptions()
    options.dt_min = 5.0
    options.method = IVPMethod.MuntheKaas

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert (curve.status == IVPStatus.SUCCESS_EVENT_BUT_TOL)
    assert (curve.success == True)

    assert (curve.t.size == 3)
    # assert (len(curve.y) == 3)
    assert (curve.ybar.shape[0] == 3)
    assert (curve.ybar.shape[1] == 4)
    # assert (len(curve.theta) == 3)
    assert (curve.thetabar.shape[0] == 3)
    assert (curve.thetabar.shape[1] == 4)

def test_solve_ivp_Homogeneous_MuntheKaas_segmented():
    """
    Tests solve_ivp exactly reporting elements in tspan on a classical problem with known solution.
    """

    from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPStatus, IVPOptions, IVPMethod
    import numpy as np

    m = -0.5
    tspan = np.linspace(0,10,41)

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]
        return CompositeAlgebra([rn.from_vector(dybar)])

    dynamics = HomogeneousIVPSystem(eoms)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert (curve.status == IVPStatus.SUCCESS)
    assert (curve.success == True)

    for ii in range(len(tspan)):
        assert tspan[ii] in curve.t

"""
Advanced features
"""

def test_solve_ivp_Homogeneous_MuntheKaas_Lorenz_GLRxRN_RN():
    """
    Tests solve_ivp against a function with custom action GLR x RN -> RN.
    """

    from lielab.domain import glr, RN, CompositeAlgebra, CompositeManifold
    from lielab.integrate import HomogeneousIVPSystem, solve_ivp, IVPOptions, IVPMethod
    import numpy as np

    def eoms(t, y):
        ybar = y[0].serialize()

        sigma = 10.0
        rho = 28.0
        b1 = 8.0
        b2 = 3.0
        beta = b1/b2

        A = np.zeros((3,3))

        A[0, 0] = -beta
        A[0, 1] = 0.0
        A[0, 2] = ybar[1]
        A[1, 0] = 0.0
        A[1, 1] = -sigma
        A[1, 2] = sigma
        A[2, 0] = -ybar[1]
        A[2, 1] = rho
        A[2, 2] = -1.0

        return CompositeAlgebra([glr(A)])

    def action(g, y):
        ghat = g[0].get_matrix()
        ybar = y[0].serialize()
        return CompositeManifold([RN.from_vector(np.dot(ghat, ybar))])

    y0bar = np.zeros((3,))
    y0bar[0] = 25.0
    y0bar[1] = 0.0
    y0bar[2] = -20.0

    tspan = np.linspace(0.0, 5.0, 20)

    dynamics = HomogeneousIVPSystem(eoms, action=action)

    options = IVPOptions()
    options.method = IVPMethod.MuntheKaas

    curve = solve_ivp(dynamics, tspan, CompositeManifold([RN.from_vector(y0bar)]), options)

    assert np.abs(curve.t[0] - 0.0) < 1e-15
    assert curve.ybar[0, 0] == y0bar[0]
    assert curve.ybar[0, 1] == y0bar[1]
    assert curve.ybar[0, 2] == y0bar[2]

    assert np.abs(curve.t[-1] - 5.0) < 1e-15
    assert np.abs(curve.ybar[-1, 0] - 15.230) < 1e-1
    assert np.abs(curve.ybar[-1, 1] + 0.797) < 1e-1
    assert np.abs(curve.ybar[-1, 2] + 1.473) < 1e-1


