"""
Unit-like tests for solve_ivp with Euclidean systems.
"""

"""
Error handling
"""

def test_solve_ivp_Euclidean_nan():
    """
    Tests solve_ivp for error handling when nans appear in the vectorfield.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp
    import numpy as np

    m = -0.5
    tf = 10.0

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]*np.nan
        return dy

    tspan = [0.0, tf]

    dynamics = EuclideanIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, y0)

    assert (curve.status == -1)
    assert (curve.success == False)

def test_solve_ivp_Euclidean_inf():
    """
    Tests solve_ivp for error handling when infs appear in the vectorfield.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp
    import numpy as np

    m = -0.5
    tf = 10.0

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]*np.inf
        return dy

    tspan = [0.0, tf]

    dynamics = EuclideanIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, y0)

    assert (curve.status == -2)
    assert (curve.success == False)

def test_solve_ivp_Euclidean_tspan_short():
    """
    Tests solve_ivp for error handling when tspan is too short.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp
    import numpy as np

    m = -0.5
    tf = 10.0

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]
        return dy

    tspan = [0.0]

    dynamics = EuclideanIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, y0)

    assert (curve.status == -3)
    assert (curve.success == False)

def test_solve_ivp_Euclidean_tspan_repeat():
    """
    Tests solve_ivp for error handling when tspan has repeated values.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp
    import numpy as np

    m = -0.5
    tf = 10.0

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]
        return dy

    tspan = [0.0, 0.0]

    dynamics = EuclideanIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, y0)

    assert (curve.status == -20)
    assert (curve.success == False)

def test_solve_ivp_Euclidean_tspan_descending():
    """
    Tests solve_ivp for error handling when tspan is in descending order.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp
    import numpy as np

    m = -0.5
    tf = 10.0

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]
        return dy

    tspan = [0.0, -tf]

    dynamics = EuclideanIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, y0)

    assert (curve.status == -20)
    assert (curve.success == False)

"""
Basic features
"""

def test_solve_ivp_Euclidean():
    """
    Tests solve_ivp against a classical problem with known solution.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp
    import numpy as np

    m = -0.5
    tf = 10.0

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]
        return dy

    tspan = [0.0, tf]

    dynamics = EuclideanIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, y0)

    assert (np.abs(curve.t[0] - 0.0) < 1e-15)
    assert (curve.ybar[0, 0] == y0[0])
    assert (curve.ybar[0, 1] == y0[1])
    assert (curve.ybar[0, 2] == y0[2])
    assert (curve.ybar[0, 3] == y0[3])


    assert (np.abs(curve.t[-1] - tf) < 1e-15)
    assert (np.abs(curve.ybar[-1, 0] - np.exp(m*tf)*y0[0]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 1] - np.exp(m*tf)*y0[1]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 2] - np.exp(m*tf)*y0[2]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 3] - np.exp(m*tf)*y0[3]) < 1e-5)

def test_solve_ivp_Euclidean_tol():
    """
    Tests solve_ivp for reporting tolerance issues when the step cannot adapt properly.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp, IVPOptions
    import numpy as np

    m = -0.5
    tf = 10.0

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]
        return dy

    tspan = [0.0, tf]

    dynamics = EuclideanIVPSystem(eoms)

    options = IVPOptions()
    options.dt_min = 5.0

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert (curve.status == 4)
    assert (curve.success == True)

    assert (curve.t.size == 3)
    # assert (len(curve.y) == 3)
    assert (curve.ybar.shape[0] == 3)
    assert (curve.ybar.shape[1] == 4)
    # assert (len(curve.theta) == 3)
    assert (curve.thetabar.shape[0] == 3)
    assert (curve.thetabar.shape[1] == 4)

def test_solve_ivp_Euclidean_event():
    """
    Tests solve_ivp handling of events with a classical problem with known solution.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp
    import numpy as np

    m = -0.5
    tf = 10.0
    y2cross = 3.5

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]
        return dy

    def event(t, y, y2cross=y2cross):
        return y[2] - y2cross

    tspan = [0.0, tf]

    dynamics = EuclideanIVPSystem(eoms, event=event)

    curve = solve_ivp(dynamics, tspan, y0)

    assert (np.abs(curve.t[0] - 0.0) < 1e-15)
    assert (curve.ybar[0, 0] == y0[0])
    assert (curve.ybar[0, 1] == y0[1])
    assert (curve.ybar[0, 2] == y0[2])
    assert (curve.ybar[0, 3] == y0[3])

    tcross = np.log(y2cross/y0[2])/m

    assert (np.abs(curve.t[-1] - tcross) < 1e-8)
    assert (np.abs(curve.ybar[-1, 0] - np.exp(m*tcross)*y0[0]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 1] - np.exp(m*tcross)*y0[1]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 2] - np.exp(m*tcross)*y0[2]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 3] - np.exp(m*tcross)*y0[3]) < 1e-5)

def test_solve_ivp_Euclidean_event_tol():
    """
    Tests solve_ivp handling of events and tolerance issues with a classical problem with known solution.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp, IVPOptions
    import numpy as np

    m = -0.5
    tf = 10.0
    y2cross = 0.25

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]
        return dy
    
    def event(t, y, y2cross=y2cross):
        return y[2] - y2cross

    tspan = [0.0, tf]

    dynamics = EuclideanIVPSystem(eoms, event=event)

    options = IVPOptions()
    options.dt_min = 5.0

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert (curve.status == 5)
    assert (curve.success == True)

    assert (curve.t.size == 3)
    # assert (len(curve.y) == 3)
    assert (curve.ybar.shape[0] == 3)
    assert (curve.ybar.shape[1] == 4)
    # assert (len(curve.theta) == 3)
    assert (curve.thetabar.shape[0] == 3)
    assert (curve.thetabar.shape[1] == 4)

def test_solve_ivp_Euclidean_segmented():
    """
    Tests solve_ivp exactly reporting elements in tspan on a classical problem with known solution.
    """

    from lielab.integrate import EuclideanIVPSystem, solve_ivp
    import numpy as np

    m = -0.5
    tspan = np.linspace(0,10,41)

    y0 = np.array([2.0, 4.0, 6.0, -8.0])

    def eoms(t, y, m=m):
        dy = np.zeros((4,))

        dy[0] = m*y[0]
        dy[1] = m*y[1]
        dy[2] = m*y[2]
        dy[3] = m*y[3]
        return dy

    dynamics = EuclideanIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, y0)

    for ii in range(len(tspan)):
        assert tspan[ii] in curve.t

"""
Advanced features
"""
