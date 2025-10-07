import numpy as np

from lielab.domain import *
from lielab.functions import *
from lielab.integrate import *
from lielab.testing import *

# Errors

def test_solve_ivp_nan():
    """
    """

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

    curve = solve_ivp(dynamics, tspan, y0)

    assert (curve.status == -1)
    assert (curve.success == False)

def test_solve_ivp_inf():
    """
    """

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

    curve = solve_ivp(dynamics, tspan, y0)

    assert (curve.status == -2)
    assert (curve.success == False)

# Exponential decay

def test_solve_ivp_1_euclidean():
    """
    Tests solve_ivp against a classical problem with known solution.
    """

    m = -0.5
    tf = 10.0

    y0 = np.zeros((4,))
    y0[0] = 2.0
    y0[1] = 4.0
    y0[2] = 6.0
    y0[3] = -8.0

    # Do not modify anything below this line.

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

    assert (np.abs(curve.t[0] - 0.0) < TOL_FINE)
    assert (curve.ybar[0, 0] == y0[0])
    assert (curve.ybar[0, 1] == y0[1])
    assert (curve.ybar[0, 2] == y0[2])
    assert (curve.ybar[0, 3] == y0[3])


    assert (np.abs(curve.t[-1] - tf) < TOL_FINE)
    assert (np.abs(curve.ybar[-1, 0] - np.exp(m*tf)*y0[0]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 1] - np.exp(m*tf)*y0[1]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 2] - np.exp(m*tf)*y0[2]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 3] - np.exp(m*tf)*y0[3]) < 1e-5)

def test_solve_ivp_1_euclidean_event():
    """
    Tests solve_ivp against a classical problem with known solution.
    """

    m = -0.5
    tf = 10.0
    y2cross = 3.5

    y0 = np.zeros((4,))
    y0[0] = 2.0
    y0[1] = 4.0
    y0[2] = 6.0
    y0[3] = -8.0

    # Do not modify anything below this line.

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

    assert (np.abs(curve.t[0] - 0.0) < TOL_FINE)
    assert (curve.ybar[0, 0] == y0[0])
    assert (curve.ybar[0, 1] == y0[1])
    assert (curve.ybar[0, 2] == y0[2])
    assert (curve.ybar[0, 3] == y0[3])

    tcross = np.log(y2cross/y0[2])/m

    assert (np.abs(curve.t[-1] - tcross) < TOL_FINE)
    assert (np.abs(curve.ybar[-1, 0] - np.exp(m*tcross)*y0[0]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 1] - np.exp(m*tcross)*y0[1]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 2] - np.exp(m*tcross)*y0[2]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 3] - np.exp(m*tcross)*y0[3]) < 1e-5)

def test_solve_ivp_1():
    """
    Tests solve_ivp against a classical problem with known solution.
    """

    m = -0.5
    tf = 10.0

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    # Do not modify anything below this line.

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

    curve = solve_ivp(dynamics, tspan, y0)

    assert (np.abs(curve.t[0] - 0.0) < TOL_FINE)
    assert (curve.ybar[0, 0] == y0bar[0])
    assert (curve.ybar[0, 1] == y0bar[1])
    assert (curve.ybar[0, 2] == y0bar[2])
    assert (curve.ybar[0, 3] == y0bar[3])


    assert (np.abs(curve.t[-1] - tf) < TOL_FINE)
    assert (np.abs(curve.ybar[-1, 0] - np.exp(m*tf)*y0bar[0]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 1] - np.exp(m*tf)*y0bar[1]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 2] - np.exp(m*tf)*y0bar[2]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 3] - np.exp(m*tf)*y0bar[3]) < 1e-5)

def test_solve_ivp_1_event():
    """
    Tests solve_ivp against a classical problem with known solution.
    """

    m = -0.5
    tf = 10.0
    y2cross = 3.5

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    # Do not modify anything below this line.

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

    curve = solve_ivp(dynamics, tspan, y0)

    assert (np.abs(curve.t[0] - 0.0) < TOL_FINE)
    assert (curve.ybar[0, 0] == y0bar[0])
    assert (curve.ybar[0, 1] == y0bar[1])
    assert (curve.ybar[0, 2] == y0bar[2])
    assert (curve.ybar[0, 3] == y0bar[3])

    tcross = np.log(y2cross/y0bar[2])/m

    assert (np.abs(curve.t[-1] - tcross) < TOL_FINE)
    assert (np.abs(curve.ybar[-1, 0] - np.exp(m*tcross)*y0bar[0]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 1] - np.exp(m*tcross)*y0bar[1]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 2] - np.exp(m*tcross)*y0bar[2]) < 1e-5)
    assert (np.abs(curve.ybar[-1, 3] - np.exp(m*tcross)*y0bar[3]) < 1e-5)

def test_solve_ivp_1_segmented():
    """
    Tests solve_ivp against a classical problem with known solution.
    """

    m = -0.5
    tspan = np.linspace(0,10,41)

    y0bar = [2.0, 4.0, 6.0, -8.0]
    y0 = CompositeManifold([RN.from_vector(y0bar)])

    # Do not modify anything below this line.

    def eoms(t, y, m=m):
        ybar = y.serialize()
        dybar = np.zeros((4,))

        dybar[0] = m*ybar[0]
        dybar[1] = m*ybar[1]
        dybar[2] = m*ybar[2]
        dybar[3] = m*ybar[3]
        return CompositeAlgebra([rn.from_vector(dybar)])

    dynamics = HomogeneousIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, y0)

    for ii in range(len(tspan)):
        assert tspan[ii] in curve.t

# Lorenz equations

def test_solve_ivp_lorenz_euclidean():
    """
    Tests solve_ivp against a function that has not been wrapped.
    """

    def eoms(t, y):
        dy = np.zeros((3,))
        sigma = 10.0
        rho = 28.0
        b1 = 8.0
        b2 = 3.0
        beta = b1/b2

        dy[0] = -beta*y[0] + y[1]*y[2]
        dy[1] = -sigma*y[1] + sigma*y[2]
        dy[2] = -y[0]*y[1] + rho*y[1] - y[2]
        return dy

    y0 = np.zeros((3,))
    y0[0] = 25.0
    y0[1] = 0.0
    y0[2] = -20.0

    tspan = [0.0, 5.0]

    dynamics = EuclideanIVPSystem(eoms)

    curve = solve_ivp(dynamics, tspan, y0)

    assert np.abs(curve.t[0] - 0.0) < TOL_FINE
    assert curve.ybar[0, 0] == y0[0]
    assert curve.ybar[0, 1] == y0[1]
    assert curve.ybar[0, 2] == y0[2]

    assert np.abs(curve.t[-1] - 5.0) < TOL_FINE
    assert np.abs(curve.ybar[-1, 0] - 15.230) < 1e-1
    assert np.abs(curve.ybar[-1, 1] + 0.797) < 1e-1
    assert np.abs(curve.ybar[-1, 2] + 1.473) < 1e-1

def test_solve_ivp_lorenz_RNxRN_RN():
    """
    Tests solve_ivp against a function that has been wrapped.
    """

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

    assert np.abs(curve.t[0] - 0.0) < TOL_FINE
    assert curve.ybar[0, 0] == y0bar[0]
    assert curve.ybar[0, 1] == y0bar[1]
    assert curve.ybar[0, 2] == y0bar[2]

    assert np.abs(curve.t[-1] - 5.0) < TOL_FINE
    assert np.abs(curve.ybar[-1, 0] - 15.230) < 1e-1
    assert np.abs(curve.ybar[-1, 1] + 0.797) < 1e-1
    assert np.abs(curve.ybar[-1, 2] + 1.473) < 1e-1


def test_solve_ivp_lorenz_GLRxRN_RN():
    """
    Tests solve_ivp against a function with custom action GLR x RN -> RN.
    """

    def eoms(t, y):
        ybar = y.space[0].serialize()

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
        ghat = g.space[0].get_matrix()
        ybar = y.space[0].serialize()
        return CompositeManifold([RN.from_vector(np.dot(ghat, ybar))])

    y0bar = np.zeros((3,))
    y0bar[0] = 25.0
    y0bar[1] = 0.0
    y0bar[2] = -20.0

    tspan = np.linspace(0.0, 5.0, 20)

    dynamics = HomogeneousIVPSystem(eoms, action=action)

    curve = solve_ivp(dynamics, tspan, CompositeManifold([RN.from_vector(y0bar)]))

    assert np.abs(curve.t[0] - 0.0) < TOL_FINE
    assert curve.ybar[0, 0] == y0bar[0]
    assert curve.ybar[0, 1] == y0bar[1]
    assert curve.ybar[0, 2] == y0bar[2]

    assert np.abs(curve.t[-1] - 5.0) < TOL_FINE
    # assert np.abs(curve.ybar[-1, 0] - 15.230) < 1e-1  # TODO:
    # assert np.abs(curve.ybar[-1, 1] + 0.797) < 1e-1
    # assert np.abs(curve.ybar[-1, 2] + 1.473) < 1e-1

def test_solve_ivp_composite1():
    """
    Tests solve_ivp against a function with custom action (SE x GLR) x (SE x RN) -> (SE x RN).
    """

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

    x0bar = np.zeros((3,))
    x0bar[0] = 0.0
    x0bar[1] = 0.0
    x0bar[2] = np.pi/2.0
    x0 = se.from_vector(x0bar)

    lambda0bar = np.zeros((3,))
    lambda0bar[0] = 1.15407533e-03
    lambda0bar[1] = -3.17495766e+01
    lambda0bar[2] = -4.41935411e+00

    tspan = [0.0, 1.0]

    y0 = CompositeManifold([exp(x0), RN.from_vector(lambda0bar)])

    options = IVPOptions()
    # method = MuntheKaasZanna(Coefficients.RKDP54_7M); # TODO: Switch to RKV87r (or other default) once better error control is developed
    # options.method = method
    dynamics = HomogeneousIVPSystem(eoms, action=action)

    curve = solve_ivp(dynamics, tspan, y0, options)

    assert np.abs(curve.t[0] - 0.0) < TOL_FINE
    assert curve.ybar[0, 1] == -1.0
    assert curve.ybar[0, 3] == 1.0
    assert curve.ybar[0, 2] == x0bar[0]
    assert curve.ybar[0, 5] == x0bar[1]
    assert curve.ybar[0, 9] == lambda0bar[0]
    assert curve.ybar[0, 10] == lambda0bar[1]
    assert curve.ybar[0, 11] == lambda0bar[2]

    assert np.abs(curve.t[-1] - 1.0) < TOL_FINE
    # assert np.abs(curve.ybar[-1, 2] - 0.12732395447351627) < 1e-1 # TODO:
    # assert np.abs(curve.ybar[-1, 5] + 0.0) < 1e-1
    # assert np.abs(curve.ybar[-1, 9] - lambda0bar[0]) < 1e-1
    # assert np.abs(curve.ybar[-1, 10] + lambda0bar[1]) < 1e-1
    # assert np.abs(curve.ybar[-1, 11] - lambda0bar[2]) < 1e-1


# def test_MuntheKaas_vfex1():
#     """
#     Tests the MuntheKaas function with vfex1.
#     """

#     from lielab.domain import SO, CompositeManifold
#     from lielab.dynamics import vfex1
#     from lielab.integrate import MuntheKaas

#     ts = MuntheKaas()
#     vf = vfex1()

#     y0 = SO(3)
#     M0 = CompositeManifold([y0])

#     out = ts(vf, M0, 1.0, 0.02)

#     y1 = out.low.space[0]

#     abs(y1(0,0) - 0.999596034819844) < TOL_FINE
#     abs(y1(0,1) - 0.020398551422611) < TOL_FINE
#     abs(y1(0,2) - 0.019790560181659) < TOL_FINE
#     abs(y1(1,0) + 0.019990512514326) < TOL_FINE
#     abs(y1(1,1) - 0.999587901244249) < TOL_FINE
#     abs(y1(1,2) + 0.020601143063711) < TOL_FINE
#     abs(y1(2,0) + 0.020202637992582) < TOL_FINE
#     abs(y1(2,1) - 0.020197197478265) < TOL_FINE
#     abs(y1(2,2) - 0.999591880035129) < TOL_FINE


# def test_MuntheKaas_vfex2():
#     """
#     Tests the MuntheKaas function with vfex2.
#     """

#     from lielab.domain import RN, CompositeManifold
#     from lielab.dynamics import vfex2
#     from lielab.integrate import MuntheKaas

#     ts = MuntheKaas()
#     vf = vfex2()

#     y0 = RN.from_vector([25.0, 0.0, -20.0])
#     M0 = CompositeManifold([y0])

#     out = ts(vf, M0, 0.0, 0.02)

#     y1 = out.low.space[0]

#     assert abs(y1(0) - 24.425986197956878) < TOL_FINE
#     assert abs(y1(1) + 3.596428324678375) < TOL_FINE
#     assert abs(y1(2) + 19.733395791914329) < TOL_FINE


# def test_Flow_fails():
#     """
#     Tests cases where Flow should fail.
#     """

#     from lielab.domain import SO, CompositeManifold
#     from lielab.dynamics import vfex1
#     from lielab.integrate import Flow

#     y0 = SO(3)
#     M0 = CompositeManifold([y0])
#     vf = vfex1()
#     tspan = []
#     f = Flow()

#     # tspan has no values
#     with pytest.raises(Exception):
#         f(vf, tspan, M0)
    
#     # tspan has 1 value
#     tspan.append(0.0)
#     with pytest.raises(Exception):
#         f(vf, tspan, M0)


# def test_solve_ivp_copy_output():
#     """
#     Tests that solve_ivp outputs are copied.
#     """

#     from lielab.integrate import solve_ivp

#     def vf(t, y):
#         return np.array([y[1], y[0]])

#     tspan = [0.0, 5.5]

#     y0_1 = [1.0, 0.0]
#     y0_2 = [2.0, 0.0]

#     curve1 = solve_ivp(vf, tspan, y0_1)
#     curve2 = solve_ivp(vf, tspan, y0_2)

#     assert abs(curve1.ybar[-1, 0] - curve2.ybar[-1, 0]) >= 1e-4
#     assert abs(curve1.ybar[-1, 1] - curve2.ybar[-1, 1]) >= 1e-4


# def test_Flow_vfex2_rk45_fixed():
#     from lielab.domain import RN, CompositeManifold
#     from lielab.dynamics import vfex2
#     from lielab.integrate import Flow

#     vf = vfex2()
#     flow = Flow()
#     flow.variable_time_step = False
#     flow.dt = 0.02

#     tspan = [0.0, 5.0]
#     y0 = RN.from_vector([25.0, 0.0, -20.0])
#     M0 = CompositeManifold([y0])

#     curve = flow(vf, tspan, M0)

#     assert abs(curve.t[0] - tspan[0]) < TOL_FINE
#     assert curve.ybar[0, 0] == y0(0)
#     assert curve.ybar[0, 1] == y0(1)
#     assert curve.ybar[0, 2] == y0(2)

#     nrows = curve.ybar.shape[0]

#     assert abs(curve.t[nrows-1] - 5.0) < TOL_FINE
#     assert abs(curve.ybar[nrows-1, 0] - 15.210570567999987) < TOL_FINE
#     assert abs(curve.ybar[nrows-1, 1] + 0.788689660918195) < TOL_FINE
#     assert abs(curve.ybar[nrows-1, 2] + 1.459476938449221) < TOL_FINE


# def test_Flow_vfex2_rk45_variable():
#     from lielab.domain import RN, CompositeManifold
#     from lielab.dynamics import vfex2
#     from lielab.integrate import Flow

#     vf = vfex2()
#     flow = Flow()
#     flow.variable_time_step = True
#     flow.dt = 0.02

#     tspan = [0.0, 5.0]
#     y0 = RN.from_vector([25.0, 0.0, -20.0])
#     M0 = CompositeManifold([y0])

#     curve = flow(vf, tspan, M0)

#     assert abs(curve.t[0] - 0.0) < TOL_FINE
#     assert abs(curve.t[1] - 0.007703769593747) < TOL_COARSE
#     assert abs(curve.t[2] - 0.015420629134474) < TOL_COARSE
#     assert abs(curve.t[3] - 0.023255332563845) < TOL_COARSE
#     assert abs(curve.t[4] - 0.031246577586516) < TOL_COARSE

#     assert abs(curve.t[0] - tspan[0]) < TOL_FINE
#     assert curve.ybar[0, 0] == y0(0)
#     assert curve.ybar[0, 1] == y0(1)
#     assert curve.ybar[0, 2] == y0(2)

#     nrows = curve.y.shape[0]

#     assert abs(curve.t[nrows-1] - 5.0) < TOL_FINE
#     assert abs(curve.ybar[nrows-1, 0] - 15.230102737555342) < TOL_COARSE
#     assert abs(curve.ybar[nrows-1, 1] + 0.796697875936802) < TOL_COARSE
#     assert abs(curve.ybar[nrows-1, 2] + 1.472989006310112) < TOL_COARSE

# def test_solve_ivp_unwrapped_rk45_variable():
#     """
#     Tests solve_ivp against a function that has not been wrapped.
#     """

#     from lielab.integrate import solve_ivp
    
#     y0 = np.array([25.0, 0.0, -20.0])

#     tspan = [0.0, 5.0]

#     curve = solve_ivp(eoms_lorenz_unwrapped, tspan, y0)

#     assert abs(curve.t[0] - 0.0) < TOL_FINE
#     assert abs(curve.t[1] - 0.007703769593747) < TOL_COARSE
#     assert abs(curve.t[2] - 0.015420629134474) < TOL_COARSE
#     assert abs(curve.t[3] - 0.023255332563845) < TOL_COARSE
#     assert abs(curve.t[4] - 0.031246577586516) < TOL_COARSE

#     assert abs(curve.t[0] - tspan[0]) < TOL_FINE
#     assert curve.ybar[0, 0] == y0[0]
#     assert curve.ybar[0, 1] == y0[1]
#     assert curve.ybar[0, 2] == y0[2]

#     nrows = curve.ybar.shape[0]

#     assert abs(curve.t[nrows-1] - 5.0) < TOL_FINE
#     assert abs(curve.ybar[nrows-1, 0] - 15.230102737555342) < TOL_COARSE
#     assert abs(curve.ybar[nrows-1, 1] + 0.796697875936802) < TOL_COARSE
#     assert abs(curve.ybar[nrows-1, 2] + 1.472989006310112) < TOL_COARSE


# def test_solve_ivp_A1_vector():
#     """
#     Problem A-1
    
#     Uses VectorXd representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """

#     from lielab.integrate import solve_ivp
#     import numpy as np

#     def vf(t, y):
#         m = 1.0/4.0
#         w = 8.0
#         k = 2.0
#         dy = np.array([y[1], (w - k*y[1])/m])
#         return dy

#     def vf_event(t, y):
#         H = 10.0
#         return H - y[0]

#     tspan = np.array([0.0, 5.0])
#     y0 = np.array([0.0, 0.0])

#     out = solve_ivp(vf, tspan, y0, vf_event)

#     assert abs(out.t[0] - 0.0) <= 1.0e-9
#     assert abs(out.t[-1] - 2.62499999990522) <= 1.0e-9

#     for ii in range(len(out.t)):
#         assert abs(out.ybar[ii, 0] - (4.0*(out.t[ii] + 1.0/8.0*np.exp(-8.0*out.t[ii]) - 1.0/8.0))) <= 1e-6
#         assert abs(out.ybar[ii, 1] - (4.0*(1.0 - np.exp(-8.0*out.t[ii])))) <= 1e-6


# def test_solve_ivp_A1_hom():
#     """
#     Problem A-1
    
#     Uses hom representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """
    
#     from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
#     from lielab.integrate import solve_ivp
#     import numpy as np

#     def vf(t, M):
#         m = 1.0/4.0
#         w = 8.0
#         k = 2.0
#         y = M.space[0]
#         dy = rn.from_vector([y(1), (w - k*y(1))/m])
#         return CompositeAlgebra([dy])

#     def vf_event(t, M):
#         H = 10.0
#         y = M.space[0]
#         return H - y(0)

#     tspan = np.array([0.0, 5.0])
#     a1M0 = CompositeManifold([RN.from_vector([0.0, 0.0])])

#     out = solve_ivp(vf, tspan, a1M0, vf_event)

#     assert abs(out.t[0] - 0.0) <= 1.0e-9
#     assert abs(out.t[-1] - 2.62499999990522) <= 1.0e-9

#     for ii in range(len(out.t)):
#         assert abs(out.ybar[ii, 0] - (4.0*(out.t[ii] + 1.0/8.0*np.exp(-8.0*out.t[ii]) - 1.0/8.0))) <= 1e-6
#         assert abs(out.ybar[ii, 1] - (4.0*(1.0 - np.exp(-8.0*out.t[ii])))) <= 1e-6


# def test_solve_ivp_A2_vector():
#     """
#     Problem A-2
    
#     Uses VectorXd representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """

#     from lielab.integrate import solve_ivp
#     import numpy as np
#     from math import acos

#     g = 32
#     R = 4000*5280
#     H = 237000*5280

#     def vf(t, y):
#         dy = np.array([-np.sqrt(2*g*R**2) * np.sqrt((H - y[0])/(H*y[0]))])
#         return dy

#     def vf_event(t, y):
#         return y[0] - R
    
#     def h_to_t(h):
#         return (H**(3/2) / (8*h))*(np.sqrt(h/H - (h/H)**2) + 1/2*acos(2*h/H - 1))

#     out = solve_ivp(vf, [0.0, 800000.0], [H - 1e-6], vf_event, dt_max=10000)

#     assert abs(out.t[0] - 0.0) <= 1.0e-14
#     assert abs(out.t[-1] - h_to_t(out.ybar[-1, 0])) <= 1e-1


# def test_solve_ivp_A2_hom():
#     """
#     Problem A-2
    
#     Uses hom representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """

#     from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
#     from lielab.integrate import solve_ivp
#     import numpy as np
#     from math import acos

#     g = 32
#     R = 4000*5280
#     H = 237000*5280

#     def vf(t, M):
#         y = M.space[0]
#         dy = rn.from_vector([-np.sqrt(2*g*R**2) * np.sqrt((H - y(0))/(H*y(0)))])
#         return CompositeAlgebra([dy])

#     def vf_event(t, M):
#         y = M.space[0]
#         return y(0) - R
    
#     def h_to_t(h):
#         return (H**(3/2) / (8*h))*(np.sqrt(h/H - (h/H)**2) + 1/2*acos(2*h/H - 1))
    
#     out = solve_ivp(vf, [0.0, 800000.0], CompositeManifold([RN.from_vector([H - 1e-6])]), vf_event, dt_max=10000)

#     assert abs(out.t[0] - 0.0) <= 1.0e-14
#     assert abs(out.t[-1] - h_to_t(out.ybar[-1, 0])) <= 1e-1


# def test_solve_ivp_A4_vector():
#     """
#     Problem A-4
    
#     Uses VectorXd representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """

#     from lielab.integrate import solve_ivp
#     import numpy as np

#     def vf(t, y):
#         y1d = y[1]
#         y2d = -(16.0 * np.pi**2 * np.exp(-2.0*t) - 1.0/4.0)*y[0]
        
#         dy = np.array([y1d, y2d])
#         return dy
    
#     def t_to_y1(t):
#         return np.exp(t/2.0)*np.cos(4*np.pi*np.exp(-t))

#     def t_to_y2(t):
#         return np.exp(t/2)*(4*np.pi*np.exp(-t)*np.sin(4*np.pi*np.exp(-t)) + 1/2*np.cos(4*np.pi*np.exp(-t)))
    
#     out = solve_ivp(vf, [0.0, np.log(8) - np.log(1)], [1.0, 0.5])

#     assert abs(out.t[0] - 0.0) <= 1.0e-14
#     assert abs(out.ybar[-1, 0] - 0.0) <= 1e-5

#     for ii in range(len(out.t)):
#         assert abs(out.ybar[ii, 0] - t_to_y1(out.t[ii])) <= 1e-5
#         assert abs(out.ybar[ii, 1] - t_to_y2(out.t[ii])) <= 1e-5


# def test_solve_ivp_A4_hom():
#     """
#     Problem A-4
    
#     Uses hom representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """

#     from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
#     from lielab.integrate import solve_ivp
#     import numpy as np

#     def vf(t, M):
#         y = M.space[0]
#         y1d = y(1)
#         y2d = -(16.0 * np.pi**2 * np.exp(-2.0*t) - 1.0/4.0)*y(0)
#         return CompositeAlgebra([rn.from_vector([y1d, y2d])])
    
#     def t_to_y1(t):
#         return np.exp(t/2.0)*np.cos(4*np.pi*np.exp(-t))

#     def t_to_y2(t):
#         return np.exp(t/2)*(4*np.pi*np.exp(-t)*np.sin(4*np.pi*np.exp(-t)) + 1/2*np.cos(4*np.pi*np.exp(-t)))
    
#     out = solve_ivp(vf, [0.0, np.log(8) - np.log(1)], CompositeManifold([RN.from_vector([1.0, 0.5])]))

#     assert abs(out.t[0] - 0.0) <= 1.0e-14
#     assert abs(out.ybar[-1, 0] - 0.0) <= 1e-5

#     for ii in range(len(out.t)):
#         assert abs(out.ybar[ii, 0] - t_to_y1(out.t[ii])) <= 1e-5
#         assert abs(out.ybar[ii, 1] - t_to_y2(out.t[ii])) <= 1e-5


# def test_solve_ivp_A5_vector():
#     """
#     Problem A-5
    
#     Uses VectorXd representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """

#     from lielab.integrate import solve_ivp
#     import numpy as np

#     def vf(t, y):
#         y1d = y[1]
#         y2d = -16.0 * np.cos(np.pi*t/2.0)*y[1] - (64.0 * np.pi**2 + 64.0*np.cos(np.pi*t/2.0)**2 - 4.0*np.pi*np.sin(np.pi*t/2.0))*y[0]
#         return np.array([y1d, y2d])
    
#     def vf_event(t, y):
#         if (t > 1.07) and (t < 1.30):
#             return -y[0]
#         return 1
    
#     def t_to_y1(t):
#         return np.exp(-16/np.pi*np.sin(np.pi*t/2))*np.cos(8*np.pi*t)
    
#     out = solve_ivp(vf, [0.0, 2.0], [1.0, -8.0], vf_event)

#     assert abs(out.t[0] - 0.0) <= 1.0e-14
#     assert abs(out.t[-1] - (2*10-1)/16.0) <= 1e-6

#     for ii in range(len(out.t)):
#         assert abs(out.ybar[ii, 0] - t_to_y1(out.t[ii])) <= 1e-6


# def test_solve_ivp_A5_hom():
#     """
#     Problem A-5
    
#     Uses hom representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """

#     from lielab.domain import rn, RN, CompositeAlgebra, CompositeManifold
#     from lielab.integrate import solve_ivp
#     import numpy as np

#     def vf(t, M):
#         y = M.space[0]
#         y1d = y(1)
#         y2d = -16.0 * np.cos(np.pi*t/2.0)*y(1) - (64.0 * np.pi**2 + 64.0*np.cos(np.pi*t/2.0)**2 - 4.0*np.pi*np.sin(np.pi*t/2.0))*y(0)
#         return CompositeAlgebra([rn.from_vector([y1d, y2d])])
    
#     def vf_event(t, M):
#         y = M.space[0]
#         if (t > 1.07) and (t < 1.30):
#             return -y(0)
#         return 1
    
#     def t_to_y1(t):
#         return np.exp(-16/np.pi*np.sin(np.pi*t/2))*np.cos(8*np.pi*t)
    
#     out = solve_ivp(vf, [0.0, 2.0], CompositeManifold([RN.from_vector([1.0, -8.0])]), vf_event)

#     assert abs(out.t[0] - 0.0) <= 1.0e-14
#     assert abs(out.t[-1] - (2*10-1)/16.0) <= 1e-6

#     for ii in range(len(out.t)):
#         assert abs(out.ybar[ii, 0] - t_to_y1(out.t[ii])) <= 1e-6


# def test_solve_ivp_B1_vector():
#     """
#     Problem B-1
    
#     Uses VectorXd representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """

#     from lielab.integrate import solve_ivp
#     import numpy as np

#     K = 5.0/3.0
#     L = 1.0/4.0

#     def vf(t, y):
#         E = y[0] - y[1]

#         dy0 = y[0]
#         dy1 = 0
#         if (E < -L/K):
#             dy1 = -L
#         elif ((-L/K <= E) and (E <= L/K)):
#             dy1 = K*E
#         elif (L/K < E):
#             dy1 = L

#         dy = np.array([dy0, dy1])
#         return dy
    
    
#     t1 = 0.1569
#     y2t1 = 1.0199
    
#     def t_to_y1(t):
#         return np.exp(t)

#     def t_to_y2(t):
#         if (t < t1):
#             return K/(K+1)*np.exp(t) + 1/(K+1)*np.exp(-K*t)
#         return L*(t - t1) + y2t1
    
#     out = solve_ivp(vf, [0.0, 0.5], [1.0, 1.0])

#     assert abs(out.t[0] - 0.0) <= 1.0e-14
#     assert abs(out.t[-1] - 0.5) <= 1.0e-14

#     for ii in range(len(out.t)):
#         assert abs(out.ybar[ii, 0] - t_to_y1(out.t[ii])) <= 1e-7
#         assert abs(out.ybar[ii, 1] - t_to_y2(out.t[ii])) <= 1e-4 # Answer given in the document is only good to 1e-4 (see t1 and y2(t1))


# def test_solve_ivp_B1_hom():
#     """
#     Problem B-1
    
#     Uses hom representation and events.
     
#     Source: Thompson, S. A collection of test problems for ordinary differential
#     equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
#     Oak Ridge National Lab., TN (USA), 1987.
#     """

#     from lielab.domain import CompositeAlgebra, CompositeManifold, rn, RN, glr
#     from lielab.integrate import solve_ivp, CustomMuntheKaas
#     import numpy as np

#     method = CustomMuntheKaas()

#     def action(G, M):
#         _G0 = G.space[0]
#         _G1 = G.space[1]
#         _Y0 = M.space[0]
#         _Y1 = M.space[1]

#         Y0next = RN.from_vector(np.dot(_G0.data, _Y0.data))
#         Y1next = _G1*_Y1
#         return CompositeManifold([Y0next, Y1next])
    
#     method.action = action

#     K = 5.0/3.0
#     L = 1.0/4.0

#     def vf(t, M):
#         y0 = M.space[0](0)
#         y1 = M.space[1](0)

#         E = y0 - y1

#         if (E < -L/K):
#             dy1 = -L
#         elif ((-L/K <= E) and (E <= L/K)):
#             dy1 = K*E
#         elif (L/K < E):
#             dy1 = L
#         else:
#             dy1 = 0

#         return CompositeAlgebra([glr.basis(0,1), dy1*rn.basis(0,2)])
    
#     t1 = 0.1569
#     y2t1 = 1.0199
    
#     def t_to_y1(t):
#         return np.exp(t)

#     def t_to_y2(t):
#         if (t < t1):
#             return K/(K+1)*np.exp(t) + 1/(K+1)*np.exp(-K*t)
#         return L*(t - t1) + y2t1
    
#     out = solve_ivp(vf, [0.0, 0.5], CompositeManifold([RN.from_vector([1.0]), RN.from_vector([1.0])]), method=method)

#     assert abs(out.t[0] - 0.0) <= 1.0e-14
#     assert abs(out.t[-1] - 0.5) <= 1.0e-14

#     for ii in range(len(out.t)):
#         assert abs(out.ybar[ii, 0] - t_to_y1(out.t[ii])) <= 1e-14 # Answer is analytic
#         assert abs(out.ybar[ii, 1] - t_to_y2(out.t[ii])) <= 1e-4 # Answer given in the document is only good to 1e-4 (see t1 and y2(t1))
