import numpy as np
import pytest

from lielab.testing import *

def eoms_lorenz_unwrapped(t, X):
    dx = np.array([0.0,0.0,0.0])
    sigma = 10.0
    rho = 28.0
    b1 = 8.0
    b2 = 3.0
    beta = b1/b2

    dx[0] = -beta*X[0] + X[1]*X[2]
    dx[1] = -sigma*X[1] + sigma*X[2]
    dx[2] = -X[0]*X[1] + rho*X[1] - X[2]
    return dx

def test_MuntheKaas_vfex1():
    """
    Tests the MuntheKaas function with vfex1.
    """

    from lielab.domain import SO, hmlie
    from lielab.dynamics import vfex1
    from lielab.topos import MuntheKaas

    ts = MuntheKaas()
    vf = vfex1()

    y0 = SO(3)
    M0 = hmlie([y0])

    out = ts(vf, M0, 1.0, 0.02)

    y1 = out.low.space[0]

    abs(y1(0,0) - 0.999596034819844) < TOL_FINE
    abs(y1(0,1) - 0.020398551422611) < TOL_FINE
    abs(y1(0,2) - 0.019790560181659) < TOL_FINE
    abs(y1(1,0) + 0.019990512514326) < TOL_FINE
    abs(y1(1,1) - 0.999587901244249) < TOL_FINE
    abs(y1(1,2) + 0.020601143063711) < TOL_FINE
    abs(y1(2,0) + 0.020202637992582) < TOL_FINE
    abs(y1(2,1) - 0.020197197478265) < TOL_FINE
    abs(y1(2,2) - 0.999591880035129) < TOL_FINE


def test_MuntheKaas_vfex2():
    """
    Tests the MuntheKaas function with vfex2.
    """

    from lielab.domain import RN, hmlie
    from lielab.dynamics import vfex2
    from lielab.topos import MuntheKaas

    ts = MuntheKaas()
    vf = vfex2()

    y0 = RN([25.0, 0.0, -20.0])
    M0 = hmlie([y0])

    out = ts(vf, M0, 0.0, 0.02)

    y1 = out.low.space[0]

    assert abs(y1(0) - 24.425986197956878) < TOL_FINE
    assert abs(y1(1) + 3.596428324678375) < TOL_FINE
    assert abs(y1(2) + 19.733395791914329) < TOL_FINE


def test_Flow_fails():
    """
    Tests cases where Flow should fail.
    """

    from lielab.domain import SO, hmlie
    from lielab.dynamics import vfex1
    from lielab.topos import Flow

    y0 = SO(3)
    M0 = hmlie([y0])
    vf = vfex1()
    tspan = []
    f = Flow()

    # tspan has no values
    with pytest.raises(Exception):
        f(vf, tspan, M0)
    
    # tspan has 1 value
    tspan.append(0.0)
    with pytest.raises(Exception):
        f(vf, tspan, M0)


def test_Flow_copy_output():
    """
    Tests that Flows outputs are copied.
    """

    from lielab.topos import Flow

    f = Flow()

    def vf(t, y):
        return np.array([y[1], y[0]])

    tspan = [0.0, 5.5]

    y0_1 = [1.0, 0.0]
    y0_2 = [2.0, 0.0]

    curve1 = f(vf, tspan, y0_1)
    curve2 = f(vf, tspan, y0_2)

    assert abs(curve1.y[-1, 0] - curve2.y[-1, 0]) >= 1e-4
    assert abs(curve1.y[-1, 1] - curve2.y[-1, 1]) >= 1e-4


def test_Flow_vfex2_rk45_fixed():
    from lielab.domain import RN, hmlie
    from lielab.dynamics import vfex2
    from lielab.topos import Flow

    vf = vfex2()
    flow = Flow()
    flow.variable_time_step = False
    flow.dt = 0.02

    tspan = [0.0, 5.0]
    y0 = RN([25.0, 0.0, -20.0])
    M0 = hmlie([y0])

    curve = flow(vf, tspan, M0)

    assert abs(curve.t[0] - tspan[0]) < TOL_FINE
    assert curve.y[0, 0] == y0(0)
    assert curve.y[0, 1] == y0(1)
    assert curve.y[0, 2] == y0(2)

    nrows = curve.y.shape[0]

    assert abs(curve.t[nrows-1] - 5.0) < TOL_FINE
    assert abs(curve.y[nrows-1, 0] - 15.210570567999987) < TOL_FINE
    assert abs(curve.y[nrows-1, 1] + 0.788689660918195) < TOL_FINE
    assert abs(curve.y[nrows-1, 2] + 1.459476938449221) < TOL_FINE


def test_Flow_vfex2_rk45_variable():
    from lielab.domain import RN, hmlie
    from lielab.dynamics import vfex2
    from lielab.topos import Flow

    vf = vfex2()
    flow = Flow()
    flow.variable_time_step = True
    flow.dt = 0.02

    tspan = [0.0, 5.0]
    y0 = RN([25.0, 0.0, -20.0])
    M0 = hmlie([y0])

    curve = flow(vf, tspan, M0)

    assert abs(curve.t[0] - 0.0) < TOL_FINE
    assert abs(curve.t[1] - 0.007703769593747) < TOL_COARSE
    assert abs(curve.t[2] - 0.015420629134474) < TOL_COARSE
    assert abs(curve.t[3] - 0.023255332563845) < TOL_COARSE
    assert abs(curve.t[4] - 0.031246577586516) < TOL_COARSE

    assert abs(curve.t[0] - tspan[0]) < TOL_FINE
    assert curve.y[0, 0] == y0(0)
    assert curve.y[0, 1] == y0(1)
    assert curve.y[0, 2] == y0(2)

    nrows = curve.y.shape[0]

    assert abs(curve.t[nrows-1] - 5.0) < TOL_FINE
    assert abs(curve.y[nrows-1, 0] - 15.230102737555342) < TOL_COARSE
    assert abs(curve.y[nrows-1, 1] + 0.796697875936802) < TOL_COARSE
    assert abs(curve.y[nrows-1, 2] + 1.472989006310112) < TOL_COARSE

def test_Flow_unwrapped_rk45_variable():
    """
    Tests Flow against a function that has not been wrapped.
    """

    from lielab.topos import Flow
    
    y0 = np.array([25.0, 0.0, -20.0])

    tspan = [0.0, 5.0]

    f = Flow()
    f.variable_time_step = True
    f.dt = 0.02
    curve = f(eoms_lorenz_unwrapped, tspan, y0)

    assert abs(curve.t[0] - 0.0) < TOL_FINE
    assert abs(curve.t[1] - 0.007703769593747) < TOL_COARSE
    assert abs(curve.t[2] - 0.015420629134474) < TOL_COARSE
    assert abs(curve.t[3] - 0.023255332563845) < TOL_COARSE
    assert abs(curve.t[4] - 0.031246577586516) < TOL_COARSE

    assert abs(curve.t[0] - tspan[0]) < TOL_FINE
    assert curve.y[0, 0] == y0[0]
    assert curve.y[0, 1] == y0[1]
    assert curve.y[0, 2] == y0[2]

    nrows = curve.y.shape[0]

    assert abs(curve.t[nrows-1] - 5.0) < TOL_FINE
    assert abs(curve.y[nrows-1, 0] - 15.230102737555342) < TOL_COARSE
    assert abs(curve.y[nrows-1, 1] + 0.796697875936802) < TOL_COARSE
    assert abs(curve.y[nrows-1, 2] + 1.472989006310112) < TOL_COARSE


def test_Flow_A1_vector():
    """
    Problem A-1
    
    Uses VectorXd representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """

    from lielab.topos import Flow
    import numpy as np

    f = Flow()

    def vf(t, y):
        m = 1.0/4.0
        w = 8.0
        k = 2.0
        dy = np.array([y[1], (w - k*y[1])/m])
        return dy

    def vf_event(t, y):
        H = 10.0
        return H - y[0]

    tspan = np.array([0.0, 5.0])
    y0 = np.array([0.0, 0.0])

    out = f(vf, tspan, y0, vf_event)

    assert abs(out.t[0] - 0.0) <= 1.0e-9
    assert abs(out.t[-1] - 2.62499999990522) <= 1.0e-9

    for ii in range(len(out.t)):
        assert abs(out.y[ii, 0] - (4.0*(out.t[ii] + 1.0/8.0*np.exp(-8.0*out.t[ii]) - 1.0/8.0))) <= 1e-6
        assert abs(out.y[ii, 1] - (4.0*(1.0 - np.exp(-8.0*out.t[ii])))) <= 1e-6


def test_Flow_A1_hom():
    """
    Problem A-1
    
    Uses hom representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """
    
    from lielab.domain import rn, RN, halie, hmlie
    from lielab.topos import Flow
    import numpy as np

    f = Flow()

    def vf(t, M):
        m = 1.0/4.0
        w = 8.0
        k = 2.0
        y = M.space[0]
        dy = rn([y(1), (w - k*y(1))/m])
        return halie([dy])

    def vf_event(t, M):
        H = 10.0
        y = M.space[0]
        return H - y(0)

    tspan = np.array([0.0, 5.0])
    a1M0 = hmlie([RN([0.0, 0.0])])

    out = f(vf, tspan, a1M0, vf_event)

    assert abs(out.t[0] - 0.0) <= 1.0e-9
    assert abs(out.t[-1] - 2.62499999990522) <= 1.0e-9

    for ii in range(len(out.t)):
        assert abs(out.y[ii, 0] - (4.0*(out.t[ii] + 1.0/8.0*np.exp(-8.0*out.t[ii]) - 1.0/8.0))) <= 1e-6
        assert abs(out.y[ii, 1] - (4.0*(1.0 - np.exp(-8.0*out.t[ii])))) <= 1e-6


def test_Flow_A2_vector():
    """
    Problem A-2
    
    Uses VectorXd representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """

    from lielab.topos import Flow
    import numpy as np
    from math import acos

    f = Flow()
    f.dt_max = 10000

    g = 32
    R = 4000*5280
    H = 237000*5280

    def vf(t, y):
        dy = np.array([-np.sqrt(2*g*R**2) * np.sqrt((H - y[0])/(H*y[0]))])
        return dy

    def vf_event(t, y):
        return y[0] - R
    
    def h_to_t(h):
        return (H**(3/2) / (8*h))*(np.sqrt(h/H - (h/H)**2) + 1/2*acos(2*h/H - 1))
    
    out = f(vf, [0.0, 800000.0], [H - 1e-6], vf_event)

    assert abs(out.t[0] - 0.0) <= 1.0e-14
    assert abs(out.t[-1] - h_to_t(out.y[-1, 0])) <= 1e-1


def test_Flow_A2_hom():
    """
    Problem A-2
    
    Uses hom representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """

    from lielab.domain import rn, RN, halie, hmlie
    from lielab.topos import Flow
    import numpy as np
    from math import acos

    f = Flow()
    f.dt_max = 10000

    g = 32
    R = 4000*5280
    H = 237000*5280

    def vf(t, M):
        y = M.space[0]
        dy = rn([-np.sqrt(2*g*R**2) * np.sqrt((H - y(0))/(H*y(0)))])
        return halie([dy])

    def vf_event(t, M):
        y = M.space[0]
        return y(0) - R
    
    def h_to_t(h):
        return (H**(3/2) / (8*h))*(np.sqrt(h/H - (h/H)**2) + 1/2*acos(2*h/H - 1))
    
    out = f(vf, [0.0, 800000.0], hmlie([RN([H - 1e-6])]), vf_event)

    assert abs(out.t[0] - 0.0) <= 1.0e-14
    assert abs(out.t[-1] - h_to_t(out.y[-1, 0])) <= 1e-1


def test_Flow_A4_vector():
    """
    Problem A-4
    
    Uses VectorXd representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """

    from lielab.topos import Flow
    import numpy as np

    f = Flow()

    def vf(t, y):
        y1d = y[1]
        y2d = -(16.0 * np.pi**2 * np.exp(-2.0*t) - 1.0/4.0)*y[0]
        
        dy = np.array([y1d, y2d])
        return dy
    
    def t_to_y1(t):
        return np.exp(t/2.0)*np.cos(4*np.pi*np.exp(-t))

    def t_to_y2(t):
        return np.exp(t/2)*(4*np.pi*np.exp(-t)*np.sin(4*np.pi*np.exp(-t)) + 1/2*np.cos(4*np.pi*np.exp(-t)))
    
    out = f(vf, [0.0, np.log(8) - np.log(1)], [1.0, 0.5])

    assert abs(out.t[0] - 0.0) <= 1.0e-14
    assert abs(out.y[-1, 0] - 0.0) <= 1e-5

    for ii in range(len(out.t)):
        assert abs(out.y[ii, 0] - t_to_y1(out.t[ii])) <= 1e-5
        assert abs(out.y[ii, 1] - t_to_y2(out.t[ii])) <= 1e-5


def test_Flow_A4_hom():
    """
    Problem A-4
    
    Uses hom representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """

    from lielab.domain import rn, RN, halie, hmlie
    from lielab.topos import Flow
    import numpy as np

    f = Flow()

    def vf(t, M):
        y = M.space[0]
        y1d = y(1)
        y2d = -(16.0 * np.pi**2 * np.exp(-2.0*t) - 1.0/4.0)*y(0)
        return halie([rn([y1d, y2d])])
    
    def t_to_y1(t):
        return np.exp(t/2.0)*np.cos(4*np.pi*np.exp(-t))

    def t_to_y2(t):
        return np.exp(t/2)*(4*np.pi*np.exp(-t)*np.sin(4*np.pi*np.exp(-t)) + 1/2*np.cos(4*np.pi*np.exp(-t)))
    
    out = f(vf, [0.0, np.log(8) - np.log(1)], hmlie([RN([1.0, 0.5])]))

    assert abs(out.t[0] - 0.0) <= 1.0e-14
    assert abs(out.y[-1, 0] - 0.0) <= 1e-5

    for ii in range(len(out.t)):
        assert abs(out.y[ii, 0] - t_to_y1(out.t[ii])) <= 1e-5
        assert abs(out.y[ii, 1] - t_to_y2(out.t[ii])) <= 1e-5


def test_Flow_A5_vector():
    """
    Problem A-5
    
    Uses VectorXd representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """

    from lielab.topos import Flow
    import numpy as np

    f = Flow()

    def vf(t, y):
        y1d = y[1]
        y2d = -16.0 * np.cos(np.pi*t/2.0)*y[1] - (64.0 * np.pi**2 + 64.0*np.cos(np.pi*t/2.0)**2 - 4.0*np.pi*np.sin(np.pi*t/2.0))*y[0]
        return np.array([y1d, y2d])
    
    def vf_event(t, y):
        if (t > 1.07) and (t < 1.30):
            return -y[0]
        return 1
    
    def t_to_y1(t):
        return np.exp(-16/np.pi*np.sin(np.pi*t/2))*np.cos(8*np.pi*t)
    
    out = f(vf, [0.0, 2.0], [1.0, -8.0], vf_event)

    assert abs(out.t[0] - 0.0) <= 1.0e-14
    assert abs(out.t[-1] - (2*10-1)/16.0) <= 1e-6

    for ii in range(len(out.t)):
        assert abs(out.y[ii, 0] - t_to_y1(out.t[ii])) <= 1e-6


def test_Flow_A5_hom():
    """
    Problem A-5
    
    Uses hom representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """

    from lielab.domain import rn, RN, halie, hmlie
    from lielab.topos import Flow
    import numpy as np

    f = Flow()

    def vf(t, M):
        y = M.space[0]
        y1d = y(1)
        y2d = -16.0 * np.cos(np.pi*t/2.0)*y(1) - (64.0 * np.pi**2 + 64.0*np.cos(np.pi*t/2.0)**2 - 4.0*np.pi*np.sin(np.pi*t/2.0))*y(0)
        return halie([rn([y1d, y2d])])
    
    def vf_event(t, M):
        y = M.space[0]
        if (t > 1.07) and (t < 1.30):
            return -y(0)
        return 1
    
    def t_to_y1(t):
        return np.exp(-16/np.pi*np.sin(np.pi*t/2))*np.cos(8*np.pi*t)
    
    out = f(vf, [0.0, 2.0], hmlie([RN([1.0, -8.0])]), vf_event)

    assert abs(out.t[0] - 0.0) <= 1.0e-14
    assert abs(out.t[-1] - (2*10-1)/16.0) <= 1e-6

    for ii in range(len(out.t)):
        assert abs(out.y[ii, 0] - t_to_y1(out.t[ii])) <= 1e-6


def test_Flow_B1_vector():
    """
    Problem B-1
    
    Uses VectorXd representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """

    from lielab.topos import Flow
    import numpy as np

    f = Flow()

    K = 5.0/3.0
    L = 1.0/4.0

    def vf(t, y):
        E = y[0] - y[1]

        dy0 = y[0]
        dy1 = 0
        if (E < -L/K):
            dy1 = -L
        elif ((-L/K <= E) and (E <= L/K)):
            dy1 = K*E
        elif (L/K < E):
            dy1 = L

        dy = np.array([dy0, dy1])
        return dy
    
    
    t1 = 0.1569
    y2t1 = 1.0199
    
    def t_to_y1(t):
        return np.exp(t)

    def t_to_y2(t):
        if (t < t1):
            return K/(K+1)*np.exp(t) + 1/(K+1)*np.exp(-K*t)
        return L*(t - t1) + y2t1
    
    out = f(vf, [0.0, 0.5], [1.0, 1.0])

    assert abs(out.t[0] - 0.0) <= 1.0e-14
    assert abs(out.t[-1] - 0.5) <= 1.0e-14

    for ii in range(len(out.t)):
        assert abs(out.y[ii, 0] - t_to_y1(out.t[ii])) <= 1e-7
        assert abs(out.y[ii, 1] - t_to_y2(out.t[ii])) <= 1e-4 # Answer given in the document is only good to 1e-4 (see t1 and y2(t1))


def test_Flow_B1_hom():
    """
    Problem B-1
    
    Uses hom representation and events.
     
    Source: Thompson, S. A collection of test problems for ordinary differential
    equation solvers which have provisions for rootfinding. No. ORNL/TM-9912.
    Oak Ridge National Lab., TN (USA), 1987.
    """

    from lielab.domain import halie, hmlie, rn, RN, gl
    from lielab.topos import Flow, CustomMuntheKaas
    from lielab.functions import exp
    import numpy as np

    f = Flow()
    f.stepper = CustomMuntheKaas()

    def left(g, G):
        _g0 = g.space[0]
        _g1 = g.space[1]
        _Y0 = G.space[0]
        _Y1 = G.space[1]

        Y0next = RN(np.dot(exp(_g0)._data, _Y0._data))
        Y1next = exp(_g1)*_Y1
        return hmlie([Y0next, Y1next])
    
    f.stepper.left = left

    K = 5.0/3.0
    L = 1.0/4.0

    def vf(t, M):
        y0 = M.space[0](0)
        y1 = M.space[1](0)

        E = y0 - y1

        if (E < -L/K):
            dy1 = -L
        elif ((-L/K <= E) and (E <= L/K)):
            dy1 = K*E
        elif (L/K < E):
            dy1 = L
        else:
            dy1 = 0

        return halie([gl.basis(0,1), dy1*rn.basis(0,2)])
    
    t1 = 0.1569
    y2t1 = 1.0199
    
    def t_to_y1(t):
        return np.exp(t)

    def t_to_y2(t):
        if (t < t1):
            return K/(K+1)*np.exp(t) + 1/(K+1)*np.exp(-K*t)
        return L*(t - t1) + y2t1
    
    out = f(vf, [0.0, 0.5], hmlie([RN([1.0]), RN([1.0])]))

    assert abs(out.t[0] - 0.0) <= 1.0e-14
    assert abs(out.t[-1] - 0.5) <= 1.0e-14

    for ii in range(len(out.t)):
        assert abs(out.y[ii, 0] - t_to_y1(out.t[ii])) <= 1e-14 # Answer is analytic
        assert abs(out.y[ii, 1] - t_to_y2(out.t[ii])) <= 1e-4 # Answer given in the document is only good to 1e-4 (see t1 and y2(t1))
