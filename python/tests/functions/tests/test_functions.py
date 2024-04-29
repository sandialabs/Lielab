"""
Test cases for the functions Python bindings.
"""

import pytest
from lielab.testing import *
import numpy as np

def test_Ad():
    """
    Tests the Ad function.
    """

    from lielab.domain import so, SO
    from lielab.functions import exp, Ad

    u = so(3)
    v = so(3)
    w = so(3)
    ansso = so(3)
    Gso = SO(3)

    u.set_vector([1,0,0])
    v.set_vector([0,1,0])
    w.set_vector([0,0,1])
    Gso = exp(v)

    # GuG^-1
    ansso = Ad(Gso, u)
    truthso = np.array([[0, 0.841470984807896, 0],
                        [-0.841470984807897, 0, -0.540302305868140],
                        [0, 0.540302305868140, 0]])
    
    assert_matrix(ansso.get_matrix(), truthso)

    # GvG^-1 = v when G = exp(v)
    ansso = Ad(Gso, v)
    
    assert_matrix(ansso.get_matrix(), v.get_matrix())

    # GwG^-1
    ansso = Ad(Gso, w)
    truthso = np.array([[0, -0.540302305868140, 0],
                        [0.540302305868140, 0, -0.841470984807897],
                        [0, 0.841470984807897, 0]])
    
    assert_matrix(ansso.get_matrix(), truthso)


def test_cayley1():
    """
    Tests the cayley1 function
    """

    from lielab.domain import so
    from lielab.functions import cayley1

    rx = so([1,0,0])
    ry = so([0,1,0])

    # Values calculated by hand
    ex1 = cayley1(rx)
    ans = np.array([[1.0, 0.0, 0.0],
                    [0.0, 0.6,-0.8],
                    [0.0, 0.8, 0.6]])

    assert_matrix(ex1.get_matrix(), ans)

    ex2 = cayley1(ry)
    ans = np.array([[0.6, 0.0, 0.8],
                    [0.0, 1.0, 0.0],
                    [-0.8, 0.0, 0.6]])

    assert_matrix(ex2.get_matrix(), ans)


def test_cayley2():
    """
    Tests the cayley2 function
    """

    from lielab.domain import so
    from lielab.functions import cayley2

    rx = so([1,0,0])
    ry = so([0,1,0])

    # Values calculated by hand
    ex1 = cayley2(rx)
    ans = np.array([[1.0, 0.0, 0.0],
                    [0.0, 0.6,-0.8],
                    [0.0, 0.8, 0.6]])

    assert_matrix(ex1.get_matrix(), ans)

    ex2 = cayley2(rx + 2*ry)
    ans = np.array([[0.0, 0.0, 1.0],
                    [0.8, 0.6, 0.0],
                    [-0.6, 0.8, 0.0]])

    assert_matrix(ex2.get_matrix(), ans)

def test_cayley1_and_2():
    """
    Tests cayley1 and cayley2 together with known identities.
    """

    from lielab.domain import so
    from lielab.functions import cayley1, cayley2

    # Identity cayley1 = cayley2 for all basis elements
    dim = so.basis(0, 10).get_dimension()
    for ii in range(dim):
        g = so.basis(ii, 10)
        assert_domain(cayley1(g), cayley2(g))


def test_Killing():
    """
    Tests the Killing function.
    """

    from lielab.domain import so
    from lielab.functions import Killing

    rx = so([1,0,0])
    ry = so([0,1,0])
    rz = so([0,0,1])

    assert abs(Killing(rx,rx) + 2.0) <= TOL_FINE
    assert abs(Killing(rx,ry) + 0.0) <= TOL_FINE
    assert abs(Killing(rx,rz) + 0.0) <= TOL_FINE
    assert abs(Killing(ry,rx) + 0.0) <= TOL_FINE
    assert abs(Killing(ry,ry) + 2.0) <= TOL_FINE
    assert abs(Killing(ry,rz) + 0.0) <= TOL_FINE
    assert abs(Killing(rz,rx) + 0.0) <= TOL_FINE
    assert abs(Killing(rz,ry) + 0.0) <= TOL_FINE
    assert abs(Killing(rz,rz) + 2.0) <= TOL_FINE

def test_Killingform():
    """
    Tests the Killingform function.
    """

    from lielab.domain import so
    from lielab.functions import Killingform

    rx = so([1.0, 0.0, 0.0])

    K = Killingform(rx)
    Id = np.identity(rx.get_dimension())

    assert abs(K.trace() + 6) <= TOL_FINE
    assert_matrix(K @ np.linalg.inv(K), Id)

    so6 = so.basis(0,6)

    K = Killingform(so6)
    Id = np.identity(so6.get_dimension())

    assert abs(K.trace() + 120) <= TOL_FINE
    assert_matrix(K @ np.linalg.inv(K), Id)


def test_dcayley1inv():
    """
    Tests the inverse of the dcayley1 function
    """
    
    from lielab.domain import so
    from lielab.functions import dcayley1inv

    u = so(3)
    v = so(3)
    ansso = so(3)

    u.set_vector([1,0,0])
    v.set_vector([0,1,0])

    ansso = dcayley1inv(u, v)
    truthso  = np.array([[0.0, 0.5, 1.0],
                        [-0.5, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    ansso = dcayley1inv(v, u)
    truthso = np.array([[0.0,-0.5, 0.0],
                        [0.5, 0.0,-1.0],
                        [0.0, 1.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)


def test_dexp():
    """
    Tests the dexp function.
    """

    from lielab.domain import so, rn
    from lielab.functions import dexp

    u = so(3)
    v = so(3)
    ansso = so(3)

    xx = np.array([1.0, 0.0, 0.0])
    yy = np.array([0.0, 1.0, 0.0])

    u.set_vector(xx)
    v.set_vector(yy)


    # order = 0
    ansso = dexp(u, v, 0)
    truthso = np.array([[0.0, 0.0, 1.0],
                        [0.0, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])
    
    assert_matrix(ansso.get_matrix(), truthso)

    # order = 1
    ansso = dexp(u, v, 1)
    truthso = np.array([[0.0, -0.5, 1.0],
                        [0.5, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 2
    ansso = dexp(u, v, 2)
    truthso = np.array([[0.0, -0.5, 0.833333333333333],
                        [0.5, 0.0, 0.0],
                        [-0.833333333333333, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 3
    ansso = dexp(u, v, 3)
    truthso = np.array([[0.0, -0.458333333333333, 0.833333333333333],
                        [0.458333333333333, 0.0, 0.0],
                        [-0.833333333333333, 0.0, 0.0]])

    # order = 4
    ansso = dexp(u, v, 4)
    truthso = np.array([[0.0, -0.458333333333333, 0.841666666666667],
                        [0.458333333333333, 0.0, 0.0],
                        [-0.841666666666667, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 8
    ansso = dexp(u, v, 8)
    truthso = np.array([[0.0, -0.459697420634921, 0.841471009700176],
                        [0.459697420634921, 0.0, 0.0],
                        [-0.841471009700176, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    x = rn(4)
    y = rn(4)
    ansrn = rn(4)

    x.set_vector(xx)
    y.set_vector(yy)

    # default order
    ansrn = dexp(x, y)
    truthrn = np.array([[0, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]])

    assert_matrix(ansrn.get_matrix(), truthrn)

    # ridiculous order (checks abelian speedhack)
    ansrn = dexp(x, y, 999999999)
    truthrn = np.array([[0, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]])

    assert_matrix(ansrn.get_matrix(), truthrn)


def test_dexpinv():
    """
    Tests the dexpinv function.
    """

    from lielab.domain import so, rn
    from lielab.functions import dexpinv

    u = so(3)
    v = so(3)
    ansso = so(3)

    xx = np.array([1.0, 0.0, 0.0])
    yy = np.array([0.0, 1.0, 0.0])

    u.set_vector(xx)
    v.set_vector(yy)


    # order = 0
    ansso = dexpinv(u, v, 0)
    truthso = np.array([[0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0]])
    
    assert_matrix(ansso.get_matrix(), truthso)

    # order = 1
    ansso = dexpinv(u, v, 1)
    truthso = np.array([[0.0, 0.0, 1.0],
                        [0.0, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 2
    ansso = dexpinv(u, v, 2)
    truthso = np.array([[0.0, 0.5, 1.0],
                        [-0.5, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 4
    ansso = dexpinv(u, v, 4)
    truthso = np.array([[0.0, 0.5, 0.916666666666667],
                        [-0.5, 0.0, 0.0],
                        [-0.916666666666667, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 12
    ansso = dexpinv(u, v, 12)
    truthso = np.array([[0.0, 0.500000000000000, 0.915243861398375],
                        [-0.500000000000000, 0.0, 0.0],
                        [-0.915243861398375, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    x = rn(4)
    y = rn(4)
    ansrn = rn(4)

    x.set_vector(xx)
    y.set_vector(yy)

    # default order
    ansrn = dexpinv(x, y)
    truthrn = np.array([[0, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]])

    assert_matrix(ansrn.get_matrix(), truthrn)

    # ridiculous order (checks abelian speedhack)
    ansrn = dexpinv(x, y, 999999999)
    truthrn = np.array([[0, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]])

    assert_matrix(ansrn.get_matrix(), truthrn)
