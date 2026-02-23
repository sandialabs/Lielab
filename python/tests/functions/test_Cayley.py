import pytest
from lielab.testing import *
import numpy as np

def test_cay():
    """
    Tests the cay function
    """

    from lielab.domain import so
    from lielab.functions import cay
    from lielab.testing import check_almost_equal_tol

    rx = so.from_vector([1,0,0])
    ry = so.from_vector([0,1,0])

    # Values calculated by hand
    ex1 = cay(rx)
    ans = np.array([[1.0, 0.0, 0.0],
                    [0.0, 0.6,-0.8],
                    [0.0, 0.8, 0.6]])

    assert check_almost_equal_tol(ex1.get_matrix(), ans)

    ex2 = cay(ry)
    ans = np.array([[0.6, 0.0, 0.8],
                    [0.0, 1.0, 0.0],
                    [-0.8, 0.0, 0.6]])

    assert check_almost_equal_tol(ex2.get_matrix(), ans)

def test_cay2():
    """
    Tests the cay2 function
    """

    from lielab.domain import so
    from lielab.functions import cay2
    from lielab.testing import check_almost_equal_tol

    rx = so.from_vector([1,0,0])
    ry = so.from_vector([0,1,0])

    # Values calculated by hand
    ex1 = cay2(rx)
    ans = np.array([[1.0, 0.0, 0.0],
                    [0.0, 0.6,-0.8],
                    [0.0, 0.8, 0.6]])

    assert check_almost_equal_tol(ex1.get_matrix(), ans)

    ex2 = cay2(rx + 2*ry)
    ans = np.array([[0.0, 0.0, 1.0],
                    [0.8, 0.6, 0.0],
                    [-0.6, 0.8, 0.0]])

    assert check_almost_equal_tol(ex2.get_matrix(), ans)

def test_cay_and_cay2():
    """
    Tests cay and cay2 together with known identities.
    """

    from lielab.domain import CompositeGroup, so
    from lielab.functions import cay, cay2
    from lielab.testing import check_almost_equal_nulp

    # Identity Cayley = Cayley2 for all basis elements
    dim = so.basis(0, 10).get_dimension()
    for ii in range(dim):
        g = so.basis(ii, 10)
        assert check_almost_equal_nulp(CompositeGroup([cay(g)]), CompositeGroup([cay2(g)]), 1, True)

def test_dcayinv():
    """
    Tests the inverse of the dcay function
    """
    
    from lielab.domain import so
    from lielab.functions import dcayinv
    from lielab.testing import check_almost_equal_tol

    u = so(3)
    v = so(3)
    ansso = so(3)

    u.set_vector([1,0,0])
    v.set_vector([0,1,0])

    ansso = dcayinv(u, v)
    truthso  = np.array([[0.0, 0.5, 1.0],
                        [-0.5, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])

    assert check_almost_equal_tol(ansso.get_matrix(), truthso)

    ansso = dcayinv(v, u)
    truthso = np.array([[0.0,-0.5, 0.0],
                        [0.5, 0.0,-1.0],
                        [0.0, 1.0, 0.0]])

    assert check_almost_equal_tol(ansso.get_matrix(), truthso)
