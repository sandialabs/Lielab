import pytest
from lielab.testing import *
import numpy as np

def test_Cayley():
    """
    Tests the Cayley function
    """

    from lielab.domain import so
    from lielab.functions import Cayley

    rx = so.from_vector([1,0,0])
    ry = so.from_vector([0,1,0])

    # Values calculated by hand
    ex1 = Cayley(rx)
    ans = np.array([[1.0, 0.0, 0.0],
                    [0.0, 0.6,-0.8],
                    [0.0, 0.8, 0.6]])

    assert_matrix(ex1.get_matrix(), ans)

    ex2 = Cayley(ry)
    ans = np.array([[0.6, 0.0, 0.8],
                    [0.0, 1.0, 0.0],
                    [-0.8, 0.0, 0.6]])

    assert_matrix(ex2.get_matrix(), ans)

def test_Cayley2():
    """
    Tests the Cayley2 function
    """

    from lielab.domain import so
    from lielab.functions import Cayley2

    rx = so.from_vector([1,0,0])
    ry = so.from_vector([0,1,0])

    # Values calculated by hand
    ex1 = Cayley2(rx)
    ans = np.array([[1.0, 0.0, 0.0],
                    [0.0, 0.6,-0.8],
                    [0.0, 0.8, 0.6]])

    assert_matrix(ex1.get_matrix(), ans)

    ex2 = Cayley2(rx + 2*ry)
    ans = np.array([[0.0, 0.0, 1.0],
                    [0.8, 0.6, 0.0],
                    [-0.6, 0.8, 0.0]])

    assert_matrix(ex2.get_matrix(), ans)

def test_Cayley_and_Cayley2():
    """
    Tests Cayley and Cayley2 together with known identities.
    """

    from lielab.domain import so
    from lielab.functions import Cayley, Cayley2

    # Identity Cayley = Cayley2 for all basis elements
    dim = so.basis(0, 10).get_dimension()
    for ii in range(dim):
        g = so.basis(ii, 10)
        assert_domain(Cayley(g), Cayley2(g))

def test_dCayleyinv():
    """
    Tests the inverse of the dCayley function
    """
    
    from lielab.domain import so
    from lielab.functions import dCayleyinv

    u = so(3)
    v = so(3)
    ansso = so(3)

    u.set_vector([1,0,0])
    v.set_vector([0,1,0])

    ansso = dCayleyinv(u, v)
    truthso  = np.array([[0.0, 0.5, 1.0],
                        [-0.5, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    ansso = dCayleyinv(v, u)
    truthso = np.array([[0.0,-0.5, 0.0],
                        [0.5, 0.0,-1.0],
                        [0.0, 1.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)
