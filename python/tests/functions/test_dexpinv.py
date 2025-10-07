import pytest
from lielab.testing import *
import numpy as np

def test_dexpinv_numerical():
    """
    Tests the dexpinv_numerical function.
    """

    from lielab.domain import so, rn
    from lielab.functions import dexpinv_numerical

    u = so(3)
    v = so(3)
    ansso = so(3)

    xx = np.array([1.0, 0.0, 0.0])
    yy = np.array([0.0, 1.0, 0.0])

    u.set_vector(xx)
    v.set_vector(yy)


    # order = 0
    ansso = dexpinv_numerical(u, v, 0)
    truthso = np.array([[0.0, 0.0, 1.0],
                        [0.0, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])
    
    assert_matrix(ansso.get_matrix(), truthso)

    # order = 1
    ansso = dexpinv_numerical(u, v, 1)
    truthso = np.array([[0.0, 0.5, 1.0],
                        [-0.5, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 3
    ansso = dexpinv_numerical(u, v, 3)
    truthso = np.array([[0.0, 0.5, 0.916666666666667],
                        [-0.5, 0.0, 0.0],
                        [-0.916666666666667, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 11
    ansso = dexpinv_numerical(u, v, 11)
    truthso = np.array([[0.0, 0.500000000000000, 0.915243861398375],
                        [-0.500000000000000, 0.0, 0.0],
                        [-0.915243861398375, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    x = rn(3)
    y = rn(3)
    ansrn = rn(3)

    x.set_vector(xx)
    y.set_vector(yy)

    # default order
    ansrn = dexpinv_numerical(x, y)
    truthrn = np.array([[0, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]])

    assert_matrix(ansrn.get_matrix(), truthrn)

    # ridiculous order (checks abelian speedhack)
    ansrn = dexpinv_numerical(x, y, 999999999)
    truthrn = np.array([[0, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]])

    assert_matrix(ansrn.get_matrix(), truthrn)
