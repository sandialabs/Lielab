import pytest
from lielab.testing import *
import numpy as np

def test_dexp_numerical():
    """
    Tests the dexp_numerical function.
    """

    from lielab.domain import so, rn
    from lielab.functions import dexp_numerical

    u = so(3)
    v = so(3)
    ansso = so(3)

    xx = np.array([1.0, 0.0, 0.0])
    yy = np.array([0.0, 1.0, 0.0])

    u.set_vector(xx)
    v.set_vector(yy)


    # order = 0
    ansso = dexp_numerical(u, v, 0)
    truthso = np.array([[0.0, 0.0, 1.0],
                        [0.0, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])
    
    assert_matrix(ansso.get_matrix(), truthso)

    # order = 1
    ansso = dexp_numerical(u, v, 1)
    truthso = np.array([[0.0, -0.5, 1.0],
                        [0.5, 0.0, 0.0],
                        [-1.0, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 2
    ansso = dexp_numerical(u, v, 2)
    truthso = np.array([[0.0, -0.5, 0.833333333333333],
                        [0.5, 0.0, 0.0],
                        [-0.833333333333333, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 3
    ansso = dexp_numerical(u, v, 3)
    truthso = np.array([[0.0, -0.458333333333333, 0.833333333333333],
                        [0.458333333333333, 0.0, 0.0],
                        [-0.833333333333333, 0.0, 0.0]])

    # order = 4
    ansso = dexp_numerical(u, v, 4)
    truthso = np.array([[0.0, -0.458333333333333, 0.841666666666667],
                        [0.458333333333333, 0.0, 0.0],
                        [-0.841666666666667, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    # order = 8
    ansso = dexp_numerical(u, v, 8)
    truthso = np.array([[0.0, -0.459697420634921, 0.841471009700176],
                        [0.459697420634921, 0.0, 0.0],
                        [-0.841471009700176, 0.0, 0.0]])

    assert_matrix(ansso.get_matrix(), truthso)

    x = rn(3)
    y = rn(3)
    ansrn = rn(3)

    x.set_vector(xx)
    y.set_vector(yy)

    # default order
    ansrn = dexp_numerical(x, y)
    truthrn = np.array([[0, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]])

    assert_matrix(ansrn.get_matrix(), truthrn)

    # ridiculous order (checks abelian speedhack)
    ansrn = dexp_numerical(x, y, 999999999)
    truthrn = np.array([[0, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]])

    assert_matrix(ansrn.get_matrix(), truthrn)
