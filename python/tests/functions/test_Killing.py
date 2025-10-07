import pytest
from lielab.testing import *
import numpy as np

def test_Killing():
    """
    Tests the Killing function.
    """

    from lielab.domain import so
    from lielab.functions import Killing

    rx = so.from_vector([1,0,0])
    ry = so.from_vector([0,1,0])
    rz = so.from_vector([0,0,1])

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

    rx = so.from_vector([1.0, 0.0, 0.0])

    K = Killingform(rx)
    Id = np.identity(rx.get_dimension())

    assert abs(K.trace() + 6) <= TOL_FINE
    assert_matrix(K @ np.linalg.inv(K), Id)

    so6 = so.basis(0,6)

    K = Killingform(so6)
    Id = np.identity(so6.get_dimension())

    assert abs(K.trace() + 120) <= TOL_FINE
    assert_matrix(K @ np.linalg.inv(K), Id)
