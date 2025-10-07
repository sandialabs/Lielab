import pytest
from lielab.testing import *
import numpy as np

def test_exp_so2():
    from lielab.domain import so, SO
    from lielab.functions import exp

    # Identity
    so2zero = so(2)
    SO2Id = exp(so2zero)
    assert (SO2Id.get_shape() == 2)
    assert (SO2Id(0, 0) == 1.0)
    assert (SO2Id(0, 1) == 0.0)
    assert (SO2Id(1, 0) == 0.0)
    assert (SO2Id(1, 1) == 1.0)

    # Other
    SO2_mx = exp(np.pi/4.0*so.basis(0, 2))
    assert (SO2_mx.get_shape() == 2)
    assert (np.abs(SO2_mx(0, 0) - np.sqrt(2.0)/2.0) < 1e-14)
    assert (np.abs(SO2_mx(0, 1) - -np.sqrt(2.0)/2.0) < 1e-14)
    assert (np.abs(SO2_mx(1, 0) - np.sqrt(2.0)/2.0) < 1e-14)
    assert (np.abs(SO2_mx(1, 1) - np.sqrt(2.0)/2.0) < 1e-14)

def test_exp_so3():
    from lielab.domain import so, SO
    from lielab.functions import exp

    # Identity and small valued tests. Tests the divide by zero is correctly handled.
    eps = 1.0e-15

    so3_zero = so(3)
    SO3_Id = exp(so3_zero)
    assert (SO3_Id.get_shape() == 3)
    assert (SO3_Id(0, 0) == 1.0)
    assert (SO3_Id(0, 1) == 0.0)
    assert (SO3_Id(0, 2) == 0.0)
    assert (SO3_Id(1, 0) == 0.0)
    assert (SO3_Id(1, 1) == 1.0)
    assert (SO3_Id(1, 2) == 0.0)
    assert (SO3_Id(2, 0) == 0.0)
    assert (SO3_Id(2, 1) == 0.0)
    assert (SO3_Id(2, 2) == 1.0)

    so3_epsx = eps*so.basis(0, 3)
    SO3_epsx = exp(so3_epsx)
    assert (SO3_epsx.get_shape() == 3)
    assert (SO3_epsx(0, 0) == 1.0)
    assert (SO3_epsx(0, 1) == 0.0)
    assert (SO3_epsx(0, 2) == 0.0)
    assert (SO3_epsx(1, 0) == 0.0)
    assert (SO3_epsx(1, 1) == 1.0)
    assert (SO3_epsx(1, 2) == -eps)
    assert (SO3_epsx(2, 0) == 0.0)
    assert (SO3_epsx(2, 1) == eps)
    assert (SO3_epsx(2, 2) == 1.0)

    so3_epsy = eps*so.basis(1, 3)
    SO3_epsy = exp(so3_epsy)
    assert (SO3_epsy.get_shape() == 3)
    assert (SO3_epsy(0, 0) == 1.0)
    assert (SO3_epsy(0, 1) == 0.0)
    assert (SO3_epsy(0, 2) == eps)
    assert (SO3_epsy(1, 0) == 0.0)
    assert (SO3_epsy(1, 1) == 1.0)
    assert (SO3_epsy(1, 2) == 0.0)
    assert (SO3_epsy(2, 0) == -eps)
    assert (SO3_epsy(2, 1) == 0.0)
    assert (SO3_epsy(2, 2) == 1.0)

    so3_epsz = eps*so.basis(2, 3)
    SO3_epsz = exp(so3_epsz)
    assert (SO3_epsz.get_shape() == 3)
    assert (SO3_epsz(0, 0) == 1.0)
    assert (SO3_epsz(0, 1) == -eps)
    assert (SO3_epsz(0, 2) == 0.0)
    assert (SO3_epsz(1, 0) == eps)
    assert (SO3_epsz(1, 1) == 1.0)
    assert (SO3_epsz(1, 2) == 0.0)
    assert (SO3_epsz(2, 0) == 0.0)
    assert (SO3_epsz(2, 1) == 0.0)
    assert (SO3_epsz(2, 2) == 1.0)

    # Other
    m = np.pi/4.0

    so3_mx = m*so.basis(0, 3)
    SO3_mx = exp(so3_mx)
    assert (SO3_mx.get_shape() == 3)
    assert (SO3_mx(0, 0) == 1.0)
    assert (SO3_mx(0, 1) == 0.0)
    assert (SO3_mx(0, 2) == 0.0)
    assert (SO3_mx(1, 0) == 0.0)
    assert (np.abs(SO3_mx(1, 1) - np.sqrt(2.0)/2.0) < 1e-14)
    assert (np.abs(SO3_mx(1, 2) - -np.sqrt(2.0)/2.0) < 1e-14)
    assert (SO3_mx(2, 0) == 0.0)
    assert (np.abs(SO3_mx(2, 1) - np.sqrt(2.0)/2.0) < 1e-14)
    assert (np.abs(SO3_mx(2, 2) - np.sqrt(2.0)/2.0) < 1e-14)
