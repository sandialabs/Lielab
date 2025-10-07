import pytest
from lielab.testing import *
import numpy as np

def test_log_SO2():
    from lielab.domain import so, SO
    from lielab.functions import log, exp

    # Identity and small valued tests.
    eps = 1.0e-15

    SO2Id = SO(2)
    zero = log(SO2Id)
    assert (zero.get_shape() == 2)
    zerobar = zero.get_vector()
    assert (zerobar.size == 1)
    assert (zerobar[0] == 0.0)

    SO2_epsx = exp(eps*so.basis(0, 2))
    epsx = log(SO2_epsx)
    assert (epsx.get_shape() == 2)
    epsxbar = epsx.get_vector()
    assert (epsxbar.size == 1)
    assert (epsxbar[0] == eps)

    # Large-ish valued tests.
    m = 2.0

    SO2_mx = exp(m*so.basis(0, 2))
    mx = log(SO2_mx)
    assert (mx.get_shape() == 2)
    mxbar = mx.get_vector()
    assert (mxbar.size == 1)
    assert (mxbar[0] == m)

def test_log_SO3():
    from lielab.domain import so, SO
    from lielab.functions import log, exp

    # Identity and small valued tests. Tests the divide by zero is correctly handled.
    eps = 1.0e-15

    SO3Id = SO(3)
    zero = log(SO3Id)
    assert (zero.get_shape() == 3)
    zerobar = zero.get_vector()
    assert (zerobar.size == 3)
    assert (zerobar[0] == 0.0)
    assert (zerobar[1] == 0.0)
    assert (zerobar[2] == 0.0)

    SO3_epsx = exp(eps*so.basis(0, 3))
    epsx = log(SO3_epsx)
    assert (epsx.get_shape() == 3)
    epsxbar = epsx.get_vector()
    assert (epsxbar.size == 3)
    assert (epsxbar[0] == eps)
    assert (epsxbar[1] == 0.0)
    assert (epsxbar[2] == 0.0)

    SO3_epsy = exp(eps*so.basis(1, 3))
    epsy = log(SO3_epsy)
    assert (epsy.get_shape() == 3)
    epsybar = epsy.get_vector()
    assert (epsybar.size == 3)
    assert (epsybar[0] == 0.0)
    assert (epsybar[1] == eps)
    assert (epsybar[2] == 0.0)

    SO3_epsz = exp(eps*so.basis(2, 3))
    epsz = log(SO3_epsz)
    assert (epsz.get_shape() == 3)
    epszbar = epsz.get_vector()
    assert (epszbar.size == 3)
    assert (epszbar[0] == 0.0)
    assert (epszbar[1] == 0.0)
    assert (epszbar[2] == eps)

    SO3_epsv = exp(eps*so.basis(0, 3) + 2*eps*so.basis(1, 3) - 1*eps*so.basis(2, 3))
    epsv = log(SO3_epsv)
    assert (epsv.get_shape() == 3)
    epsvbar = epsv.get_vector()
    assert (epsvbar.size == 3)
    assert (epsvbar[0] == eps)
    assert (epsvbar[1] == 2*eps)
    assert (epsvbar[2] == -eps)

    # Large-ish valued tests. Tests angle > pi/2 is correctly handled.
    m = 2.0

    SO3_mx = exp(m*so.basis(0, 3))
    mx = log(SO3_mx)
    assert (mx.get_shape() == 3)
    mxbar = mx.get_vector()
    assert (mxbar.size == 3)
    assert (mxbar[0] == m)
    assert (mxbar[1] == 0.0)
    assert (mxbar[2] == 0.0)

    SO3_my = exp(m*so.basis(1, 3))
    my = log(SO3_my)
    assert (my.get_shape() == 3)
    mybar = my.get_vector()
    assert (mybar.size == 3)
    assert (mybar[0] == 0.0)
    assert (mybar[1] == m)
    assert (mybar[2] == 0.0)

    SO3_mz = exp(m*so.basis(2, 3))
    mz = log(SO3_mz)
    assert (mz.get_shape() == 3)
    mzbar = mz.get_vector()
    assert (mzbar.size == 3)
    assert (mzbar[0] == 0.0)
    assert (mzbar[1] == 0.0)
    assert (mzbar[2] == m)

    SO3_mv = exp(m/3.0*so.basis(0, 3) + 2*m/3.0*so.basis(1, 3) - 1*m/3.0*so.basis(2, 3))
    mv = log(SO3_mv)
    assert (mv.get_shape() == 3)
    mvbar = mv.get_vector()
    assert (mvbar.size == 3)
    assert (np.abs(mvbar[0] - m/3.0) < 1e-14)
    assert (np.abs(mvbar[1] - 2*m/3.0) < 1e-14)
    assert (np.abs(mvbar[2] - -m/3.0) < 1e-14)

    # Gimbal lock points

    SO3_pix = exp(np.pi*so.basis(0, 3))
    pix = log(SO3_pix)
    assert (pix.get_shape() == 3)
    pixbar = pix.get_vector()
    assert (pixbar.size == 3)
    assert (np.abs(pixbar[0] - np.pi) < 1e-14)
    assert (pixbar[1] == 0.0)
    assert (pixbar[2] == 0.0)

    SO3_piy = exp(np.pi*so.basis(1, 3))
    piy = log(SO3_piy)
    assert (piy.get_shape() == 3)
    piybar = piy.get_vector()
    assert (piybar.size == 3)
    assert (piybar[0] == 0.0)
    assert (np.abs(piybar[1] - np.pi) < 1e-14)
    assert (piybar[2] == 0.0)

    SO3_piz = exp(np.pi*so.basis(2, 3))
    piz = log(SO3_piz)
    assert (piz.get_shape() == 3)
    pizbar = piz.get_vector()
    assert (pizbar.size == 3)
    assert (pizbar[0] == 0.0)
    assert (pizbar[1] == 0.0)
    assert (np.abs(pizbar[2] - np.pi) < 1e-14)
