from lielab.testing import *

def test_from_SO3():
    """
    Tests the from_SO3 function
    TODO: update to new SU serialization
    """

    from lielab.domain import so, SO, SU
    from lielab.functions import exp

    # Build 90 degree rotations in x, y, and z
    rx = exp(np.pi/2.0*so.basis(0,3))
    ry = exp(np.pi/2.0*so.basis(1,3))
    rz = exp(np.pi/2.0*so.basis(2,3))

    # Test 90 degree x rotation
    _qx = SU.from_SO3(rx)
    qx = _qx.to_quaternion()

    assert abs(qx[0] - np.sqrt(2)/2.0) <= TOL_FINE
    assert abs(qx[1] - np.sqrt(2)/2.0) <= TOL_FINE
    assert abs(qx[2] - 0.0) <= TOL_FINE
    assert abs(qx[3] - 0.0) <= TOL_FINE

    # Test 90 degree y rotation
    _qy = SU.from_SO3(ry)
    qy = _qy.to_quaternion()

    assert abs(qy[0] - np.sqrt(2)/2.0) <= TOL_FINE
    assert abs(qy[1] - 0.0) <= TOL_FINE
    assert abs(qy[2] - np.sqrt(2)/2.0) <= TOL_FINE
    assert abs(qy[3] - 0.0) <= TOL_FINE

    # Test 90 degree z rotation
    _qz = SU.from_SO3(rz)
    qz = _qz.to_quaternion()

    assert abs(qz[0] - np.sqrt(2)/2.0) <= TOL_FINE
    assert abs(qz[1] - 0.0) <= TOL_FINE
    assert abs(qz[2] - 0.0) <= TOL_FINE
    assert abs(qz[3] - np.sqrt(2)/2.0) <= TOL_FINE


def test_from_quaternion():
    """
    Tests quaternions against well-known identities.
    """

    from lielab.domain import CompositeGroup, SU
    from lielab.testing import check_almost_equal_nulp

    qm1 = SU.from_quaternion(-1, 0, 0, 0)
    qi = SU.from_quaternion(0, 1, 0, 0)
    qj = SU.from_quaternion(0, 0, 1, 0)
    qk = SU.from_quaternion(0, 0, 0, 1)

    # Hamilton's identities
    # i^2 = j^2 = k^2 = -1
    assert check_almost_equal_nulp(CompositeGroup([qi*qi]), CompositeGroup([qm1]), 1, True)
    assert check_almost_equal_nulp(CompositeGroup([qj*qj]), CompositeGroup([qm1]), 1, True)
    assert check_almost_equal_nulp(CompositeGroup([qk*qk]), CompositeGroup([qm1]), 1, True)

    # ij = -ji = -k
    assert check_almost_equal_nulp(CompositeGroup([qi*qj]), CompositeGroup([(qj*qi).inverse()]), 1, True)
    assert check_almost_equal_nulp(CompositeGroup([qi*qj]), CompositeGroup([qk.inverse()]), 1, True)

    # jk = -kj = -i
    assert check_almost_equal_nulp(CompositeGroup([qj*qk]), CompositeGroup([(qk*qj).inverse()]), 1, True)
    assert check_almost_equal_nulp(CompositeGroup([qj*qk]), CompositeGroup([qi.inverse()]), 1, True)

    # ki = -ik = -j
    assert check_almost_equal_nulp(CompositeGroup([qk*qi]), CompositeGroup([(qi*qk).inverse()]), 1, True)
    assert check_almost_equal_nulp(CompositeGroup([qk*qi]), CompositeGroup([qj.inverse()]), 1, True)
