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

    from lielab.domain import SU

    qm1 = SU.from_quaternion(-1, 0, 0, 0)
    qi = SU.from_quaternion(0, 1, 0, 0)
    qj = SU.from_quaternion(0, 0, 1, 0)
    qk = SU.from_quaternion(0, 0, 0, 1)

    # Hamilton's identities
    # i^2 = j^2 = k^2 = -1
    assert_matrix(qi*qi, qm1)
    assert_matrix(qj*qj, qm1)
    assert_matrix(qk*qk, qm1)

    # ij = -ji = -k
    assert_matrix(qi*qj, (qj*qi).inverse())
    assert_matrix(qi*qj, qk.inverse())

    # jk = -kj = -i
    assert_matrix(qj*qk, (qk*qj).inverse())
    assert_matrix(qj*qk, qi.inverse())

    # ki = -ik = -j
    assert_matrix(qk*qi, (qi*qk).inverse())
    assert_matrix(qk*qi, qj.inverse())
