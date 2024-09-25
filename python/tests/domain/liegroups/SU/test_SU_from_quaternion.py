from lielab.testing import *

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
