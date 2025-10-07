# from math import sqrt, np.pi
from lielab.testing import *
import lielab

DCMId = lielab.domain.SO(3)
DCMrotx = lielab.functions.exp(np.pi/2.0*lielab.domain.so.basis(0,3))
DCMroty = lielab.functions.exp(np.pi/2.0*lielab.domain.so.basis(1,3))
DCMrotz = lielab.functions.exp(np.pi/2.0*lielab.domain.so.basis(2,3))
some_angle = np.pi/2.0*5.0/7.0

def test_SO():
    """
    Tests SO against well-known identities.
    """

    from lielab.domain import so
    from lielab.functions import exp

    x = exp(so.basis(0,3))
    y = exp(so.basis(1,3))
    z = exp(so.basis(2,3))

    # Group product
    x*y
    assert_domain(x*y, (y.inverse()*x.inverse()).inverse())

    # Group inverse
    x.inverse()

def test_from_SU2():
    """
    Tests the from_SU2 function.
    """

    from lielab.domain import su, SU, SO
    from lielab.functions import exp

    # Build 90 degree rotations in x, y, and z
    # TODO: simplify with basis()
    u = su(2)
    v = su(2)
    w = su(2)
    u.set_vector([0, 0, 1/2])
    v.set_vector([0, 1/2, 0])
    w.set_vector([1/2, 0, 0])

    qx = exp(np.pi/2.0*u)
    qy = exp(np.pi/2.0*v)
    qz = exp(np.pi/2.0*w)

    # Test 90 degree x rotation
    # rx = quaternion_to_dcm(qx)
    rx = SO.from_SU2(qx)

    assert_domain(rx, DCMrotx)

    # Test 90 degree y rotation
    ry = SO.from_SU2(qy)

    assert_domain(ry, DCMroty)

    # Test 90 degree z rotation
    rz = SO.from_SU2(qz)

    assert_domain(rz, DCMrotz)


def test_from_eulerangles_body123():
    """
    Tests the from_eulerangles_body123 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body123(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body123(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body123(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body123(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body123(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body231():
    """
    Tests the from_eulerangles_body231 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body231(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body231(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body231(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body231(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body231(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body312():
    """
    Tests the from_eulerangles_body312 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body312(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body312(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body312(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body312(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body312(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body132():
    """
    Tests the from_eulerangles_body132 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body132(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body132(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body132(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body132(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body132(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body213():
    """
    Tests the from_eulerangles_body213 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body213(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body213(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body213(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body213(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body213(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body321():
    """
    Tests the from_eulerangles_body321 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body321(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body321(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body321(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body321(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body321(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body121():
    """
    Tests the from_eulerangles_body121 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body121(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body121(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body121(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body121(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body121(np.pi/2.0, np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body121(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body131():
    """
    Tests the from_eulerangles_body131 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body131(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body131(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body131(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body131(np.pi/2.0, -np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body131(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body131(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body212():
    """
    Tests the from_eulerangles_body212 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body212(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body212(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body212(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body212(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body212(np.pi/2.0, -np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body212(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body232():
    """
    Tests the from_eulerangles_body232 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body232(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body232(np.pi/2.0, np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body232(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body232(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body232(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body232(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body313():
    """
    Tests the from_eulerangles_body313 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body313(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body313(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body313(-np.pi/2.0, -np.pi/2.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body313(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body313(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body313(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_body323():
    """
    Tests the from_eulerangles_body323 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_body323(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_body323(np.pi/2.0, -np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_body323(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body323(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_body323(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_body323(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space123():
    """
    Tests the from_eulerangles_space123 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space123(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space123(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space123(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space123(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space123(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space231():
    """
    Tests the from_eulerangles_space231 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space231(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space231(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space231(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space231(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space231(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space312():
    """
    Tests the from_eulerangles_space312 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space312(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space312(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space312(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space312(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space312(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space132():
    """
    Tests the from_eulerangles_space132 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space132(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space132(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space132(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space132(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space132(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space213():
    """
    Tests the from_eulerangles_space213 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space213(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space213(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space213(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space213(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space213(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space321():
    """
    Tests the from_eulerangles_space321 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space321(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space321(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space321(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space321(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space321(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space121():
    """
    Tests the from_eulerangles_space121 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space121(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space121(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space121(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space121(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space121(np.pi/2.0, -np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space121(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space131():
    """
    Tests the from_eulerangles_space131 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space131(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space131(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space131(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space131(np.pi/2.0, np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space131(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space131(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space212():
    """
    Tests the from_eulerangles_space212 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space212(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space212(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space212(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space212(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space212(np.pi/2.0, np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space212(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space232():
    """
    Tests the from_eulerangles_space232 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space232(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space232(np.pi/2.0, -np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space232(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space232(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space232(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space232(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space313():
    """
    Tests the from_eulerangles_space313 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space313(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space313(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space313(-np.pi/2.0, np.pi/2.0, np.pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space313(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space313(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space313(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_from_eulerangles_space323():
    """
    Tests the from_eulerangles_space323 function.
    """

    from lielab.domain import SO

    # Identity
    dcm = SO.from_eulerangles_space323(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = SO.from_eulerangles_space323(np.pi/2.0, np.pi/2.0, -np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = SO.from_eulerangles_space323(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space323(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate 90 degrees by z-axis
    dcm = SO.from_eulerangles_space323(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = SO.from_eulerangles_space323(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_matrix()) - 1) < TOL_FINE
    assert abs((dcm.get_matrix() @ dcm.get_matrix().transpose()).trace() - 3) < TOL_FINE


def test_to_eulerangles_body123():
    """
    Tests the to_eulerangles_body123 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body123()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body123()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body123()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body123()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body123(0.5, np.pi/2.0, 0.5).to_eulerangles_body123()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body123(some_angle, some_angle, some_angle).to_eulerangles_body123()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body231():
    """
    Tests the to_eulerangles_body231 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body231()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body231()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body231()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body231()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body231(0.5, np.pi/2.0, 0.5).to_eulerangles_body231()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body231(some_angle, some_angle, some_angle).to_eulerangles_body231()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body312():
    """
    Tests the to_eulerangles_body312 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body312()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body312()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body312()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body312()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body312(0.5, np.pi/2.0, 0.5).to_eulerangles_body312()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body312(some_angle, some_angle, some_angle).to_eulerangles_body312()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body132():
    """
    Tests the to_eulerangles_body132 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body132()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body132()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body132()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body132()
    assert abs(theta1- 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3- 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body132(0.5, np.pi/2.0, 0.5).to_eulerangles_body132()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body132(some_angle, some_angle, some_angle).to_eulerangles_body132()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body213():
    """
    Tests the to_eulerangles_body213 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body213()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body213()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body213()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body213()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body213(0.5, np.pi/2.0, 0.5).to_eulerangles_body213()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body213(some_angle, some_angle, some_angle).to_eulerangles_body213()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body321():
    """
    Tests the to_eulerangles_body321 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body321()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body321()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body321()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body321()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body321(0.5, np.pi/2.0, 0.5).to_eulerangles_body321()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body321(some_angle, some_angle, some_angle).to_eulerangles_body321()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body121():
    """
    Tests the to_eulerangles_body121 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body121()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body121()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body121()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body121()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body121(0.5, 0.0, 0.5).to_eulerangles_body121()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body121(some_angle, some_angle, some_angle).to_eulerangles_body121()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body131():
    """
    Tests the to_eulerangles_body131 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body131()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body131()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body131()
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body131()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body131(0.5, 0.0, 0.5).to_eulerangles_body131()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body131(some_angle, some_angle, some_angle).to_eulerangles_body131()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body212():
    """
    Tests the to_eulerangles_body212 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body212()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body212()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body212()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body212()
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body212(0.5, 0.0, 0.5).to_eulerangles_body212()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body212(some_angle, some_angle, some_angle).to_eulerangles_body212()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body232():
    """
    Tests the to_eulerangles_body232 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body232()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body232()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body232()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body232()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body232(0.5, 0.0, 0.5).to_eulerangles_body232()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body232(some_angle, some_angle, some_angle).to_eulerangles_body232()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body313():
    """
    Tests the to_eulerangles_body313 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body313()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body313()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body313()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body313()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body313(0.5, 0.0, 0.5).to_eulerangles_body313()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body313(some_angle, some_angle, some_angle).to_eulerangles_body313()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_body323():
    """
    Tests the to_eulerangles_body323 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_body323()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_body323()
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_body323()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_body323()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_body323(0.5, 0.0, 0.5).to_eulerangles_body323()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_body323(some_angle, some_angle, some_angle).to_eulerangles_body323()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space123():
    """
    Tests the to_eulerangles_space123 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space123()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space123()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space123()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space123()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space123(0.5, np.pi/2.0, 0.5).to_eulerangles_space123()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space123(some_angle, some_angle, some_angle).to_eulerangles_space123()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space231():
    """
    Tests the to_eulerangles_space231 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space231()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space231()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space231()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space231()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space231(0.5, np.pi/2.0, 0.5).to_eulerangles_space231()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space231(some_angle, some_angle, some_angle).to_eulerangles_space231()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space312():
    """
    Tests the to_eulerangles_space312 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space312()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space312()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space312()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space312()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space312(0.5, np.pi/2.0, 0.5).to_eulerangles_space312()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space312(some_angle, some_angle, some_angle).to_eulerangles_space312()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space132():
    """
    Tests the to_eulerangles_space132 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space132()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space132()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space132()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space132()
    assert abs(theta1- 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3- 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space132(0.5, np.pi/2.0, 0.5).to_eulerangles_space132()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space132(some_angle, some_angle, some_angle).to_eulerangles_space132()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space213():
    """
    Tests the to_eulerangles_space213 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space213()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space213()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space213()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space213()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space213(0.5, np.pi/2.0, 0.5).to_eulerangles_space213()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space213(some_angle, some_angle, some_angle).to_eulerangles_space213()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space321():
    """
    Tests the to_eulerangles_space321 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space321()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space321()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space321()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space321()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space321(0.5, np.pi/2.0, 0.5).to_eulerangles_space321()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space321(some_angle, some_angle, some_angle).to_eulerangles_space321()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space121():
    """
    Tests the to_eulerangles_space121 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space121()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space121()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space121()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space121()
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space121(0.5, 0.0, 0.5).to_eulerangles_space121()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space121(some_angle, some_angle, some_angle).to_eulerangles_space121()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space131():
    """
    Tests the to_eulerangles_space131 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space131()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space131()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space131()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space131()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space131(0.5, 0.0, 0.5).to_eulerangles_space131()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space131(some_angle, some_angle, some_angle).to_eulerangles_space131()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space212():
    """
    Tests the to_eulerangles_space212 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space212()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space212()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space212()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space212()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space212(0.5, 0.0, 0.5).to_eulerangles_space212()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space212(some_angle, some_angle, some_angle).to_eulerangles_space212()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space232():
    """
    Tests the to_eulerangles_space232 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space232()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space232()
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space232()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space232()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space232(0.5, 0.0, 0.5).to_eulerangles_space232()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space232(some_angle, some_angle, some_angle).to_eulerangles_space232()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space313():
    """
    Tests the to_eulerangles_space313 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space313()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space313()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space313()
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space313()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space313(0.5, 0.0, 0.5).to_eulerangles_space313()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space313(some_angle, some_angle, some_angle).to_eulerangles_space313()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_to_eulerangles_space323():
    """
    Tests the to_eulerangles_space323 function.
    """

    from lielab.domain import SO

    # Identity
    theta1, theta2, theta3 = DCMId.to_eulerangles_space323()
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotx.to_eulerangles_space323()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = DCMroty.to_eulerangles_space323()
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = DCMrotz.to_eulerangles_space323()
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lock condition
    theta1, theta2, theta3 = SO.from_eulerangles_space323(0.5, 0.0, 0.5).to_eulerangles_space323()
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = SO.from_eulerangles_space323(some_angle, some_angle, some_angle).to_eulerangles_space323()
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle

