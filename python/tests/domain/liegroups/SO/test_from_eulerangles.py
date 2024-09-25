# from math import sqrt, np.pi
from lielab.testing import *
import lielab
# from lielab.domain.SO import *
# from lielab.domain import so, SO
# from lielab.functions import exp

DCMId = lielab.domain.SO(3)
DCMrotx = lielab.functions.exp(np.pi/2.0*lielab.domain.so.basis(0,3))
DCMroty = lielab.functions.exp(np.pi/2.0*lielab.domain.so.basis(1,3))
DCMrotz = lielab.functions.exp(np.pi/2.0*lielab.domain.so.basis(2,3))
some_angle = np.pi/2.0*5.0/7.0

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
