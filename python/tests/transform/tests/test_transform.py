from math import sqrt, pi
from lielab.testing import *
from lielab.transform import *
from lielab.domain import so, SO
from lielab.functions import exp

DCMId = SO(3)
DCMrotx = exp(np.pi/2.0*so.basis(0,3))
DCMroty = exp(np.pi/2.0*so.basis(1,3))
DCMrotz = exp(np.pi/2.0*so.basis(2,3))
some_angle = np.pi/2.0*5.0/7.0

def test_dcm_to_quaternion():
    """
    Tests the dcm_to_quaternion_function
    """

    from lielab.domain import so, SO
    from lielab.functions import exp
    from lielab.transform import dcm_to_quaternion

    # Build 90 degree rotations in x, y, and z
    rx = exp(np.pi/2.0*so.basis(0,3))
    ry = exp(np.pi/2.0*so.basis(1,3))
    rz = exp(np.pi/2.0*so.basis(2,3))

    # Test 90 degree x rotation
    _qx = dcm_to_quaternion(rx)
    qx = _qx.serialize()

    assert abs(qx[0] - sqrt(2)/2.0) <= TOL_FINE
    assert abs(qx[1] - sqrt(2)/2.0) <= TOL_FINE
    assert abs(qx[2] - 0.0) <= TOL_FINE
    assert abs(qx[3] - 0.0) <= TOL_FINE

    # Test 90 degree y rotation
    _qy = dcm_to_quaternion(ry)
    qy = _qy.serialize()

    assert abs(qy[0] - sqrt(2)/2.0) <= TOL_FINE
    assert abs(qy[1] - 0.0) <= TOL_FINE
    assert abs(qy[2] - sqrt(2)/2.0) <= TOL_FINE
    assert abs(qy[3] - 0.0) <= TOL_FINE

    # Test 90 degree z rotation
    _qz = dcm_to_quaternion(rz)
    qz = _qz.serialize()

    assert abs(qz[0] - sqrt(2)/2.0) <= TOL_FINE
    assert abs(qz[1] - 0.0) <= TOL_FINE
    assert abs(qz[2] - 0.0) <= TOL_FINE
    assert abs(qz[3] - sqrt(2)/2.0) <= TOL_FINE


def test_dcm_to_eanglebody123():
    """
    Tests the dcm_to_eanglebody123 function.
    """

    from lielab.transform import dcm_to_eanglebody123, eanglebody123_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody123(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody123(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody123(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody123(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody123(eanglebody123_to_dcm(0.5, np.pi/2.0, 0.5));
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody123(eanglebody123_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody231():
    """
    Tests the dcm_to_eanglebody231 function.
    """

    from lielab.transform import dcm_to_eanglebody231, eanglebody231_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody231(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody231(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody231(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody231(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody231(eanglebody231_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody231(eanglebody231_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody312():
    """
    Tests the dcm_to_eanglebody312 function.
    """

    from lielab.transform import dcm_to_eanglebody312, eanglebody312_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody312(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody312(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody312(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody312(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody312(eanglebody312_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody312(eanglebody312_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody132():
    """
    Tests the dcm_to_eanglebody132 function.
    """

    from lielab.transform import dcm_to_eanglebody132, eanglebody132_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody132(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody132(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody132(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody132(DCMrotz)
    assert abs(theta1- 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3- 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody132(eanglebody132_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody132(eanglebody132_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody213():
    """
    Tests the dcm_to_eanglebody213 function.
    """

    from lielab.transform import dcm_to_eanglebody213, eanglebody213_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody213(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody213(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody213(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody213(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody213(eanglebody213_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody213(eanglebody213_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody321():
    """
    Tests the dcm_to_eanglebody321 function.
    """

    from lielab.transform import dcm_to_eanglebody321, eanglebody321_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody321(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody321(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody321(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody321(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody321(eanglebody321_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody321(eanglebody321_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody121():
    """
    Tests the dcm_to_eanglebody121 function.
    """

    from lielab.transform import dcm_to_eanglebody121, eanglebody121_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody121(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody121(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody121(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody121(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody121(eanglebody121_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody121(eanglebody121_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody131():
    """
    Tests the dcm_to_eanglebody131 function.
    """

    from lielab.transform import dcm_to_eanglebody131, eanglebody131_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody131(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody131(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody131(DCMroty)
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody131(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody131(eanglebody131_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody131(eanglebody131_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody212():
    """
    Tests the dcm_to_eanglebody212 function.
    """

    from lielab.transform import dcm_to_eanglebody212, eanglebody212_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody212(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody212(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody212(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody212(DCMrotz)
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody212(eanglebody212_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody212(eanglebody212_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody232():
    """
    Tests the dcm_to_eanglebody232 function.
    """

    from lielab.transform import dcm_to_eanglebody232, eanglebody232_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody232(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody232(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody232(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody232(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody232(eanglebody232_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody232(eanglebody232_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody313():
    """
    Tests the dcm_to_eanglebody313 function.
    """

    from lielab.transform import dcm_to_eanglebody313, eanglebody313_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody313(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody313(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody313(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody313(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody313(eanglebody313_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody313(eanglebody313_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody323():
    """
    Tests the dcm_to_eanglebody323 function.
    """

    from lielab.transform import dcm_to_eanglebody323, eanglebody323_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglebody323(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody323(DCMrotx)
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglebody323(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglebody323(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglebody323(eanglebody323_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody323(eanglebody323_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace123():
    """
    Tests the dcm_to_eanglespace123 function.
    """

    from lielab.transform import dcm_to_eanglespace123, eanglespace123_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace123(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace123(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace123(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace123(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace123(eanglespace123_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace123(eanglespace123_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace231():
    """
    Tests the dcm_to_eanglespace231 function.
    """

    from lielab.transform import dcm_to_eanglespace231, eanglespace231_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace231(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace231(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace231(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace231(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace231(eanglespace231_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace231(eanglespace231_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace312():
    """
    Tests the dcm_to_eanglespace312 function.
    """

    from lielab.transform import dcm_to_eanglespace312, eanglespace312_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace312(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace312(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace312(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace312(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace312(eanglespace312_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace312(eanglespace312_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace132():
    """
    Tests the dcm_to_eanglespace132 function.
    """

    from lielab.transform import dcm_to_eanglespace132, eanglespace132_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace132(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace132(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace132(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace132(DCMrotz)
    assert abs(theta1- 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3- 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace132(eanglespace132_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace132(eanglespace132_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace213():
    """
    Tests the dcm_to_eanglespace213 function.
    """

    from lielab.transform import dcm_to_eanglespace213, eanglespace213_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace213(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace213(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace213(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace213(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace213(eanglespace213_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace213(eanglespace213_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace321():
    """
    Tests the dcm_to_eanglespace321 function.
    """

    from lielab.transform import dcm_to_eanglespace321, eanglespace321_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace321(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace321(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace321(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace321(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace321(eanglespace321_to_dcm(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace321(eanglespace321_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace121():
    """
    Tests the dcm_to_eanglespace121 function.
    """

    from lielab.transform import dcm_to_eanglespace121, eanglespace121_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace121(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace121(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace121(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace121(DCMrotz)
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace121(eanglespace121_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace121(eanglespace121_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace131():
    """
    Tests the dcm_to_eanglespace131 function.
    """

    from lielab.transform import dcm_to_eanglespace131, eanglespace131_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace131(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace131(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace131(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace131(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace131(eanglespace131_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace131(eanglespace131_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace212():
    """
    Tests the dcm_to_eanglespace212 function.
    """

    from lielab.transform import dcm_to_eanglespace212, eanglespace212_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace212(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace212(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace212(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace212(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace212(eanglespace212_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace212(eanglespace212_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace232():
    """
    Tests the dcm_to_eanglespace232 function.
    """

    from lielab.transform import dcm_to_eanglespace232, eanglespace232_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace232(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace232(DCMrotx)
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace232(DCMroty)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace232(DCMrotz)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace232(eanglespace232_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace232(eanglespace232_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace313():
    """
    Tests the dcm_to_eanglespace313 function.
    """

    from lielab.transform import dcm_to_eanglespace313, eanglespace313_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace313(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace313(DCMrotx)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace313(DCMroty)
    assert abs(theta1 + np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace313(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace313(eanglespace313_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace313(eanglespace313_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace323():
    """
    Tests the dcm_to_eanglespace323 function.
    """

    from lielab.transform import dcm_to_eanglespace323, eanglespace323_to_dcm

    # Identity
    theta1, theta2, theta3 = dcm_to_eanglespace323(DCMId)
    assert theta1 == 0.0
    assert theta2 == 0.0
    assert theta3 == 0.0

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace323(DCMrotx)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 + np.pi/2.0) <= TOL_FINE

    # Rotate 90 degrees by y-axis
    theta1, theta2, theta3 = dcm_to_eanglespace323(DCMroty)
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Rotate 90 degrees by x-axis
    theta1, theta2, theta3 = dcm_to_eanglespace323(DCMrotz)
    assert abs(theta1 - np.pi/2.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert abs(theta3 - 0.0) <= TOL_FINE

    # Check gimbal lockc condition
    theta1, theta2, theta3 = dcm_to_eanglespace323(eanglespace323_to_dcm(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace323(eanglespace323_to_dcm(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_eanglebody123_to_dcm():
    """
    Tests the eanglebody123_to_dcm function.
    """

    from lielab.transform import eanglebody123_to_dcm

    # Identity
    dcm = eanglebody123_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody123_to_dcm(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody123_to_dcm(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody123_to_dcm(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody123_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody231_to_dcm():
    """
    Tests the eanglebody231_to_dcm function.
    """

    from lielab.transform import eanglebody231_to_dcm

    # Identity
    dcm = eanglebody231_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody231_to_dcm(0.0, 0.0, np.pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody231_to_dcm(np.pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody231_to_dcm(0.0, np.pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody231_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody312_to_dcm():
    """
    Tests the eanglebody312_to_dcm function.
    """

    from lielab.transform import eanglebody312_to_dcm

    # Identity
    dcm = eanglebody312_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody312_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody312_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody312_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody312_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody132_to_dcm():
    """
    Tests the eanglebody132_to_dcm function.
    """

    from lielab.transform import eanglebody132_to_dcm

    # Identity
    dcm = eanglebody132_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody132_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody132_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody132_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody132_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody213_to_dcm():
    """
    Tests the eanglebody213_to_dcm function.
    """

    from lielab.transform import eanglebody213_to_dcm

    # Identity
    dcm = eanglebody213_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody213_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody213_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody213_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody213_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody321_to_dcm():
    """
    Tests the eanglebody321_to_dcm function.
    """

    from lielab.transform import eanglebody321_to_dcm

    # Identity
    dcm = eanglebody321_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody321_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody321_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody321_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody321_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody121_to_dcm():
    """
    Tests the eanglebody121_to_dcm function.
    """

    from lielab.transform import eanglebody121_to_dcm

    # Identity
    dcm = eanglebody121_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody121_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody121_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody121_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody121_to_dcm(pi/2.0, pi/2.0, -pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody121_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody131_to_dcm():
    """
    Tests the eanglebody131_to_dcm function.
    """

    from lielab.transform import eanglebody131_to_dcm

    # Identity
    dcm = eanglebody131_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody131_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody131_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody131_to_dcm(pi/2.0, -pi/2.0, -pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody131_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody131_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody212_to_dcm():
    """
    Tests the eanglebody212_to_dcm function.
    """

    from lielab.transform import eanglebody212_to_dcm

    # Identity
    dcm = eanglebody212_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody212_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody212_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody212_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody212_to_dcm(pi/2.0, -pi/2.0, -pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody212_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody232_to_dcm():
    """
    Tests the eanglebody232_to_dcm function.
    """

    from lielab.transform import eanglebody232_to_dcm

    # Identity
    dcm = eanglebody232_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody232_to_dcm(pi/2.0, pi/2.0, -pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody232_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody232_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody232_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody232_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody313_to_dcm():
    """
    Tests the eanglebody313_to_dcm function.
    """

    from lielab.transform import eanglebody313_to_dcm

    # Identity
    dcm = eanglebody313_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody313_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody313_to_dcm(-pi/2.0, -pi/2.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody313_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody313_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody313_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglebody323_to_dcm():
    """
    Tests the eanglebody323_to_dcm function.
    """

    from lielab.transform import eanglebody323_to_dcm

    # Identity
    dcm = eanglebody323_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglebody323_to_dcm(pi/2.0, -pi/2.0, -pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglebody323_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody323_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglebody323_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglebody323_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace123_to_dcm():
    """
    Tests the eanglespace123_to_dcm function.
    """

    from lielab.transform import eanglespace123_to_dcm

    # Identity
    dcm = eanglespace123_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace123_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace123_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace123_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace123_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace231_to_dcm():
    """
    Tests the eanglespace231_to_dcm function.
    """

    from lielab.transform import eanglespace231_to_dcm

    # Identity
    dcm = eanglespace231_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace231_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace231_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace231_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace231_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace312_to_dcm():
    """
    Tests the eanglespace312_to_dcm function.
    """

    from lielab.transform import eanglespace312_to_dcm

    # Identity
    dcm = eanglespace312_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace312_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace312_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace312_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace312_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace132_to_dcm():
    """
    Tests the eanglespace132_to_dcm function.
    """

    from lielab.transform import eanglespace132_to_dcm

    # Identity
    dcm = eanglespace132_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace132_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace132_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace132_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace132_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace213_to_dcm():
    """
    Tests the eanglespace213_to_dcm function.
    """

    from lielab.transform import eanglespace213_to_dcm

    # Identity
    dcm = eanglespace213_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace213_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace213_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace213_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace213_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace321_to_dcm():
    """
    Tests the eanglespace321_to_dcm function.
    """

    from lielab.transform import eanglespace321_to_dcm

    # Identity
    dcm = eanglespace321_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace321_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace321_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace321_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace321_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace121_to_dcm():
    """
    Tests the eanglespace121_to_dcm function.
    """

    from lielab.transform import eanglespace121_to_dcm

    # Identity
    dcm = eanglespace121_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace121_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace121_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace121_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace121_to_dcm(pi/2.0, -pi/2.0, -pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace121_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace131_to_dcm():
    """
    Tests the eanglespace131_to_dcm function.
    """

    from lielab.transform import eanglespace131_to_dcm

    # Identity
    dcm = eanglespace131_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace131_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace131_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace131_to_dcm(pi/2.0, pi/2.0, -pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace131_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace131_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace212_to_dcm():
    """
    Tests the eanglespace212_to_dcm function.
    """

    from lielab.transform import eanglespace212_to_dcm

    # Identity
    dcm = eanglespace212_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace212_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace212_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace212_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace212_to_dcm(pi/2.0, pi/2.0, -pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace212_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace232_to_dcm():
    """
    Tests the eanglespace232_to_dcm function.
    """

    from lielab.transform import eanglespace232_to_dcm

    # Identity
    dcm = eanglespace232_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace232_to_dcm(pi/2.0, -pi/2.0, -pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace232_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace232_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace232_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace232_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace313_to_dcm():
    """
    Tests the eanglespace313_to_dcm function.
    """

    from lielab.transform import eanglespace313_to_dcm

    # Identity
    dcm = eanglespace313_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace313_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace313_to_dcm(-pi/2.0, pi/2.0, pi/2.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace313_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace313_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace313_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_eanglespace323_to_dcm():
    """
    Tests the eanglespace323_to_dcm function.
    """

    from lielab.transform import eanglespace323_to_dcm

    # Identity
    dcm = eanglespace323_to_dcm(0.0, 0.0, 0.0)
    assert_domain(DCMId, dcm)

    # Rotate 90 degrees by x-axis
    dcm = eanglespace323_to_dcm(pi/2.0, pi/2.0, -pi/2.0)
    assert_domain(DCMrotx, dcm)

    # Rotate 90 degrees by y-axis
    dcm = eanglespace323_to_dcm(0.0, pi/2.0, 0.0)
    assert_domain(DCMroty, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace323_to_dcm(pi/2.0, 0.0, 0.0)
    assert_domain(DCMrotz, dcm)

    # Rotate 90 degrees by z-axis
    dcm = eanglespace323_to_dcm(0.0, 0.0, pi/2.0)
    assert_domain(DCMrotz, dcm)

    # Rotate by a random angle
    dcm = eanglespace323_to_dcm(some_angle, some_angle, some_angle)
    assert abs(np.linalg.det(dcm.get_ados_representation()) - 1) < TOL_FINE
    assert abs((dcm.get_ados_representation() @ dcm.get_ados_representation().transpose()).trace() - 3) < TOL_FINE


def test_quaternion_to_dcm():
    """
    Tests the quaternion_to_dcm function.
    """

    from lielab.domain import su, SU
    from lielab.transform import quaternion_to_dcm
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
    rx = quaternion_to_dcm(qx)

    assert_domain(rx, DCMrotx)

    # Test 90 degree y rotation
    ry = quaternion_to_dcm(qy)

    assert_domain(ry, DCMroty)

    # Test 90 degree z rotation
    rz = quaternion_to_dcm(qz)

    assert_domain(rz, DCMrotz)
