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

# def test_dcm_to_quaternion():
#     """
#     Tests the dcm_to_quaternion_function
#     TODO: update to new SU serialization
#     """

#     from lielab.domain import so, SO
#     from lielab.functions import exp
#     from lielab.transform import dcm_to_quaternion

#     # Build 90 degree rotations in x, y, and z
#     rx = exp(np.pi/2.0*so.basis(0,3))
#     ry = exp(np.pi/2.0*so.basis(1,3))
#     rz = exp(np.pi/2.0*so.basis(2,3))

#     # Test 90 degree x rotation
#     _qx = dcm_to_quaternion(rx)
#     qx = _qx.serialize()

#     assert abs(qx[0] - sqrt(2)/2.0) <= TOL_FINE
#     assert abs(qx[1] - sqrt(2)/2.0) <= TOL_FINE
#     assert abs(qx[2] - 0.0) <= TOL_FINE
#     assert abs(qx[3] - 0.0) <= TOL_FINE

#     # Test 90 degree y rotation
#     _qy = dcm_to_quaternion(ry)
#     qy = _qy.serialize()

#     assert abs(qy[0] - sqrt(2)/2.0) <= TOL_FINE
#     assert abs(qy[1] - 0.0) <= TOL_FINE
#     assert abs(qy[2] - sqrt(2)/2.0) <= TOL_FINE
#     assert abs(qy[3] - 0.0) <= TOL_FINE

#     # Test 90 degree z rotation
#     _qz = dcm_to_quaternion(rz)
#     qz = _qz.serialize()

#     assert abs(qz[0] - sqrt(2)/2.0) <= TOL_FINE
#     assert abs(qz[1] - 0.0) <= TOL_FINE
#     assert abs(qz[2] - 0.0) <= TOL_FINE
#     assert abs(qz[3] - sqrt(2)/2.0) <= TOL_FINE


def test_dcm_to_eanglebody123():
    """
    Tests the dcm_to_eanglebody123 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody123

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
    theta1, theta2, theta3 = dcm_to_eanglebody123(SO.from_eulerangles_body123(0.5, np.pi/2.0, 0.5));
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody123(SO.from_eulerangles_body123(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody231():
    """
    Tests the dcm_to_eanglebody231 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody231

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
    theta1, theta2, theta3 = dcm_to_eanglebody231(SO.from_eulerangles_body231(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody231(SO.from_eulerangles_body231(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody312():
    """
    Tests the dcm_to_eanglebody312 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody312

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
    theta1, theta2, theta3 = dcm_to_eanglebody312(SO.from_eulerangles_body312(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody312(SO.from_eulerangles_body312(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody132():
    """
    Tests the dcm_to_eanglebody132 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody132

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
    theta1, theta2, theta3 = dcm_to_eanglebody132(SO.from_eulerangles_body132(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody132(SO.from_eulerangles_body132(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody213():
    """
    Tests the dcm_to_eanglebody213 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody213

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
    theta1, theta2, theta3 = dcm_to_eanglebody213(SO.from_eulerangles_body213(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody213(SO.from_eulerangles_body213(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody321():
    """
    Tests the dcm_to_eanglebody321 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody321

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
    theta1, theta2, theta3 = dcm_to_eanglebody321(SO.from_eulerangles_body321(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody321(SO.from_eulerangles_body321(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody121():
    """
    Tests the dcm_to_eanglebody121 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody121

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
    theta1, theta2, theta3 = dcm_to_eanglebody121(SO.from_eulerangles_body121(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody121(SO.from_eulerangles_body121(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody131():
    """
    Tests the dcm_to_eanglebody131 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody131

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
    theta1, theta2, theta3 = dcm_to_eanglebody131(SO.from_eulerangles_body131(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody131(SO.from_eulerangles_body131(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody212():
    """
    Tests the dcm_to_eanglebody212 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody212

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
    theta1, theta2, theta3 = dcm_to_eanglebody212(SO.from_eulerangles_body212(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody212(SO.from_eulerangles_body212(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody232():
    """
    Tests the dcm_to_eanglebody232 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody232

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
    theta1, theta2, theta3 = dcm_to_eanglebody232(SO.from_eulerangles_body232(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody232(SO.from_eulerangles_body232(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody313():
    """
    Tests the dcm_to_eanglebody313 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody313

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
    theta1, theta2, theta3 = dcm_to_eanglebody313(SO.from_eulerangles_body313(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody313(SO.from_eulerangles_body313(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglebody323():
    """
    Tests the dcm_to_eanglebody323 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglebody323

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
    theta1, theta2, theta3 = dcm_to_eanglebody323(SO.from_eulerangles_body323(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglebody323(SO.from_eulerangles_body323(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace123():
    """
    Tests the dcm_to_eanglespace123 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace123

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
    theta1, theta2, theta3 = dcm_to_eanglespace123(SO.from_eulerangles_space123(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace123(SO.from_eulerangles_space123(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace231():
    """
    Tests the dcm_to_eanglespace231 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace231

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
    theta1, theta2, theta3 = dcm_to_eanglespace231(SO.from_eulerangles_space231(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace231(SO.from_eulerangles_space231(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace312():
    """
    Tests the dcm_to_eanglespace312 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace312

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
    theta1, theta2, theta3 = dcm_to_eanglespace312(SO.from_eulerangles_space312(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 0.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace312(SO.from_eulerangles_space312(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace132():
    """
    Tests the dcm_to_eanglespace132 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace132

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
    theta1, theta2, theta3 = dcm_to_eanglespace132(SO.from_eulerangles_space132(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace132(SO.from_eulerangles_space132(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace213():
    """
    Tests the dcm_to_eanglespace213 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace213

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
    theta1, theta2, theta3 = dcm_to_eanglespace213(SO.from_eulerangles_space213(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace213(SO.from_eulerangles_space213(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace321():
    """
    Tests the dcm_to_eanglespace321 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace321

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
    theta1, theta2, theta3 = dcm_to_eanglespace321(SO.from_eulerangles_space321(0.5, np.pi/2.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - np.pi/2.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace321(SO.from_eulerangles_space321(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace121():
    """
    Tests the dcm_to_eanglespace121 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace121

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
    theta1, theta2, theta3 = dcm_to_eanglespace121(SO.from_eulerangles_space121(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace121(SO.from_eulerangles_space121(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace131():
    """
    Tests the dcm_to_eanglespace131 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace131

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
    theta1, theta2, theta3 = dcm_to_eanglespace131(SO.from_eulerangles_space131(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace131(SO.from_eulerangles_space131(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace212():
    """
    Tests the dcm_to_eanglespace212 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace212

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
    theta1, theta2, theta3 = dcm_to_eanglespace212(SO.from_eulerangles_space212(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace212(SO.from_eulerangles_space212(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace232():
    """
    Tests the dcm_to_eanglespace232 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace232

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
    theta1, theta2, theta3 = dcm_to_eanglespace232(SO.from_eulerangles_space232(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace232(SO.from_eulerangles_space232(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace313():
    """
    Tests the dcm_to_eanglespace313 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace313

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
    theta1, theta2, theta3 = dcm_to_eanglespace313(SO.from_eulerangles_space313(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace313(SO.from_eulerangles_space313(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


def test_dcm_to_eanglespace323():
    """
    Tests the dcm_to_eanglespace323 function.
    """

    from lielab.domain import SO
    from lielab.transform import dcm_to_eanglespace323

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
    theta1, theta2, theta3 = dcm_to_eanglespace323(SO.from_eulerangles_space323(0.5, 0.0, 0.5))
    assert abs(theta1 - 1.0) <= TOL_FINE
    assert abs(theta2 - 0.0) <= TOL_FINE
    assert theta3 == 0.0

    # Rotate by a random angle (also asserts the inverse function)
    theta1, theta2, theta3 = dcm_to_eanglespace323(SO.from_eulerangles_space323(some_angle, some_angle, some_angle))
    assert theta1 == some_angle
    assert theta2 == some_angle
    assert theta3 == some_angle


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
