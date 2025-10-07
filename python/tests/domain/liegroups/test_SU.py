from lielab.testing import *
import pytest

def complex(a,b):
    return a + b*1j

def test_SU_to_string():
    from lielab.domain import SU

    x0 = SU(0)
    assert (x0.to_string() == "SU(0)")
    x1 = SU(1)
    assert (x1.to_string() == "SU(1)")
    x10 = SU(10)
    assert (x10.to_string() == "SU(10)")

def test_SU_main_initializer():
    from lielab.domain import SU

    xblank = SU()
    assert (xblank.get_dimension() == 0)

    x0 = SU(0)
    assert (x0.get_dimension() == 0)
    x1 = SU(1)
    assert (x1.get_dimension() == 0)
    x10 = SU(2)
    assert (x10.get_dimension() == 3)

def test_SU_matrix_initializer():
    from lielab.domain import SU

    x0 = SU(np.random.rand(0, 0) + 1j*np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = SU(np.random.rand(1, 1) + 1j*np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = SU(np.random.rand(2, 2) + 1j*np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        SU(np.random.rand(2, 3) + 1j*np.random.rand(2, 3))
    
    with pytest.raises(RuntimeError):
        SU(np.random.rand(3, 2) + 1j*np.random.rand(3, 2))

def test_SU_from_shape_initializer():
    from lielab.domain import SU

    x0 = SU.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = SU.from_shape(1)
    assert (x1.get_dimension() == 0)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = SU.from_shape(2)
    assert (x2.get_dimension() == 3)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_SU_get_dimension():
    from lielab.domain import SU

    zero = SU(0)
    one = SU(1)
    two = SU(2)
    three = SU(3)
    four = SU(4)
    five = SU(5)
    six = SU(6)
    seven = SU(7)
    eight = SU(8)

    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 0)
    assert (two.get_dimension() == 3)
    assert (three.get_dimension() == 8)
    assert (four.get_dimension() == 15)
    assert (five.get_dimension() == 24)
    assert (six.get_dimension() == 35)
    assert (seven.get_dimension() == 48)
    assert (eight.get_dimension() == 63)

def test_SU_serialize_unserialize():
    """
    * Tests the serialize/unserialize operation.
    """

    from lielab.domain import SU

    x0 = SU(0)
    x0.unserialize([])
    x0bar = x0.serialize()

    assert (x0bar.size == 0)

    x1 = SU(1)
    x1.unserialize([1.0, 2.0, 3.0])
    x1bar = x1.serialize()

    assert (x1bar.size == 2)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 2.0)

    x1.unserialize([4.0, 5.0])
    x1bar = x1.serialize()

    assert (x1bar.size == 2)
    assert (x1bar[0] == 4.0)
    assert (x1bar[1] == 5.0)

    x1.unserialize([6.0])
    x1bar = x1.serialize()

    assert (x1bar.size == 2)
    assert (x1bar[0] == 6.0)
    assert (x1bar[1] == 5.0)

    x2 = SU(2)
    x2.unserialize([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 8)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 3.0)
    assert (x2bar[3] == 4.0)
    assert (x2bar[4] == 5.0)
    assert (x2bar[5] == 6.0)
    assert (x2bar[6] == 7.0)
    assert (x2bar[7] == 0.0)

    x2.unserialize([8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 8)
    assert (x2bar[0] == 8.0)
    assert (x2bar[1] == 9.0)
    assert (x2bar[2] == 10.0)
    assert (x2bar[3] == 11.0)
    assert (x2bar[4] == 12.0)
    assert (x2bar[5] == 13.0)
    assert (x2bar[6] == 14.0)
    assert (x2bar[7] == 15.0)

    x2.unserialize([16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 8)
    assert (x2bar[0] == 16.0)
    assert (x2bar[1] == 17.0)
    assert (x2bar[2] == 18.0)
    assert (x2bar[3] == 19.0)
    assert (x2bar[4] == 20.0)
    assert (x2bar[5] == 21.0)
    assert (x2bar[6] == 22.0)
    assert (x2bar[7] == 23.0)

def test_SU_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import SU

    x0 = SU(0)
    x0.unserialize([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = SU(1)
    x1.unserialize([1.0, 2.0, 3.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == complex(1.0, 2.0))

    x1.unserialize([4.0, 5.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == complex(4.0, 5.0))

    x1.unserialize([6.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == complex(6.0, 5.0))

    x2 = SU(2)
    x2.unserialize([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == complex(1.0, 2.0))
    assert (x2hat[0, 1] == complex(3.0, 4.0))
    assert (x2hat[1, 0] == complex(5.0, 6.0))
    assert (x2hat[1, 1] == complex(7.0, 0.0))

    x2.unserialize([8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == complex(8.0, 9.0))
    assert (x2hat[0, 1] == complex(10.0, 11.0))
    assert (x2hat[1, 0] == complex(12.0, 13.0))
    assert (x2hat[1, 1] == complex(14.0, 15.0))

    x2.unserialize([16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == complex(16.0, 17.0))
    assert (x2hat[0, 1] == complex(18.0, 19.0))
    assert (x2hat[1, 0] == complex(20.0, 21.0))
    assert (x2hat[1, 1] == complex(22.0, 23.0))

def test_SU_operator_parenthesis():
    from lielab.domain import SU

    x0 = SU.from_shape(0)
    x0.unserialize([])

    # Out of bounds
    # assert (np.isnan(x0(-1)))
    # assert (np.isnan(x0(0)))
    # assert (np.isnan(x0(1)))

    # Out of bounds
    assert (np.isnan(np.real(x0(0, -1))))
    assert (np.isnan(np.imag(x0(0, -1))))
    assert (np.isnan(np.real(x0(-1, 0))))
    assert (np.isnan(np.imag(x0(-1, 0))))
    assert (np.isnan(np.real(x0(-1, -1))))
    assert (np.isnan(np.imag(x0(-1, -1))))
    assert (np.isnan(np.real(x0(0, 0))))
    assert (np.isnan(np.imag(x0(0, 0))))
    assert (np.isnan(np.real(x0(0, 1))))
    assert (np.isnan(np.imag(x0(0, 1))))
    assert (np.isnan(np.real(x0(1, 0))))
    assert (np.isnan(np.imag(x0(1, 0))))
    assert (np.isnan(np.real(x0(1, 1))))
    assert (np.isnan(np.imag(x0(1, 1))))

    x1 = SU(1)
    x1.unserialize([1.0, 2.0])

    # In bounds
    # assert (x1(0) == 1.0)
    # assert (x1(1) == 2.0)

    # Out of bounds
    # assert (np.isnan(x1(-1)))
    # assert (np.isnan(x1(2)))

    # In bounds
    assert (x1(0, 0) == complex(1.0, 2.0))
    assert (x1(-1, -1) == complex(1.0, 2.0))

    # Out of bounds
    assert (np.isnan(np.real(x1(0, -2))))
    assert (np.isnan(np.imag(x1(0, -2))))
    assert (np.isnan(np.real(x1(-2, 0))))
    assert (np.isnan(np.imag(x1(-2, 0))))
    assert (np.isnan(np.real(x1(-2, -2))))
    assert (np.isnan(np.imag(x1(-2, -2))))
    assert (np.isnan(np.real(x1(0, 1))))
    assert (np.isnan(np.imag(x1(0, 1))))
    assert (np.isnan(np.real(x1(1, 0))))
    assert (np.isnan(np.imag(x1(1, 0))))
    assert (np.isnan(np.real(x1(1, 1))))
    assert (np.isnan(np.imag(x1(1, 1))))

    x2 = SU(2)
    x2.unserialize([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])

    # In bounds
    # assert (x2(0) == 1.0)
    # assert (x2(1) == 2.0)
    # assert (x2(2) == 3.0)
    # assert (x2(3) == 4.0)
    # assert (x2(4) == 5.0)
    # assert (x2(5) == 6.0)
    # assert (x2(6) == 7.0)
    # assert (x2(7) == 8.0)

    # Out of bounds
    # assert (np.isnan(x2(-1)))
    # assert (np.isnan(x2(8)))

    # In bounds
    assert (x2(0, 0) == complex(1.0, 2.0))
    assert (x2(0, 1) == complex(3.0, 4.0))
    assert (x2(1, 0) == complex(5.0, 6.0))
    assert (x2(1, 1) == complex(7.0, 8.0))
    assert (x2(-1, -1) == complex(7.0, 8.0))
    assert (x2(-1, -2) == complex(5.0, 6.0))
    assert (x2(-2, -1) == complex(3.0, 4.0))
    assert (x2(-2, -2) == complex(1.0, 2.0))

    # Out of bounds
    assert (np.isnan(np.real(x2(0, -3))))
    assert (np.isnan(np.imag(x2(0, -3))))
    assert (np.isnan(np.real(x2(-3, 0))))
    assert (np.isnan(np.imag(x2(-3, 0))))
    assert (np.isnan(np.real(x2(-3, -3))))
    assert (np.isnan(np.imag(x2(-3, -3))))
    assert (np.isnan(np.real(x2(0, 2))))
    assert (np.isnan(np.imag(x2(0, 2))))
    assert (np.isnan(np.real(x2(2, 0))))
    assert (np.isnan(np.imag(x2(2, 0))))
    assert (np.isnan(np.real(x2(2, 2))))
    assert (np.isnan(np.imag(x2(2, 2))))

def test_SU_math_ops_SU():
    from lielab.domain import SU

    x1 = SU(2)
    x2 = SU(2)
    x1.unserialize([1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0])
    x2.unserialize([1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0])

    x1_prod_x2 = x1*x2
    x1_prod_x2bar = x1_prod_x2.serialize()
    assert (x1_prod_x2bar.size == 8)
    assert (x1_prod_x2bar[0] == -4.75)
    assert (x1_prod_x2bar[1] == 12.0)
    assert (x1_prod_x2bar[2] == 0.75)
    assert (x1_prod_x2bar[3] == 15.5)
    assert (x1_prod_x2bar[4] == -1.25)
    assert (x1_prod_x2bar[5] == 5.75)
    assert (x1_prod_x2bar[6] == 2.75)
    assert (x1_prod_x2bar[7] == 6.75)

    x1 *= x2
    x1bar = x1.serialize()
    assert (x1bar.size == 8)
    assert (x1bar[0] == -4.75)
    assert (x1bar[1] == 12.0)
    assert (x1bar[2] == 0.75)
    assert (x1bar[3] == 15.5)
    assert (x1bar[4] == -1.25)
    assert (x1bar[5] == 5.75)
    assert (x1bar[6] == 2.75)
    assert (x1bar[7] == 6.75)

    x1.unserialize([1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0])

    x1_inv = x1.inverse()
    x1_invbar = x1_inv.serialize()
    assert (x1_invbar.size == 8)
    assert (np.abs(x1_invbar[0] - -0.25) <= 1e-14)
    assert (np.abs(x1_invbar[1] - 0.25) <= 1e-14)
    assert (np.abs(x1_invbar[2] - 1.0) <= 1e-14)
    assert (np.abs(x1_invbar[3] - -0.75) <= 1e-14)
    assert (np.abs(x1_invbar[4] - 0.25) <= 1e-14)
    assert (np.abs(x1_invbar[5] - -0.25) <= 1e-14)
    assert (np.abs(x1_invbar[6] - -0.5) <= 1e-14)
    assert (np.abs(x1_invbar[7] - 0.25) <= 1e-14)

# def test_SU_project():
#     from lielab.domain import SU

#     rand_2_2 = Eigen::MatrixXcd::Random(2, 2)
#     proj_2_2 = SU.project(rand_2_2)

#     assert (proj_2_2.shape[0] == 2)
#     assert (proj_2_2.shape[1] == 2)
#     assert (proj_2_2(0, 0) == rand_2_2(0, 0))
#     assert (proj_2_2(0, 1) == rand_2_2(0, 1))
#     assert (proj_2_2(1, 0) == rand_2_2(1, 0))
#     assert (proj_2_2(1, 1) == rand_2_2(1, 1))

#     rand_3_3 = Eigen::MatrixXcd::Random(3, 3)
#     proj_3_3 = SU.project(rand_3_3)

#     assert (proj_3_3.shape[0] == 3)
#     assert (proj_3_3.shape[1] == 3)
#     assert (proj_3_3(0, 0) == rand_3_3(0, 0))
#     assert (proj_3_3(0, 1) == rand_3_3(0, 1))
#     assert (proj_3_3(0, 2) == rand_3_3(0, 2))
#     assert (proj_3_3(1, 0) == rand_3_3(1, 0))
#     assert (proj_3_3(1, 1) == rand_3_3(1, 1))
#     assert (proj_3_3(1, 2) == rand_3_3(1, 2))
#     assert (proj_3_3(2, 0) == rand_3_3(2, 0))
#     assert (proj_3_3(2, 1) == rand_3_3(2, 1))
#     assert (proj_3_3(2, 2) == rand_3_3(2, 2))

#     rand_2_3 = Eigen::MatrixXcd::Random(2, 3)
#     proj_2_3 = SU.project(rand_2_3)

#     assert (proj_2_3.shape[0] == 2)
#     assert (proj_2_3.shape[1] == 2)
#     assert (proj_2_3(0, 0) == rand_2_3(0, 0))
#     assert (proj_2_3(0, 1) == rand_2_3(0, 1))
#     assert (proj_2_3(1, 0) == rand_2_3(1, 0))
#     assert (proj_2_3(1, 1) == rand_2_3(1, 1))

#     rand_3_2 = Eigen::MatrixXcd::Random(3, 2)
#     proj_3_2 = SU.project(rand_3_2)

#     assert (proj_3_2.shape[0] == 2)
#     assert (proj_3_2.shape[1] == 2)
#     assert (proj_3_2(0, 0) == rand_3_2(0, 0))
#     assert (proj_3_2(0, 1) == rand_3_2(0, 1))
#     assert (proj_3_2(1, 0) == rand_3_2(1, 0))
#     assert (proj_3_2(1, 1) == rand_3_2(1, 1))

