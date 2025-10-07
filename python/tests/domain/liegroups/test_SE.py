from lielab.testing import *
import pytest

def test_SE_to_string():
    from lielab.domain import SE

    x0 = SE(0)
    assert (x0.to_string() == "SE(0)")
    x1 = SE(1)
    assert (x1.to_string() == "SE(1)")
    x10 = SE(10)
    assert (x10.to_string() == "SE(10)")

def test_SE_main_initializer():
    from lielab.domain import SE

    xblank = SE()
    assert (xblank.get_dimension() == 0)

    x0 = SE(0)
    assert (x0.get_dimension() == 0)
    x1 = SE(1)
    assert (x1.get_dimension() == 1)
    x10 = SE(10)
    assert (x10.get_dimension() == 55)

def test_SE_matrix_initializer():
    from lielab.domain import SE

    x0 = SE(np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = SE(np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = SE(np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        SE(np.random.rand(2, 3))
    
    with pytest.raises(RuntimeError):
        SE(np.random.rand(3, 2))

def test_SE_from_shape_initializer():
    from lielab.domain import SE

    x0 = SE.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = SE.from_shape(1)
    assert (x1.get_dimension() == 0)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = SE.from_shape(2)
    assert (x2.get_dimension() == 1)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_SE_get_dimension():
    from lielab.domain import SE

    veryzero = SE.from_shape(0)
    zero = SE(0)
    one = SE(1)
    two = SE(2)
    three = SE(3)
    four = SE(4)
    five = SE(5)
    six = SE(6)
    seven = SE(7)
    eight = SE(8)

    # Dimensions
    assert (veryzero.get_dimension() == 0)
    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 1)
    assert (two.get_dimension() == 3)
    assert (three.get_dimension() == 6)
    assert (four.get_dimension() == 10)
    assert (five.get_dimension() == 15)
    assert (six.get_dimension() == 21)
    assert (seven.get_dimension() == 28)
    assert (eight.get_dimension() == 36)

def test_SE_serialize_unserialize():
    """
    * Tests the serialize/unserialize operation.
    """

    from lielab.domain import SE

    xzero = SE.from_shape(0)
    xzero.unserialize([])
    xzerobar = xzero.serialize()

    assert (xzerobar.size == 0)

    x0 = SE(0)
    x0.unserialize([])
    x0bar = x0.serialize()

    assert (x0bar.size == 1)
    assert (x0bar[0] == 1.0)

    x2 = SE(2)
    x2.unserialize([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 9)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 3.0)
    assert (x2bar[3] == 4.0)
    assert (x2bar[4] == 5.0)
    assert (x2bar[5] == 6.0)
    assert (x2bar[6] == 7.0)
    assert (x2bar[7] == 8.0)
    assert (x2bar[8] == 9.0)

    x2.unserialize([10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 9)
    assert (x2bar[0] == 10.0)
    assert (x2bar[1] == 11.0)
    assert (x2bar[2] == 12.0)
    assert (x2bar[3] == 13.0)
    assert (x2bar[4] == 14.0)
    assert (x2bar[5] == 15.0)
    assert (x2bar[6] == 16.0)
    assert (x2bar[7] == 17.0)
    assert (x2bar[8] == 18.0)

    x2.unserialize([19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 9)
    assert (x2bar[0] == 19.0)
    assert (x2bar[1] == 20.0)
    assert (x2bar[2] == 21.0)
    assert (x2bar[3] == 22.0)
    assert (x2bar[4] == 23.0)
    assert (x2bar[5] == 24.0)
    assert (x2bar[6] == 25.0)
    assert (x2bar[7] == 26.0)
    assert (x2bar[8] == 18.0)

    x3 = SE(3)
    x3.unserialize([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0])
    x3bar = x3.serialize()

    assert (x3bar.size == 16)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)
    assert (x3bar[3] == 4.0)
    assert (x3bar[4] == 5.0)
    assert (x3bar[5] == 6.0)
    assert (x3bar[6] == 7.0)
    assert (x3bar[7] == 8.0)
    assert (x3bar[8] == 9.0)
    assert (x3bar[9] == 10.0)
    assert (x3bar[10] == 11.0)
    assert (x3bar[11] == 12.0)
    assert (x3bar[12] == 13.0)
    assert (x3bar[13] == 14.0)
    assert (x3bar[14] == 15.0)
    assert (x3bar[15] == 1.0)

    x3.unserialize([16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0])
    x3bar = x3.serialize()

    assert (x3bar.size == 16)
    assert (x3bar[0] == 16.0)
    assert (x3bar[1] == 17.0)
    assert (x3bar[2] == 18.0)
    assert (x3bar[3] == 19.0)
    assert (x3bar[4] == 20.0)
    assert (x3bar[5] == 21.0)
    assert (x3bar[6] == 22.0)
    assert (x3bar[7] == 23.0)
    assert (x3bar[8] == 24.0)
    assert (x3bar[9] == 25.0)
    assert (x3bar[10] == 26.0)
    assert (x3bar[11] == 27.0)
    assert (x3bar[12] == 28.0)
    assert (x3bar[13] == 29.0)
    assert (x3bar[14] == 30.0)
    assert (x3bar[15] == 31.0)

    x3.unserialize([32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0])
    x3bar = x3.serialize()

    assert (x3bar.size == 16)
    assert (x3bar[0] == 32.0)
    assert (x3bar[1] == 33.0)
    assert (x3bar[2] == 34.0)
    assert (x3bar[3] == 35.0)
    assert (x3bar[4] == 36.0)
    assert (x3bar[5] == 37.0)
    assert (x3bar[6] == 38.0)
    assert (x3bar[7] == 39.0)
    assert (x3bar[8] == 40.0)
    assert (x3bar[9] == 41.0)
    assert (x3bar[10] == 42.0)
    assert (x3bar[11] == 43.0)
    assert (x3bar[12] == 44.0)
    assert (x3bar[13] == 45.0)
    assert (x3bar[14] == 46.0)
    assert (x3bar[15] == 47.0)

def test_SE_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import SE

    xzero = SE.from_shape(0)
    xzero.unserialize([])
    xzerohat = xzero.get_matrix()

    assert (xzerohat.shape[0] == 0)
    assert (xzerohat.shape[1] == 0)

    x0 = SE(0)
    x0.unserialize([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 1)
    assert (x0hat.shape[1] == 1)
    assert (x0hat[0, 0] == 1.0)

    x1 = SE(1)
    x1.unserialize([1.0, 2.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 0] == 1.0)
    assert (x1hat[0, 1] == 2.0)
    assert (x1hat[1, 0] == 0.0)
    assert (x1hat[1, 1] == 1.0)

    x1.unserialize([3.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 0] == 3.0)
    assert (x1hat[0, 1] == 2.0)
    assert (x1hat[1, 0] == 0.0)
    assert (x1hat[1, 1] == 1.0)

    x1.unserialize([])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 0] == 3.0)
    assert (x1hat[0, 1] == 2.0)
    assert (x1hat[1, 0] == 0.0)
    assert (x1hat[1, 1] == 1.0)

    x2 = SE(2)
    x2.unserialize([1.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 0] == 1.0)
    assert (x2hat[0, 1] == 0.0)
    assert (x2hat[0, 2] == 0.0)
    assert (x2hat[1, 0] == 0.0)
    assert (x2hat[1, 1] == 1.0)
    assert (x2hat[1, 2] == 0.0)
    assert (x2hat[2, 0] == 0.0)
    assert (x2hat[2, 1] == 0.0)
    assert (x2hat[2, 2] == 1.0)

    x2.unserialize([2.0, 3.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 0] == 2.0)
    assert (x2hat[0, 1] == 3.0)
    assert (x2hat[0, 2] == 0.0)
    assert (x2hat[1, 0] == 0.0)
    assert (x2hat[1, 1] == 1.0)
    assert (x2hat[1, 2] == 0.0)
    assert (x2hat[2, 0] == 0.0)
    assert (x2hat[2, 1] == 0.0)
    assert (x2hat[2, 2] == 1.0)

    x2.unserialize([4.0, 5.0, 6.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 0] == 4.0)
    assert (x2hat[0, 1] == 5.0)
    assert (x2hat[0, 2] == 6.0)
    assert (x2hat[1, 0] == 0.0)
    assert (x2hat[1, 1] == 1.0)
    assert (x2hat[1, 2] == 0.0)
    assert (x2hat[2, 0] == 0.0)
    assert (x2hat[2, 1] == 0.0)
    assert (x2hat[2, 2] == 1.0)

def test_SE_operator_parenthesis():
    from lielab.domain import SE

    xzero = SE.from_shape(0)
    xzero.unserialize([])

    # Out of bounds
    # assert (np.isnan(xzero(-1)))
    # assert (np.isnan(xzero(0)))
    # assert (np.isnan(xzero(1)))

    # Out of bounds
    assert (np.isnan(xzero(0, -1)))
    assert (np.isnan(xzero(-1, 0)))
    assert (np.isnan(xzero(-1, -1)))
    assert (np.isnan(xzero(0, 0)))
    assert (np.isnan(xzero(0, 1)))
    assert (np.isnan(xzero(1, 0)))
    assert (np.isnan(xzero(1, 1)))

    x1 = SE(1)
    x1.unserialize([1.0])

    # In bounds
    # assert (x1(0) == 1.0)

    # Out of bounds
    # assert (np.isnan(x1(-1)))
    # assert (np.isnan(x1(1)))

    # In bounds
    assert (x1(0, 0) == 1.0)
    assert (x1(0, 1) == 0.0)
    assert (x1(1, 0) == 0.0)
    assert (x1(1, 1) == 1.0)
    assert (x1(-1, -1) == 1.0)
    assert (x1(-1, -2) == 0.0)
    assert (x1(-2, -1) == 0.0)
    assert (x1(-2, -2) == 1.0)

    # Out of bounds
    assert (np.isnan(x1(0, -3)))
    assert (np.isnan(x1(-3, 0)))
    assert (np.isnan(x1(-3, -3)))
    assert (np.isnan(x1(0, 2)))
    assert (np.isnan(x1(2, 0)))
    assert (np.isnan(x1(2, 2)))

    x2 = SE(2)
    x2.unserialize([1.0, 2.0, 3.0])

    # In bounds
    # assert (x2(0) == 1.0)
    # assert (x2(1) == 2.0)
    # assert (x2(2) == 3.0)

    # Out of bounds
    # assert (np.isnan(x2(-1)))
    # assert (np.isnan(x2(3)))

    # In bounds
    assert (x2(0, 0) == 1.0)
    assert (x2(0, 1) == 2.0)
    assert (x2(0, 2) == 3.0)
    assert (x2(1, 0) == 0.0)
    assert (x2(1, 1) == 1.0)
    assert (x2(1, 2) == 0.0)
    assert (x2(2, 0) == 0.0)
    assert (x2(2, 1) == 0.0)
    assert (x2(2, 2) == 1.0)
    assert (x2(-1, -1) == 1.0)
    assert (x2(-1, -2) == 0.0)
    assert (x2(-1, -3) == 0.0)
    assert (x2(-2, -1) == 0.0)
    assert (x2(-2, -2) == 1.0)
    assert (x2(-2, -3) == 0.0)
    assert (x2(-3, -1) == 3.0)
    assert (x2(-3, -2) == 2.0)
    assert (x2(-3, -3) == 1.0)

    # Out of bounds
    assert (np.isnan(x2(0, -4)))
    assert (np.isnan(x2(-4, 0)))
    assert (np.isnan(x2(-4, -4)))
    assert (np.isnan(x2(0, 3)))
    assert (np.isnan(x2(3, 0)))
    assert (np.isnan(x2(3, 3)))

def test_SE_math_ops_SE():
    from lielab.domain import SE

    x1 = SE(2)
    x2 = SE(2)
    x1.unserialize([1.0, 2.0, 3.0])
    x2.unserialize([1.25, 2.5, 3.75])

    x1_prod_x2 = x1*x2
    x1_prod_x2bar = x1_prod_x2.serialize()
    assert (x1_prod_x2bar.size == 9)
    assert (x1_prod_x2bar[0] == 1.25)
    assert (x1_prod_x2bar[1] == 4.5)
    assert (x1_prod_x2bar[2] == 6.75)
    assert (x1_prod_x2bar[3] == 0.0)
    assert (x1_prod_x2bar[4] == 1.0)
    assert (x1_prod_x2bar[5] == 0.0)
    assert (x1_prod_x2bar[6] == 0.0)
    assert (x1_prod_x2bar[7] == 0.0)
    assert (x1_prod_x2bar[8] == 1.0)

    x1 *= x2
    x1bar = x1.serialize()
    assert (x1bar.size == 9)
    assert (x1bar[0] == 1.25)
    assert (x1bar[1] == 4.5)
    assert (x1bar[2] == 6.75)
    assert (x1bar[3] == 0.0)
    assert (x1bar[4] == 1.0)
    assert (x1bar[5] == 0.0)
    assert (x1bar[6] == 0.0)
    assert (x1bar[7] == 0.0)
    assert (x1bar[8] == 1.0)

    x1.unserialize([1.0, 2.0, 3.0])

    x1_inv = x1.inverse()
    x1_invbar = x1_inv.serialize()
    assert (x1_invbar.size == 9)
    assert (np.abs(x1_invbar[0] - 1.0) <= 1e-14)
    assert (np.abs(x1_invbar[1] - -2.0) <= 1e-14)
    assert (np.abs(x1_invbar[2] - -3.0) <= 1e-14)
    assert (np.abs(x1_invbar[3] - 0.0) <= 1e-14)
    assert (np.abs(x1_invbar[4] - 1.0) <= 1e-14)
    assert (np.abs(x1_invbar[5] - 0.0) <= 1e-14)
    assert (np.abs(x1_invbar[6] - 0.0) <= 1e-14)
    assert (np.abs(x1_invbar[7] - 0.0) <= 1e-14)
    assert (np.abs(x1_invbar[8] - 1.0) <= 1e-14)

# def test_SE_project():
#     from lielab.domain import SE

#     rand_2_2 = np.random.rand(2, 2)
#     proj_2_2 = SO::project(rand_2_2)

#     assert (proj_2_2.shape[0] == 2)
#     assert (proj_2_2.shape[1] == 2)
#     assert (std::abs(proj_2_2.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14))

#     rand_3_3 = np.random.rand(3, 3)
#     proj_3_3 = SO::project(rand_3_3)

#     assert (proj_3_3.shape[0] == 3)
#     assert (proj_3_3.shape[1] == 3)
#     assert (std::abs(proj_3_3.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14))

#     rand_2_3 = np.random.rand(2, 3)
#     proj_2_3 = SO::project(rand_2_3)

#     assert (proj_2_3.shape[0] == 2)
#     assert (proj_2_3.shape[1] == 2)
#     assert (std::abs(proj_2_3.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14))

#     rand_3_2 = np.random.rand(3, 2)
#     proj_3_2 = SO::project(rand_3_2)

#     assert (proj_3_2.shape[0] == 2)
#     assert (proj_3_2.shape[1] == 2)
#     assert (std::abs(proj_3_2.determinant()), Catch::Matchers::WithinAbs(1.0, 1e-14))
