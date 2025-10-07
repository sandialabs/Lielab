from lielab.testing import *
import pytest

def test_RN_to_string():
    from lielab.domain import RN

    xzero = RN.from_shape(0)
    assert (xzero.to_string() == "R^nan")
    x0 = RN(0)
    assert (x0.to_string() == "R^0")
    x1 = RN(1)
    assert (x1.to_string() == "R^1")
    x10 = RN(10)
    assert (x10.to_string() == "R^10")

def test_RN_main_initializer():
    from lielab.domain import RN

    xblank = RN()
    assert (xblank.get_dimension() == 0)

    x0 = RN(0)
    assert (x0.get_dimension() == 0)
    x1 = RN(1)
    assert (x1.get_dimension() == 1)
    x10 = RN(10)
    assert (x10.get_dimension() == 10)

def test_RN_matrix_initializer():
    from lielab.domain import RN

    x0 = RN(np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = RN(np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = RN(np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        RN(np.random.rand(2, 3))
    
    with pytest.raises(RuntimeError):
        RN(np.random.rand(3, 2))

def test_RN_from_shape_initializer():
    from lielab.domain import RN

    x0 = RN.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = RN.from_shape(1)
    assert (x1.get_dimension() == 0)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = RN.from_shape(2)
    assert (x2.get_dimension() == 1)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_RN_get_dimension():
    from lielab.domain import RN

    veryzero = RN.from_shape(0)
    zero = RN(0)
    one = RN(1)
    two = RN(2)
    three = RN(3)
    four = RN(4)
    five = RN(5)
    six = RN(6)
    seven = RN(7)
    eight = RN(8)

    assert (veryzero.get_dimension() == 0)
    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 1)
    assert (two.get_dimension() == 2)
    assert (three.get_dimension() == 3)
    assert (four.get_dimension() == 4)
    assert (five.get_dimension() == 5)
    assert (six.get_dimension() == 6)
    assert (seven.get_dimension() == 7)
    assert (eight.get_dimension() == 8)

def test_RN_serialize_unserialize():
    """
    * Tests the serialize/unserialize operation.
    """

    from lielab.domain import RN

    xzero = RN.from_shape(0)
    xzero.unserialize([])
    xzerobar = xzero.serialize()

    assert (xzerobar.size == 0)

    x0 = RN(0)
    x0.unserialize([])
    x0bar = x0.serialize()

    assert (x0bar.size == 0)

    x2 = RN(2)
    x2.unserialize([1.0, 2.0, 3.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 2)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)

    x2.unserialize([4.0, 5.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 2)
    assert (x2bar[0] == 4.0)
    assert (x2bar[1] == 5.0)

    x2.unserialize([6.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 2)
    assert (x2bar[0] == 6.0)
    assert (x2bar[1] == 5.0)

    x3 = RN(3)
    x3.unserialize([1.0, 2.0])
    x3bar = x3.serialize()

    assert (x3bar.size == 3)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 0.0)

    x3.unserialize([3.0, 4.0, 5.0])
    x3bar = x3.serialize()

    assert (x3bar.size == 3)
    assert (x3bar[0] == 3.0)
    assert (x3bar[1] == 4.0)
    assert (x3bar[2] == 5.0)

    x3.unserialize([6.0, 7.0, 8.0, 9.0])
    x3bar = x3.serialize()

    assert (x3bar.size == 3)
    assert (x3bar[0] == 6.0)
    assert (x3bar[1] == 7.0)
    assert (x3bar[2] == 8.0)

def test_RN_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import RN

    xzero = RN.from_shape(0)
    xzero.unserialize([])
    xzerohat = xzero.get_matrix()

    assert (xzerohat.shape[0] == 0)
    assert (xzerohat.shape[1] == 0)

    x0 = RN(0)
    x0.unserialize([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 1)
    assert (x0hat.shape[1] == 1)
    assert (x0hat[0, 0] == 1.0)

    x1 = RN(1)
    x1.unserialize([1.0, 2.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == 1.0)

    assert (x1hat[0, 0] == 1.0)
    assert (x1hat[1, 0] == 0.0)
    assert (x1hat[1, 1] == 1.0)

    x1.unserialize([3.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == 3.0)

    assert (x1hat[0, 0] == 1.0)
    assert (x1hat[1, 0] == 0.0)
    assert (x1hat[1, 1] == 1.0)

    x1.unserialize([])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == 3.0)

    assert (x1hat[0, 0] == 1.0)
    assert (x1hat[1, 0] == 0.0)
    assert (x1hat[1, 1] == 1.0)

    x2 = RN(2)
    x2.unserialize([1.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == 1.0)
    assert (x2hat[1, 2] == 0.0)

    assert (x2hat[0, 0] == 1.0)
    assert (x2hat[0, 1] == 0.0)
    assert (x2hat[1, 0] == 0.0)
    assert (x2hat[1, 1] == 1.0)
    assert (x2hat[2, 0] == 0.0)
    assert (x2hat[2, 1] == 0.0)
    assert (x2hat[2, 2] == 1.0)

    x2.unserialize([2.0, 3.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == 2.0)
    assert (x2hat[1, 2] == 3.0)

    assert (x2hat[0, 0] == 1.0)
    assert (x2hat[0, 1] == 0.0)
    assert (x2hat[1, 0] == 0.0)
    assert (x2hat[1, 1] == 1.0)
    assert (x2hat[2, 0] == 0.0)
    assert (x2hat[2, 1] == 0.0)
    assert (x2hat[2, 2] == 1.0)

    x2.unserialize([4.0, 5.0, 6.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == 4.0)
    assert (x2hat[1, 2] == 5.0)

    assert (x2hat[0, 0] == 1.0)
    assert (x2hat[0, 1] == 0.0)
    assert (x2hat[1, 0] == 0.0)
    assert (x2hat[1, 1] == 1.0)
    assert (x2hat[2, 0] == 0.0)
    assert (x2hat[2, 1] == 0.0)
    assert (x2hat[2, 2] == 1.0)

def test_RN_operator_parenthesis():
    from lielab.domain import RN

    xzero = RN.from_shape(0)
    xzero.unserialize([])

    # Out of bounds
    assert (np.isnan(xzero(-1)))
    assert (np.isnan(xzero(0)))
    assert (np.isnan(xzero(1)))

    # Out of bounds
    assert (np.isnan(xzero(0, -1)))
    assert (np.isnan(xzero(-1, 0)))
    assert (np.isnan(xzero(-1, -1)))
    assert (np.isnan(xzero(0, 0)))
    assert (np.isnan(xzero(0, 1)))
    assert (np.isnan(xzero(1, 0)))
    assert (np.isnan(xzero(1, 1)))

    x1 = RN(1)
    x1.unserialize([1.0])

    # In bounds
    assert (x1(0) == 1.0)
    assert (x1(-1) == 1.0)

    # Out of bounds
    assert (np.isnan(x1(-2)))
    assert (np.isnan(x1(1)))

    # In bounds
    assert (x1(0, 0) == 1.0)
    assert (x1(0, 1) == 1.0)
    assert (x1(1, 0) == 0.0)
    assert (x1(1, 1) == 1.0)
    assert (x1(-1, -1) == 1.0)
    assert (x1(-1, -2) == 0.0)
    assert (x1(-2, -1) == 1.0)
    assert (x1(-2, -2) == 1.0)

    # Out of bounds
    assert (np.isnan(x1(0, -3)))
    assert (np.isnan(x1(-3, 0)))
    assert (np.isnan(x1(-3, -3)))
    assert (np.isnan(x1(0, 2)))
    assert (np.isnan(x1(2, 0)))
    assert (np.isnan(x1(2, 2)))

    x2 = RN(2)
    x2.unserialize([1.0, 2.0])

    # In bounds
    assert (x2(0) == 1.0)
    assert (x2(1) == 2.0)
    assert (x2(-1) == 2.0)
    assert (x2(-2) == 1.0)

    # Out of bounds
    assert (np.isnan(x2(-3)))
    assert (np.isnan(x2(2)))

    # In bounds
    assert (x2(0, 0) == 1.0)
    assert (x2(0, 1) == 0.0)
    assert (x2(0, 2) == 1.0)
    assert (x2(1, 0) == 0.0)
    assert (x2(1, 1) == 1.0)
    assert (x2(1, 2) == 2.0)
    assert (x2(2, 0) == 0.0)
    assert (x2(2, 1) == 0.0)
    assert (x2(2, 2) == 1.0)
    assert (x2(-1, -1) == 1.0)
    assert (x2(-1, -2) == 0.0)
    assert (x2(-1, -3) == 0.0)
    assert (x2(-2, -1) == 2.0)
    assert (x2(-2, -2) == 1.0)
    assert (x2(-2, -3) == 0.0)
    assert (x2(-3, -1) == 1.0)
    assert (x2(-3, -2) == 0.0)
    assert (x2(-3, -3) == 1.0)

    # Out of bounds
    assert (np.isnan(x2(0, -4)))
    assert (np.isnan(x2(-4, 0)))
    assert (np.isnan(x2(-4, -4)))
    assert (np.isnan(x2(0, 3)))
    assert (np.isnan(x2(3, 0)))
    assert (np.isnan(x2(3, 3)))

def test_RN_math_ops_RN():
    from lielab.domain import RN

    x1 = RN(4)
    x2 = RN(4)
    x1.unserialize([1.0, 2.0, 3.0, 4.0])
    x2.unserialize([1.25, 2.5, 3.75, 1.0])

    x1_prod_x2 = x1*x2
    x1_prod_x2bar = x1_prod_x2.serialize()
    assert (x1_prod_x2bar.size == 4)
    assert (x1_prod_x2bar[0] == 2.25)
    assert (x1_prod_x2bar[1] == 4.5)
    assert (x1_prod_x2bar[2] == 6.75)
    assert (x1_prod_x2bar[3] == 5.0)

    x1 *= x2
    x1bar = x1.serialize()
    assert (x1bar.size == 4)
    assert (x1bar[0] == 2.25)
    assert (x1bar[1] == 4.5)
    assert (x1bar[2] == 6.75)
    assert (x1bar[3] == 5.0)

    x1.unserialize([1.0, 2.0, 3.0, 4.0])

    x1_inv = x1.inverse()
    x1_invbar = x1_inv.serialize()
    assert (x1_invbar.size == 4)
    assert (x1_invbar[0] == -1.0)
    assert (x1_invbar[1] == -2.0)
    assert (x1_invbar[2] == -3.0)
    assert (x1_invbar[3] == -4.0)

def test_RN_from_vector():
    """
    * Tests the from_vector operation.
    """

    from lielab.domain import RN

    x0 = RN.from_vector([])
    x0bar = x0.serialize()

    assert (x0.get_shape() == 1)
    assert (x0bar.size == 0)

    x1 = RN.from_vector([1.0])
    x1bar = x1.serialize()

    assert (x1.get_shape() == 2)
    assert (x1bar.size == 1)
    assert (x1bar[0] == 1.0)

    x2 = RN.from_vector([1.0, 2.0])
    x2bar = x2.serialize()

    assert (x2.get_shape() == 3)
    assert (x2bar.size == 2)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)

    x3 = RN.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.serialize()

    assert (x3.get_shape() == 4)
    assert (x3bar.size == 3)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)

def test_RN_project():
    from lielab.domain import RN

    rand_2_2 = np.random.rand(2, 2)
    proj_2_2 = RN.project(rand_2_2)

    assert (proj_2_2.shape[0] == 2)
    assert (proj_2_2.shape[1] == 2)
    assert (proj_2_2[0, 0] == 1.0)
    assert (proj_2_2[0, 1] == rand_2_2[0, 1])
    assert (proj_2_2[1, 0] == 0.0)
    assert (proj_2_2[1, 1] == 1.0)

    rand_3_3 = np.random.rand(3, 3)
    proj_3_3 = RN.project(rand_3_3)

    assert (proj_3_3.shape[0] == 3)
    assert (proj_3_3.shape[1] == 3)
    assert (proj_3_3[0, 0] == 1.0)
    assert (proj_3_3[0, 1] == 0.0)
    assert (proj_3_3[0, 2] == rand_3_3[0, 2])
    assert (proj_3_3[1, 0] == 0.0)
    assert (proj_3_3[1, 1] == 1.0)
    assert (proj_3_3[1, 2] == rand_3_3[1, 2])
    assert (proj_3_3[2, 0] == 0.0)
    assert (proj_3_3[2, 1] == 0.0)
    assert (proj_3_3[2, 2] == 1.0)

    rand_2_3 = np.random.rand(2, 3)
    proj_2_3 = RN.project(rand_2_3)

    assert (proj_2_3.shape[0] == 2)
    assert (proj_2_3.shape[1] == 2)
    assert (proj_2_3[0, 0] == 1.0)
    assert (proj_2_3[0, 1] == rand_2_3[0, 1])
    assert (proj_2_3[1, 0] == 0.0)
    assert (proj_2_3[1, 1] == 1.0)

    rand_3_2 = np.random.rand(3, 2)
    proj_3_2 = RN.project(rand_3_2)

    assert (proj_3_2.shape[0] == 2)
    assert (proj_3_2.shape[1] == 2)
    assert (proj_3_2[0, 0] == 1.0)
    assert (proj_3_2[0, 1] == rand_3_2[0, 1])
    assert (proj_3_2[1, 0] == 0.0)
    assert (proj_3_2[1, 1] == 1.0)
