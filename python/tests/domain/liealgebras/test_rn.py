from lielab.testing import *
import pytest

def test_rn_to_string():
    from lielab.domain import rn

    xzero = rn.from_shape(0)
    assert (xzero.to_string() == "r^nan")
    x0 = rn(0)
    assert (x0.to_string() == "r^0")
    x1 = rn(1)
    assert (x1.to_string() == "r^1")
    x10 = rn(10)
    assert (x10.to_string() == "r^10")

def test_rn_main_initializer():
    from lielab.domain import rn

    xblank = rn()
    assert (xblank.get_dimension() == 0)

    x0 = rn(0)
    assert (x0.get_dimension() == 0)
    x1 = rn(1)
    assert (x1.get_dimension() == 1)
    x10 = rn(10)
    assert (x10.get_dimension() == 10)

def test_rn_matrix_initializer():
    from lielab.domain import rn

    x0 = rn(np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = rn(np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = rn(np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        rn(np.random.rand(2, 3))

    with pytest.raises(RuntimeError):
        rn(np.random.rand(3, 2))

def test_rn_basis_initializer():
    from lielab.domain import rn

    xm10 = rn.basis(-1, 0)
    assert (xm10.get_dimension() == 0)
    xm10bar = xm10.get_vector()
    assert (xm10bar.size == 0)

    x00 = rn.basis(0, 0)
    assert (x00.get_dimension() == 0)
    x00bar = x00.get_vector()
    assert (x00bar.size == 0)

    x10 = rn.basis(1, 0)
    assert (x10.get_dimension() == 0)
    x10bar = x10.get_vector()
    assert (x10bar.size == 0)

    x01 = rn.basis(0, 1)
    assert (x01.get_dimension() == 1)
    x01bar = x01.get_vector()
    assert (x01bar.size == 1)
    assert (x01bar[0] == 1.0)

    x11 = rn.basis(1, 1)
    assert (x11.get_dimension() == 1)
    x11bar = x11.get_vector()
    assert (x11bar.size == 1)
    assert (x11bar[0] == 0.0)

    x21 = rn.basis(2, 1)
    assert (x21.get_dimension() == 1)
    x21bar = x21.get_vector()
    assert (x21bar.size == 1)
    assert (x21bar[0] == 0.0)

    x02 = rn.basis(0, 2)
    assert (x02.get_dimension() == 2)
    x02bar = x02.get_vector()
    assert (x02bar.size == 2)
    assert (x02bar[0] == 1.0)
    assert (x02bar[1] == 0.0)

def test_rn_from_shape_initializer():
    from lielab.domain import rn

    x0 = rn.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = rn.from_shape(1)
    assert (x1.get_dimension() == 0)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = rn.from_shape(2)
    assert (x2.get_dimension() == 1)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_rn_get_dimension():
    from lielab.domain import rn

    veryzero = rn.from_shape(0)
    zero = rn(0)
    one = rn(1)
    two = rn(2)
    three = rn(3)
    four = rn(4)
    five = rn(5)
    six = rn(6)
    seven = rn(7)
    eight = rn(8)

    # Dimensions
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

def test_rn_set_get_vector():
    """
    * Tests the set/get_vector operation.
    """

    from lielab.domain import rn

    xzero = rn.from_shape(0)
    xzero.set_vector([])
    xzerobar = xzero.get_vector()

    assert (xzerobar.size == 0)

    x0 = rn(0)
    x0.set_vector([])
    x0bar = x0.get_vector()

    assert (x0bar.size == 0)

    x2 = rn(2)
    x2.set_vector([1.0, 2.0, 3.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 2)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)

    x2.set_vector([4.0, 5.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 2)
    assert (x2bar[0] == 4.0)
    assert (x2bar[1] == 5.0)

    x2.set_vector([6.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 2)
    assert (x2bar[0] == 6.0)
    assert (x2bar[1] == 5.0)

    x3 = rn(3)
    x3.set_vector([1.0, 2.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 3)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 0.0)

    x3.set_vector([3.0, 4.0, 5.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 3)
    assert (x3bar[0] == 3.0)
    assert (x3bar[1] == 4.0)
    assert (x3bar[2] == 5.0)

    x3.set_vector([6.0, 7.0, 8.0, 9.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 3)
    assert (x3bar[0] == 6.0)
    assert (x3bar[1] == 7.0)
    assert (x3bar[2] == 8.0)

def test_rn_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import rn

    xzero = rn.from_shape(0)
    xzero.set_vector([])
    xzerohat = xzero.get_matrix()

    assert (xzerohat.shape[0] == 0)
    assert (xzerohat.shape[1] == 0)

    x0 = rn(0)
    x0.set_vector([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 1)
    assert (x0hat.shape[1] == 1)
    assert (x0hat[0, 0] == 0.0)

    x1 = rn(1)
    x1.set_vector([1.0, 2.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == 1.0)

    assert (x1hat[0, 0] == 0.0)
    assert (x1hat[1, 0] == 0.0)
    assert (x1hat[1, 1] == 0.0)

    x1.set_vector([3.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == 3.0)

    assert (x1hat[0, 0] == 0.0)
    assert (x1hat[1, 0] == 0.0)
    assert (x1hat[1, 1] == 0.0)

    x1.set_vector([])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == 3.0)

    assert (x1hat[0, 0] == 0.0)
    assert (x1hat[1, 0] == 0.0)
    assert (x1hat[1, 1] == 0.0)

    x2 = rn(2)
    x2.set_vector([1.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == 1.0)
    assert (x2hat[1, 2] == 0.0)

    assert (x2hat[0, 0] == 0.0)
    assert (x2hat[0, 1] == 0.0)
    assert (x2hat[1, 0] == 0.0)
    assert (x2hat[1, 1] == 0.0)
    assert (x2hat[2, 0] == 0.0)
    assert (x2hat[2, 1] == 0.0)
    assert (x2hat[2, 2] == 0.0)

    x2.set_vector([2.0, 3.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == 2.0)
    assert (x2hat[1, 2] == 3.0)

    assert (x2hat[0, 0] == 0.0)
    assert (x2hat[0, 1] == 0.0)
    assert (x2hat[1, 0] == 0.0)
    assert (x2hat[1, 1] == 0.0)
    assert (x2hat[2, 0] == 0.0)
    assert (x2hat[2, 1] == 0.0)
    assert (x2hat[2, 2] == 0.0)

    x2.set_vector([4.0, 5.0, 6.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == 4.0)
    assert (x2hat[1, 2] == 5.0)

    assert (x2hat[0, 0] == 0.0)
    assert (x2hat[0, 1] == 0.0)
    assert (x2hat[1, 0] == 0.0)
    assert (x2hat[1, 1] == 0.0)
    assert (x2hat[2, 0] == 0.0)
    assert (x2hat[2, 1] == 0.0)
    assert (x2hat[2, 2] == 0.0)

def test_rn_operator_parenthesis():
    from lielab.domain import rn

    xzero = rn.from_shape(0)
    xzero.set_vector([])

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

    x1 = rn(1)
    x1.set_vector([1.0])

    # In bounds
    assert (x1(0) == 1.0)
    assert (x1(-1) == 1.0)

    # Out of bounds
    assert (np.isnan(x1(-2)))
    assert (np.isnan(x1(1)))

    # In bounds
    assert (x1(0, 0) == 0.0)
    assert (x1(0, 1) == 1.0)
    assert (x1(1, 0) == 0.0)
    assert (x1(1, 1) == 0.0)
    assert (x1(-1, -1) == 0.0)
    assert (x1(-1, -2) == 0.0)
    assert (x1(-2, -1) == 1.0)
    assert (x1(-2, -2) == 0.0)

    # Out of bounds
    assert (np.isnan(x1(0, -3)))
    assert (np.isnan(x1(-3, 0)))
    assert (np.isnan(x1(-3, -3)))
    assert (np.isnan(x1(0, 2)))
    assert (np.isnan(x1(2, 0)))
    assert (np.isnan(x1(2, 2)))

    x2 = rn(2)
    x2.set_vector([1.0, 2.0])

    # In bounds
    assert (x2(0) == 1.0)
    assert (x2(1) == 2.0)
    assert (x2(-1) == 2.0)
    assert (x2(-2) == 1.0)

    # Out of bounds
    assert (np.isnan(x2(-3)))
    assert (np.isnan(x2(2)))

    # In bounds
    assert (x2(0, 0) == 0.0)
    assert (x2(0, 1) == 0.0)
    assert (x2(0, 2) == 1.0)
    assert (x2(1, 0) == 0.0)
    assert (x2(1, 1) == 0.0)
    assert (x2(1, 2) == 2.0)
    assert (x2(2, 0) == 0.0)
    assert (x2(2, 1) == 0.0)
    assert (x2(2, 2) == 0.0)
    assert (x2(-1, -1) == 0.0)
    assert (x2(-1, -2) == 0.0)
    assert (x2(-1, -3) == 0.0)
    assert (x2(-2, -1) == 2.0)
    assert (x2(-2, -2) == 0.0)
    assert (x2(-2, -3) == 0.0)
    assert (x2(-3, -1) == 1.0)
    assert (x2(-3, -2) == 0.0)
    assert (x2(-3, -3) == 0.0)

    # Out of bounds
    assert (np.isnan(x2(0, -4)))
    assert (np.isnan(x2(-4, 0)))
    assert (np.isnan(x2(-4, -4)))
    assert (np.isnan(x2(0, 3)))
    assert (np.isnan(x2(3, 0)))
    assert (np.isnan(x2(3, 3)))

# TODO: math ops int

def test_rn_math_ops_double():
    from lielab.domain import rn

    x1 = rn(2)
    x1.set_vector([1.25, 2.5])

    x1_lm_2 = 2.0*x1
    assert (x1_lm_2(0) == 2.5)
    assert (x1_lm_2(1) == 5.0)

    x1_rm_2 = x1*2.0
    assert (x1_rm_2(0) == 2.5)
    assert (x1_rm_2(1) == 5.0)

    x1 *= 2.0
    assert (x1(0) == 2.5)
    assert (x1(1) == 5.0)

    x1.set_vector([1.25, 2.5])

    x1_d_2 = x1/2.0
    assert (x1_d_2(0) == 0.625)
    assert (x1_d_2(1) == 1.25)

    x1 /= 2.0
    assert (x1(0) == 0.625)
    assert (x1(1) == 1.25)

def test_rn_math_ops_rn():
    from lielab.domain import rn

    x1 = rn(2)
    x2 = rn(2)
    x1.set_vector([1.0, 2.0])
    x2.set_vector([1.25, 2.5])

    x1_add_x2 = x1 + x2
    assert (x1_add_x2(0) == 2.25)
    assert (x1_add_x2(1) == 4.5)

    x1 += x2
    assert (x1(0) == 2.25)
    assert (x1(1) == 4.5)

    x1.set_vector([1.0, 2.0])

    x1_sub_x2 = x1 - x2
    assert (x1_sub_x2(0) == -0.25)
    assert (x1_sub_x2(1) == -0.5)

    x1 -= x2
    assert (x1(0) == -0.25)
    assert (x1(1) == -0.5)

    x1.set_vector([1.0, 2.0])

    x1_unary_sub = (-x1)
    assert (x1_unary_sub(0) == -1.0)
    assert (x1_unary_sub(1) == -2.0)

def test_rn_from_vector():
    """
    * Tests the from_vector operation.
    """

    from lielab.domain import rn

    x0 = rn.from_vector([])
    x0bar = x0.get_vector()

    assert (x0.get_shape() == 1)
    assert (x0bar.size == 0)

    x1 = rn.from_vector([1.0])
    x1bar = x1.get_vector()

    assert (x1.get_shape() == 2)
    assert (x1bar.size == 1)
    assert (x1bar[0] == 1.0)

    x2 = rn.from_vector([1.0, 2.0])
    x2bar = x2.get_vector()

    assert (x2.get_shape() == 3)
    assert (x2bar.size == 2)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)

    x3 = rn.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.get_vector()

    assert (x3.get_shape() == 4)
    assert (x3bar.size == 3)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)

def test_rn_project():
    from lielab.domain import rn

    rand_2_2 = np.random.rand(2, 2)
    proj_2_2 = rn.project(rand_2_2)

    assert (proj_2_2.shape[0] == 2)
    assert (proj_2_2.shape[1] == 2)
    assert (proj_2_2[0, 0] == 0.0)
    assert (proj_2_2[0, 1] == rand_2_2[0, 1])
    assert (proj_2_2[1, 0] == 0.0)
    assert (proj_2_2[1, 1] == 0.0)

    rand_3_3 = np.random.rand(3, 3)
    proj_3_3 = rn.project(rand_3_3)

    assert (proj_3_3.shape[0] == 3)
    assert (proj_3_3.shape[1] == 3)
    assert (proj_3_3[0, 0] == 0.0)
    assert (proj_3_3[0, 1] == 0.0)
    assert (proj_3_3[0, 2] == rand_3_3[0, 2])
    assert (proj_3_3[1, 0] == 0.0)
    assert (proj_3_3[1, 1] == 0.0)
    assert (proj_3_3[1, 2] == rand_3_3[1, 2])
    assert (proj_3_3[2, 0] == 0.0)
    assert (proj_3_3[2, 1] == 0.0)
    assert (proj_3_3[2, 2] == 0.0)

    rand_2_3 = np.random.rand(2, 3)
    proj_2_3 = rn.project(rand_2_3)

    assert (proj_2_3.shape[0] == 2)
    assert (proj_2_3.shape[1] == 2)
    assert (proj_2_3[0, 0] == 0.0)
    assert (proj_2_3[0, 1] == rand_2_3[0, 1])
    assert (proj_2_3[1, 0] == 0.0)
    assert (proj_2_3[1, 1] == 0.0)

    rand_3_2 = np.random.rand(3, 2)
    proj_3_2 = rn.project(rand_3_2)

    assert (proj_3_2.shape[0] == 2)
    assert (proj_3_2.shape[1] == 2)
    assert (proj_3_2[0, 0] == 0.0)
    assert (proj_3_2[0, 1] == rand_3_2[0, 1])
    assert (proj_3_2[1, 0] == 0.0)
    assert (proj_3_2[1, 1] == 0.0)
