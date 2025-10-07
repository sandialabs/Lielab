from lielab.testing import *
import pytest

def test_glr_to_string():
    from lielab.domain import glr

    x0 = glr(0)
    assert (x0.to_string() == "gl(0, R)")
    x1 = glr(1)
    assert (x1.to_string() == "gl(1, R)")
    x10 = glr(10)
    assert (x10.to_string() == "gl(10, R)")

def test_glr_main_initializer():
    from lielab.domain import glr

    xblank = glr()
    assert (xblank.get_dimension() == 0)

    x0 = glr(0)
    assert (x0.get_dimension() == 0)
    x1 = glr(1)
    assert (x1.get_dimension() == 1)
    x10 = glr(10)
    assert (x10.get_dimension() == 100)

def test_glr_matrix_initializer():
    from lielab.domain import glr

    x0 = glr(np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = glr(np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = glr(np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        glr(np.random.rand(2, 3))

    with pytest.raises(RuntimeError):
        glr(np.random.rand(3, 2))

def test_glr_basis_initializer():
    from lielab.domain import glr

    xm10 = glr.basis(-1, 0)
    assert (xm10.get_dimension() == 0)
    xm10bar = xm10.get_vector()
    assert (xm10bar.size == 0)

    x00 = glr.basis(0, 0)
    assert (x00.get_dimension() == 0)
    x00bar = x00.get_vector()
    assert (x00bar.size == 0)

    x10 = glr.basis(1, 0)
    assert (x10.get_dimension() == 0)
    x10bar = x10.get_vector()
    assert (x10bar.size == 0)

    x01 = glr.basis(0, 1)
    assert (x01.get_dimension() == 1)
    x01bar = x01.get_vector()
    assert (x01bar.size == 1)
    assert (x01bar[0] == 1.0)

    x11 = glr.basis(1, 1)
    assert (x11.get_dimension() == 1)
    x11bar = x11.get_vector()
    assert (x11bar.size == 1)
    assert (x11bar[0] == 0.0)

    x02 = glr.basis(0, 2)
    assert (x02.get_dimension() == 4)
    x02bar = x02.get_vector()
    assert (x02bar.size == 4)
    assert (x02bar[0] == 1.0)
    assert (x02bar[1] == 0.0)
    assert (x02bar[2] == 0.0)
    assert (x02bar[3] == 0.0)

def test_glr_from_shape_initializer():
    from lielab.domain import glr

    x0 = glr.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = glr.from_shape(1)
    assert (x1.get_dimension() == 1)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = glr.from_shape(2)
    assert (x2.get_dimension() == 4)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_glr_get_dimension():
    from lielab.domain import glr

    zero = glr(0)
    one = glr(1)
    two = glr(2)
    three = glr(3)
    four = glr(4)
    five = glr(5)
    six = glr(6)
    seven = glr(7)
    eight = glr(8)

    # Dimensions
    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 1)
    assert (two.get_dimension() == 4)
    assert (three.get_dimension() == 9)
    assert (four.get_dimension() == 16)
    assert (five.get_dimension() == 25)
    assert (six.get_dimension() == 36)
    assert (seven.get_dimension() == 49)
    assert (eight.get_dimension() == 64)

def test_glr_set_get_vector():
    """
    * Tests the set/get_vector operation.
    """

    from lielab.domain import glr

    x0 = glr(0)
    x0.set_vector([])
    x0bar = x0.get_vector()

    assert (x0bar.size == 0)

    x1 = glr(1)
    x1.set_vector([1.0, 2.0])
    x1bar = x1.get_vector()

    assert (x1bar.size == 1)
    assert (x1bar[0] == 1.0)

    x1.set_vector([3.0])
    x1bar = x1.get_vector()

    assert (x1bar.size == 1)
    assert (x1bar[0] == 3.0)

    x1.set_vector([])
    x1bar = x1.get_vector()

    assert (x1bar.size == 1)
    assert (x1bar[0] == 3.0)

    x2 = glr(2)
    x2.set_vector([1.0, 2.0, 3.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 3.0)
    assert (x2bar[3] == 0.0)

    x2.set_vector([4.0, 5.0, 6.0, 7.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 4.0)
    assert (x2bar[1] == 5.0)
    assert (x2bar[2] == 6.0)
    assert (x2bar[3] == 7.0)

    x2.set_vector([8.0, 9.0, 10.0, 11.0, 12.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 8.0)
    assert (x2bar[1] == 9.0)
    assert (x2bar[2] == 10.0)
    assert (x2bar[3] == 11.0)

def test_glr_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import glr

    x0 = glr(0)
    x0.set_vector([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = glr(1)
    x1.set_vector([])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == 0.0)

    x1.set_vector([1.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == 1.0)

    x1.set_vector([2.0, 3.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == 2.0)

    x2 = glr(2)
    x2.set_vector([1.0, 2.0, 3.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 1.0)
    assert (x2hat[0, 1] == 2.0)
    assert (x2hat[1, 0] == 3.0)
    assert (x2hat[1, 1] == 0.0)

    x2.set_vector([4.0, 5.0, 6.0, 7.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 4.0)
    assert (x2hat[0, 1] == 5.0)
    assert (x2hat[1, 0] == 6.0)
    assert (x2hat[1, 1] == 7.0)

    x2.set_vector([8.0, 9.0, 10.0, 11.0, 12.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 8.0)
    assert (x2hat[0, 1] == 9.0)
    assert (x2hat[1, 0] == 10.0)
    assert (x2hat[1, 1] == 11.0)

def test_glr_operator_parenthesis():
    from lielab.domain import glr

    x0 = glr.from_shape(0)
    x0.set_vector([])

    # Out of bounds
    assert (np.isnan(x0(-1)))
    assert (np.isnan(x0(0)))
    assert (np.isnan(x0(1)))

    # Out of bounds
    assert (np.isnan(x0(0, -1)))
    assert (np.isnan(x0(-1, 0)))
    assert (np.isnan(x0(-1, -1)))
    assert (np.isnan(x0(0, 0)))
    assert (np.isnan(x0(0, 1)))
    assert (np.isnan(x0(1, 0)))
    assert (np.isnan(x0(1, 1)))

    x1 = glr(1)
    x1.set_vector([1.0])

    # In bounds
    assert (x1(0) == 1.0)
    assert (x1(-1) == 1.0)

    # Out of bounds
    assert (np.isnan(x1(-2)))
    assert (np.isnan(x1(1)))

    # In bounds
    assert (x1(0, 0) == 1.0)
    assert (x1(-1, -1) == 1.0)

    # Out of bounds
    assert (np.isnan(x1(0, -2)))
    assert (np.isnan(x1(-2, 0)))
    assert (np.isnan(x1(-2, -2)))
    assert (np.isnan(x1(0, 1)))
    assert (np.isnan(x1(1, 0)))
    assert (np.isnan(x1(1, 1)))

    x2 = glr(2)
    x2.set_vector([1.0, 2.0, 3.0, 4.0])

    # In bounds
    assert (x2(0) == 1.0)
    assert (x2(1) == 2.0)
    assert (x2(2) == 3.0)
    assert (x2(3) == 4.0)
    assert (x2(-1) == 4.0)
    assert (x2(-2) == 3.0)
    assert (x2(-3) == 2.0)
    assert (x2(-4) == 1.0)

    # Out of bounds
    assert (np.isnan(x2(-5)))
    assert (np.isnan(x2(4)))

    # In bounds
    assert (x2(0, 0) == 1.0)
    assert (x2(0, 1) == 2.0)
    assert (x2(1, 0) == 3.0)
    assert (x2(1, 1) == 4.0)
    assert (x2(-1, -1) == 4.0)
    assert (x2(-1, -2) == 3.0)
    assert (x2(-2, -1) == 2.0)
    assert (x2(-2, -2) == 1.0)

    # Out of bounds
    assert (np.isnan(x2(0, -3)))
    assert (np.isnan(x2(-3, 0)))
    assert (np.isnan(x2(-3, -3)))
    assert (np.isnan(x2(0, 2)))
    assert (np.isnan(x2(2, 0)))
    assert (np.isnan(x2(2, 2)))

# TODO: math ops int

def test_glr_math_ops_double():
    from lielab.domain import glr

    x1 = glr(2)
    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_lm_2 = 2.0*x1
    assert (x1_lm_2(0, 0) == 2.5)
    assert (x1_lm_2(0, 1) == 5.0)
    assert (x1_lm_2(1, 0) == 7.5)
    assert (x1_lm_2(1, 1) == 2.0)

    x1_rm_2 = x1*2.0
    assert (x1_rm_2(0, 0) == 2.5)
    assert (x1_rm_2(0, 1) == 5.0)
    assert (x1_rm_2(1, 0) == 7.5)
    assert (x1_rm_2(1, 1) == 2.0)

    x1 *= 2.0
    assert (x1(0, 0) == 2.5)
    assert (x1(0, 1) == 5.0)
    assert (x1(1, 0) == 7.5)
    assert (x1(1, 1) == 2.0)

    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_d_2 = x1/2.0
    assert (x1_d_2(0, 0) == 0.625)
    assert (x1_d_2(0, 1) == 1.25)
    assert (x1_d_2(1, 0) == 1.875)
    assert (x1_d_2(1, 1) == 0.5)

    x1 /= 2.0
    assert (x1(0, 0) == 0.625)
    assert (x1(0, 1) == 1.25)
    assert (x1(1, 0) == 1.875)
    assert (x1(1, 1) == 0.5)

def test_glr_math_ops_cn():
    from lielab.domain import glr

    x1 = glr(2)
    x2 = glr(2)
    x1.set_vector([1.0, 2.0, 3.0, 4.0])
    x2.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_add_x2 = x1 + x2
    assert (x1_add_x2(0, 0) == 2.25)
    assert (x1_add_x2(0, 1) == 4.5)
    assert (x1_add_x2(1, 0) == 6.75)
    assert (x1_add_x2(1, 1) == 5.0)

    x1 += x2
    assert (x1(0, 0) == 2.25)
    assert (x1(0, 1) == 4.5)
    assert (x1(1, 0) == 6.75)
    assert (x1(1, 1) == 5.0)

    x1.set_vector([1.0, 2.0, 3.0, 4.0])

    x1_sub_x2 = x1 - x2
    assert (x1_sub_x2(0, 0) == -0.25)
    assert (x1_sub_x2(0, 1) == -0.5)
    assert (x1_sub_x2(1, 0) == -0.75)
    assert (x1_sub_x2(1, 1) == 3.0)

    x1 -= x2
    assert (x1(0, 0) == -0.25)
    assert (x1(0, 1) == -0.5)
    assert (x1(1, 0) == -0.75)
    assert (x1(1, 1) == 3.0)

    x1.set_vector([1.0, 2.0, 3.0, 4.0])

    x1_unary_sub = (-x1)
    assert (x1_unary_sub(0, 0) == -1.0)
    assert (x1_unary_sub(0, 1) == -2.0)
    assert (x1_unary_sub(1, 0) == -3.0)
    assert (x1_unary_sub(1, 1) == -4.0)

def test_glr_get_from_vector():
    """
    * Tests the get/from_vector operation.
    """

    from lielab.domain import glr

    x0 = glr.from_vector([])
    x0bar = x0.get_vector()

    assert (x0.get_shape() == 0)
    assert (x0bar.size == 0)

    x1 = glr.from_vector([1.0])
    x1bar = x1.get_vector()

    assert (x1.get_shape() == 1)
    assert (x1bar.size == 1)
    assert (x1bar[0] == 1.0)

    x2 = glr.from_vector([1.0, 2.0])
    x2bar = x2.get_vector()

    assert (x2.get_shape() == 2)
    assert (x2bar.size == 4)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 0.0)
    assert (x2bar[3] == 0.0)

    x3 = glr.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.get_vector()

    assert (x3.get_shape() == 2)
    assert (x3bar.size == 4)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)
    assert (x3bar[3] == 0.0)

    x5 = glr.from_vector([1.0, 2.0, 3.0, 4.0, 5.0])
    x5bar = x5.get_vector()

    assert (x5.get_shape() == 3)
    assert (x5bar.size == 9)
    assert (x5bar[0] == 1.0)
    assert (x5bar[1] == 2.0)
    assert (x5bar[2] == 3.0)
    assert (x5bar[3] == 4.0)
    assert (x5bar[4] == 5.0)
    assert (x5bar[5] == 0.0)
    assert (x5bar[6] == 0.0)
    assert (x5bar[7] == 0.0)
    assert (x5bar[8] == 0.0)

def test_glr_project():
    from lielab.domain import glr

    rand_2_2 = np.random.rand(2, 2)
    proj_2_2 = glr.project(rand_2_2)

    assert (proj_2_2.shape[0] == 2)
    assert (proj_2_2.shape[1] == 2)
    assert (proj_2_2[0, 0] == rand_2_2[0, 0])
    assert (proj_2_2[0, 1] == rand_2_2[0, 1])
    assert (proj_2_2[1, 0] == rand_2_2[1, 0])
    assert (proj_2_2[1, 1] == rand_2_2[1, 1])

    rand_3_3 = np.random.rand(3, 3)
    proj_3_3 = glr.project(rand_3_3)

    assert (proj_3_3.shape[0] == 3)
    assert (proj_3_3.shape[1] == 3)
    assert (proj_3_3[0, 0] == rand_3_3[0, 0])
    assert (proj_3_3[0, 1] == rand_3_3[0, 1])
    assert (proj_3_3[0, 2] == rand_3_3[0, 2])
    assert (proj_3_3[1, 0] == rand_3_3[1, 0])
    assert (proj_3_3[1, 1] == rand_3_3[1, 1])
    assert (proj_3_3[1, 2] == rand_3_3[1, 2])
    assert (proj_3_3[2, 0] == rand_3_3[2, 0])
    assert (proj_3_3[2, 1] == rand_3_3[2, 1])
    assert (proj_3_3[2, 2] == rand_3_3[2, 2])

    rand_2_3 = np.random.rand(2, 3)
    proj_2_3 = glr.project(rand_2_3)

    assert (proj_2_3.shape[0] == 2)
    assert (proj_2_3.shape[1] == 2)
    assert (proj_2_3[0, 0] == rand_2_3[0, 0])
    assert (proj_2_3[0, 1] == rand_2_3[0, 1])
    assert (proj_2_3[1, 0] == rand_2_3[1, 0])
    assert (proj_2_3[1, 1] == rand_2_3[1, 1])

    rand_3_2 = np.random.rand(3, 2)
    proj_3_2 = glr.project(rand_3_2)

    assert (proj_3_2.shape[0] == 2)
    assert (proj_3_2.shape[1] == 2)
    assert (proj_3_2[0, 0] == rand_3_2[0, 0])
    assert (proj_3_2[0, 1] == rand_3_2[0, 1])
    assert (proj_3_2[1, 0] == rand_3_2[1, 0])
    assert (proj_3_2[1, 1] == rand_3_2[1, 1])
