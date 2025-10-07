from lielab.testing import *
import pytest

def test_sp_to_string():
    from lielab.domain import sp

    x0 = sp(0)
    assert (x0.to_string() == "sp(0, R)")
    x1 = sp(1)
    assert (x1.to_string() == "sp(0, R)")
    x2 = sp(2)
    assert (x2.to_string() == "sp(2, R)")
    x10 = sp(10)
    assert (x10.to_string() == "sp(10, R)")

def test_sp_main_initializer():
    from lielab.domain import sp

    xblank = sp()
    assert (xblank.get_dimension() == 0)

    x0 = sp(0)
    assert (x0.get_dimension() == 0)
    x1 = sp(1)
    assert (x1.get_dimension() == 0)
    x3 = sp(3)
    assert (x3.get_dimension() == 3)

def test_sp_matrix_initializer():
    from lielab.domain import sp

    x0 = sp(np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    with pytest.raises(RuntimeError):
        sp(np.random.rand(1, 1))

    x2 = sp(np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        sp(np.random.rand(2, 4))

    with pytest.raises(RuntimeError):
        sp(np.random.rand(4, 2))

def test_sp_basis_initializer():
    from lielab.domain import sp

    xm10 = sp.basis(-1, 0)
    assert (xm10.get_dimension() == 0)
    xm10bar = xm10.get_vector()
    assert (xm10bar.size == 0)

    x00 = sp.basis(0, 0)
    assert (x00.get_dimension() == 0)
    x00bar = x00.get_vector()
    assert (x00bar.size == 0)

    x10 = sp.basis(1, 0)
    assert (x10.get_dimension() == 0)
    x10bar = x10.get_vector()
    assert (x10bar.size == 0)

    x02 = sp.basis(0, 2)
    assert (x02.get_dimension() == 3)
    x02bar = x02.get_vector()
    assert (x02bar.size == 3)
    assert (x02bar[0] == 1.0)
    assert (x02bar[1] == 0.0)
    assert (x02bar[2] == 0.0)

    x12 = sp.basis(1, 2)
    assert (x12.get_dimension() == 3)
    x12bar = x12.get_vector()
    assert (x12bar.size == 3)
    assert (x12bar[0] == 0.0)
    assert (x12bar[1] == 1.0)
    assert (x12bar[2] == 0.0)

    x22 = sp.basis(2, 2)
    assert (x22.get_dimension() == 3)
    x22bar = x22.get_vector()
    assert (x22bar.size == 3)
    assert (x22bar[0] == 0.0)
    assert (x22bar[1] == 0.0)
    assert (x22bar[2] == 1.0)

    x04 = sp.basis(0, 4)
    assert (x04.get_dimension() == 10)
    x04bar = x04.get_vector()
    assert (x04bar.size == 10)
    assert (x04bar[0] == 1.0)
    assert (x04bar[1] == 0.0)
    assert (x04bar[2] == 0.0)
    assert (x04bar[3] == 0.0)
    assert (x04bar[4] == 0.0)
    assert (x04bar[5] == 0.0)
    assert (x04bar[6] == 0.0)
    assert (x04bar[7] == 0.0)
    assert (x04bar[8] == 0.0)
    assert (x04bar[9] == 0.0)

def test_sp_from_shape_initializer():
    from lielab.domain import sp

    x0 = sp.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = sp.from_shape(1)
    assert (x1.get_dimension() == 0)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 0)
    assert (x1hat.shape[1] == 0)

    x2 = sp.from_shape(2)
    assert (x2.get_dimension() == 3)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_sp_get_dimension():
    from lielab.domain import sp

    zero = sp(0)
    two = sp(2)
    four = sp(4)
    six = sp(6)
    eight = sp(8)

    # Dimensions
    assert (zero.get_dimension() == 0)
    assert (two.get_dimension() == 3)
    assert (four.get_dimension() == 10)
    assert (six.get_dimension() == 21)
    assert (eight.get_dimension() == 36)

def test_sp_set_get_vector():
    """
    * Tests the set/get_vector operation.
    """

    from lielab.domain import sp

    x0 = sp(0)
    x0.set_vector([])
    x0bar = x0.get_vector()

    assert (x0bar.size == 0)

    x2 = sp(2)
    x2.set_vector([1.0, 2.0, 3.0, 4.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 3)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 3.0)

    x2.set_vector([5.0, 6.0, 7.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 3)
    assert (x2bar[0] == 5.0)
    assert (x2bar[1] == 6.0)
    assert (x2bar[2] == 7.0)

    x2.set_vector([8.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 3)
    assert (x2bar[0] == 8.0)
    assert (x2bar[1] == 6.0)
    assert (x2bar[2] == 7.0)

    x4 = sp(4)
    x4.set_vector([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
    x4bar = x4.get_vector()

    assert (x4bar.size == 10)
    assert (x4bar[0] == 1.0)
    assert (x4bar[1] == 2.0)
    assert (x4bar[2] == 3.0)
    assert (x4bar[3] == 4.0)
    assert (x4bar[4] == 5.0)
    assert (x4bar[5] == 6.0)
    assert (x4bar[6] == 7.0)
    assert (x4bar[7] == 8.0)
    assert (x4bar[8] == 9.0)
    assert (x4bar[9] == 0.0)

    x4.set_vector([10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0])
    x4bar = x4.get_vector()

    assert (x4bar.size == 10)
    assert (x4bar[0] == 10.0)
    assert (x4bar[1] == 11.0)
    assert (x4bar[2] == 12.0)
    assert (x4bar[3] == 13.0)
    assert (x4bar[4] == 14.0)
    assert (x4bar[5] == 15.0)
    assert (x4bar[6] == 16.0)
    assert (x4bar[7] == 17.0)
    assert (x4bar[8] == 18.0)
    assert (x4bar[9] == 19.0)

    x4.set_vector([20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0])
    x4bar = x4.get_vector()

    assert (x4bar.size == 10)
    assert (x4bar[0] == 20.0)
    assert (x4bar[1] == 21.0)
    assert (x4bar[2] == 22.0)
    assert (x4bar[3] == 23.0)
    assert (x4bar[4] == 24.0)
    assert (x4bar[5] == 25.0)
    assert (x4bar[6] == 26.0)
    assert (x4bar[7] == 27.0)
    assert (x4bar[8] == 28.0)
    assert (x4bar[9] == 29.0)

def test_sp_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import sp

    x0 = sp(0)
    x0.set_vector([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x2 = sp(2)
    x2.set_vector([1.0, 2.0, 3.0, 4.0])
    x2hat = x2.get_matrix()

    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 1.0)
    assert (x2hat[0, 1] == 2.0)
    assert (x2hat[1, 0] == 3.0)
    assert (x2hat[1, 1] == -1.0)

    x2.set_vector([5.0, 6.0, 7.0])
    x2hat = x2.get_matrix()

    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 5.0)
    assert (x2hat[0, 1] == 6.0)
    assert (x2hat[1, 0] == 7.0)
    assert (x2hat[1, 1] == -5.0)

    x2.set_vector([8.0, 9.0])
    x2hat = x2.get_matrix()

    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 8.0)
    assert (x2hat[0, 1] == 9.0)
    assert (x2hat[1, 0] == 7.0)
    assert (x2hat[1, 1] == -8.0)

def test_sp_operator_parenthesis():
    from lielab.domain import sp

    x0 = sp(0)
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

    x2 = sp(2)
    x2.set_vector([1.0, 2.0, 3.0])

    # In bounds
    assert (x2(0) == 1.0)
    assert (x2(1) == 2.0)
    assert (x2(2) == 3.0)
    assert (x2(-1) == 3.0)
    assert (x2(-2) == 2.0)
    assert (x2(-3) == 1.0)

    # Out of bounds
    assert (np.isnan(x2(-4)))
    assert (np.isnan(x2(3)))

    # In bounds
    assert (x2(0, 0) == 1.0)
    assert (x2(0, 1) == 2.0)
    assert (x2(1, 0) == 3.0)
    assert (x2(1, 1) == -1.0)
    assert (x2(-1, -1) == -1.0)
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

# TODO: Math ops int

def test_sp_math_ops_double():
    from lielab.domain import sp

    x1 = sp(2)
    x1.set_vector([1.25, 2.5, 3.75])

    x1_lm_2 = 2.0*x1
    assert (x1_lm_2(0) == 2.5)
    assert (x1_lm_2(1) == 5.0)
    assert (x1_lm_2(2) == 7.5)

    x1_rm_2 = x1*2.0
    assert (x1_rm_2(0) == 2.5)
    assert (x1_rm_2(1) == 5.0)
    assert (x1_rm_2(2) == 7.5)

    x1 *= 2.0
    assert (x1(0) == 2.5)
    assert (x1(1) == 5.0)
    assert (x1(2) == 7.5)

    x1.set_vector([1.25, 2.5, 3.75])

    x1_d_2 = x1/2.0
    assert (x1_d_2(0) == 0.625)
    assert (x1_d_2(1) == 1.25)
    assert (x1_d_2(2) == 1.875)

    x1 /= 2.0
    assert (x1(0) == 0.625)
    assert (x1(1) == 1.25)
    assert (x1(2) == 1.875)

def test_sp_math_ops_sp():
    from lielab.domain import sp

    x1 = sp(2)
    x2 = sp(2)
    x1.set_vector([1.0, 2.0, 3.0])
    x2.set_vector([1.25, 2.5, 3.75])

    x1_add_x2 = x1 + x2
    assert (x1_add_x2(0) == 2.25)
    assert (x1_add_x2(1) == 4.5)
    assert (x1_add_x2(2) == 6.75)

    x1 += x2
    assert (x1(0) == 2.25)
    assert (x1(1) == 4.5)
    assert (x1(2) == 6.75)

    x1.set_vector([1.0, 2.0, 3.0])

    x1_sub_x2 = x1 - x2
    assert (x1_sub_x2(0) == -0.25)
    assert (x1_sub_x2(1) == -0.5)
    assert (x1_sub_x2(2) == -0.75)

    x1 -= x2
    assert (x1(0) == -0.25)
    assert (x1(1) == -0.5)
    assert (x1(2) == -0.75)

    x1.set_vector([1.0, 2.0, 3.0])

    x1_unary_sub = (-x1)
    assert (x1_unary_sub(0) == -1.0)
    assert (x1_unary_sub(1) == -2.0)
    assert (x1_unary_sub(2) == -3.0)

def test_sp_from_vector():
    """
    * Tests the from_vector operation.
    """

    from lielab.domain import sp

    x0 = sp.from_vector([])
    x0bar = x0.get_vector()

    assert (x0.get_shape() == 0)
    assert (x0bar.size == 0)

    x1 = sp.from_vector([1.0])
    x1bar = x1.get_vector()

    assert (x1.get_shape() == 2)
    assert (x1bar.size == 3)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 0.0)
    assert (x1bar[2] == 0.0)

    x2 = sp.from_vector([1.0, 2.0])
    x2bar = x2.get_vector()

    assert (x2.get_shape() == 2)
    assert (x2bar.size == 3)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 0.0)

    x3 = sp.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.get_vector()

    assert (x3.get_shape() == 2)
    assert (x3bar.size == 3)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)

    x4 = sp.from_vector([1.0, 2.0, 3.0, 4.0])
    x4bar = x4.get_vector()

    assert (x4.get_shape() == 4)
    assert (x4bar.size == 10)
    assert (x4bar[0] == 1.0)
    assert (x4bar[1] == 2.0)
    assert (x4bar[2] == 3.0)
    assert (x4bar[3] == 4.0)
    assert (x4bar[4] == 0.0)
    assert (x4bar[5] == 0.0)
    assert (x4bar[6] == 0.0)
    assert (x4bar[7] == 0.0)
    assert (x4bar[8] == 0.0)
    assert (x4bar[9] == 0.0)

def test_sp_project():
    from lielab.domain import sp

    rand_2_2 = np.random.rand(2, 2)
    proj_2_2 = sp.project(rand_2_2)

    assert (proj_2_2.shape[0] == 2)
    assert (proj_2_2.shape[1] == 2)
    assert (proj_2_2[0, 1] == rand_2_2[0, 1])
    assert (proj_2_2[1, 0] == rand_2_2[1, 0])
    assert (proj_2_2.trace() == 0.0)

    rand_2_3 = np.random.rand(2, 3)
    proj_2_3 = sp.project(rand_2_3)

    assert (proj_2_3.shape[0] == 2)
    assert (proj_2_3.shape[1] == 2)
    assert (proj_2_3[0, 1] == rand_2_3[0, 1])
    assert (proj_2_3[1, 0] == rand_2_3[1, 0])
    assert (proj_2_3.trace() == 0.0)

    rand_3_2 = np.random.rand(3, 2)
    proj_3_2 = sp.project(rand_3_2)

    assert (proj_3_2.shape[0] == 2)
    assert (proj_3_2.shape[1] == 2)
    assert (proj_3_2[0, 1] == rand_3_2[0, 1])
    assert (proj_3_2[1, 0] == rand_3_2[1, 0])
    assert (proj_3_2.trace() == 0.0)

    rand_3_3 = np.random.rand(3, 3)
    proj_3_3 = sp.project(rand_3_3)

    assert (proj_3_3.shape[0] == 2)
    assert (proj_3_3.shape[1] == 2)
    assert (proj_3_3[0, 1] == rand_3_3[0, 1])
    assert (proj_3_3[1, 0] == rand_3_3[1, 0])
    assert (proj_3_3.trace() == 0.0)
