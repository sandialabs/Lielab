from lielab.testing import *
import pytest

def test_so_to_string():
    from lielab.domain import so

    x0 = so(0)
    assert (x0.to_string() == "so(0)")
    x1 = so(1)
    assert (x1.to_string() == "so(1)")
    x10 = so(10)
    assert (x10.to_string() == "so(10)")

def test_so_main_initializer():
    from lielab.domain import so

    xblank = so()
    assert (xblank.get_dimension() == 0)

    x0 = so(0)
    assert (x0.get_dimension() == 0)
    x1 = so(1)
    assert (x1.get_dimension() == 0)
    x3 = so(3)
    assert (x3.get_dimension() == 3)

def test_so_matrix_initializer():
    from lielab.domain import so

    x0 = so(np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = so(np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = so(np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        so(np.random.rand(2, 3))
    
    with pytest.raises(RuntimeError):
        so(np.random.rand(3, 2))

def test_so_basis_initializer():
    from lielab.domain import so

    xm10 = so.basis(-1, 0)
    assert (xm10.get_dimension() == 0)
    xm10bar = xm10.get_vector()
    assert (xm10bar.size == 0)

    x00 = so.basis(0, 0)
    assert (x00.get_dimension() == 0)
    x00bar = x00.get_vector()
    assert (x00bar.size == 0)

    x10 = so.basis(1, 0)
    assert (x10.get_dimension() == 0)
    x10bar = x10.get_vector()
    assert (x10bar.size == 0)

    x01 = so.basis(0, 2)
    assert (x01.get_dimension() == 1)
    x01bar = x01.get_vector()
    assert (x01bar.size == 1)
    assert (x01bar[0] == 1.0)

    x11 = so.basis(1, 2)
    assert (x11.get_dimension() == 1)
    x11bar = x11.get_vector()
    assert (x11bar.size == 1)
    assert (x11bar[0] == 0.0)

    x21 = so.basis(2, 2)
    assert (x21.get_dimension() == 1)
    x21bar = x21.get_vector()
    assert (x21bar.size == 1)
    assert (x21bar[0] == 0.0)

    x02 = so.basis(0, 3)
    assert (x02.get_dimension() == 3)
    x02bar = x02.get_vector()
    assert (x02bar.size == 3)
    assert (x02bar[0] == 1.0)
    assert (x02bar[1] == 0.0)
    assert (x02bar[2] == 0.0)

def test_so_from_shape_initializer():
    from lielab.domain import so

    x0 = so.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = so.from_shape(1)
    assert (x1.get_dimension() == 0)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = so.from_shape(2)
    assert (x2.get_dimension() == 1)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_so_get_dimension():
    from lielab.domain import so

    zero = so(0)
    one = so(1)
    two = so(2)
    three = so(3)
    four = so(4)
    five = so(5)
    six = so(6)
    seven = so(7)
    eight = so(8)

    # Dimensions
    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 0)
    assert (two.get_dimension() == 1)
    assert (three.get_dimension() == 3)
    assert (four.get_dimension() == 6)
    assert (five.get_dimension() == 10)
    assert (six.get_dimension() == 15)
    assert (seven.get_dimension() == 21)
    assert (eight.get_dimension() == 28)

def test_so_set_get_vector():
    """
    * Tests the set/get_vector operation.
    """

    from lielab.domain import so

    x0 = so(0)
    x0.set_vector([])
    x0bar = x0.get_vector()

    assert (x0bar.size == 0)

    x3 = so(3)
    x3.set_vector([1.0, 2.0, 3.0, 4.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 3)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)

    x3.set_vector([5.0, 6.0, 7.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 3)
    assert (x3bar[0] == 5.0)
    assert (x3bar[1] == 6.0)
    assert (x3bar[2] == 7.0)

    x3.set_vector([8.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 3)
    assert (x3bar[0] == 8.0)
    assert (x3bar[1] == 6.0)
    assert (x3bar[2] == 7.0)

    x4 = so(4)
    x4.set_vector([1.0, 2.0, 3.0, 4.0, 5.0])
    x4bar = x4.get_vector()

    assert (x4bar.size == 6)
    assert (x4bar[0] == 1.0)
    assert (x4bar[1] == 2.0)
    assert (x4bar[2] == 3.0)
    assert (x4bar[3] == 4.0)
    assert (x4bar[4] == 5.0)
    assert (x4bar[5] == 0.0)

    x4.set_vector([7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    x4bar = x4.get_vector()

    assert (x4bar.size == 6)
    assert (x4bar[0] == 7.0)
    assert (x4bar[1] == 8.0)
    assert (x4bar[2] == 9.0)
    assert (x4bar[3] == 10.0)
    assert (x4bar[4] == 11.0)
    assert (x4bar[5] == 12.0)

    x4.set_vector([13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0])
    x4bar = x4.get_vector()

    assert (x4bar.size == 6)
    assert (x4bar[0] == 13.0)
    assert (x4bar[1] == 14.0)
    assert (x4bar[2] == 15.0)
    assert (x4bar[3] == 16.0)
    assert (x4bar[4] == 17.0)
    assert (x4bar[5] == 18.0)

def test_so_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import so

    x0 = so(0)
    x0.set_vector([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = so(1)
    x1.set_vector([])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == 0.0)

    x2 = so(2)
    x2.set_vector([1.0, 2.0])
    x2hat = x2.get_matrix()

    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 0.0)
    assert (x2hat[0, 1] == -1.0)
    assert (x2hat[1, 0] == 1.0)
    assert (x2hat[1, 1] == 0.0)

    x2.set_vector([3.0])
    x2hat = x2.get_matrix()

    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 0.0)
    assert (x2hat[0, 1] == -3.0)
    assert (x2hat[1, 0] == 3.0)
    assert (x2hat[1, 1] == 0.0)

    x2.set_vector([])
    x2hat = x2.get_matrix()

    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 0.0)
    assert (x2hat[0, 1] == -3.0)
    assert (x2hat[1, 0] == 3.0)
    assert (x2hat[1, 1] == 0.0)

    x3 = so(3)
    x3.set_vector([1.0])
    x3hat = x3.get_matrix()
    assert (x3hat.shape[0] == 3)
    assert (x3hat.shape[1] == 3)
    assert (x3hat[0, 0] == 0.0)
    assert (x3hat[0, 1] == 0.0)
    assert (x3hat[0, 2] == 0.0)
    assert (x3hat[1, 0] == 0.0)
    assert (x3hat[1, 1] == 0.0)
    assert (x3hat[1, 2] == -1.0)
    assert (x3hat[2, 0] == 0.0)
    assert (x3hat[2, 1] == 1.0)
    assert (x3hat[2, 2] == 0.0)

    x3.set_vector([2.0, 3.0])
    x3hat = x3.get_matrix()
    assert (x3hat.shape[0] == 3)
    assert (x3hat.shape[1] == 3)
    assert (x3hat[0, 0] == 0.0)
    assert (x3hat[0, 1] == 0.0)
    assert (x3hat[0, 2] == 3.0)
    assert (x3hat[1, 0] == 0.0)
    assert (x3hat[1, 1] == 0.0)
    assert (x3hat[1, 2] == -2.0)
    assert (x3hat[2, 0] == -3.0)
    assert (x3hat[2, 1] == 2.0)
    assert (x3hat[2, 2] == 0.0)

    x3.set_vector([4.0, 5.0, 6.0])
    x3hat = x3.get_matrix()
    assert (x3hat.shape[0] == 3)
    assert (x3hat.shape[1] == 3)
    assert (x3hat[0, 0] == 0.0)
    assert (x3hat[0, 1] == -6.0)
    assert (x3hat[0, 2] == 5.0)
    assert (x3hat[1, 0] == 6.0)
    assert (x3hat[1, 1] == 0.0)
    assert (x3hat[1, 2] == -4.0)
    assert (x3hat[2, 0] == -5.0)
    assert (x3hat[2, 1] == 4.0)
    assert (x3hat[2, 2] == 0.0)

def test_so_operator_parenthesis():
    from lielab.domain import so

    x0 = so(0)
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

    x1 = so(1)
    x1.set_vector([])

    # Out of bounds
    assert (np.isnan(x1(-1)))
    assert (np.isnan(x1(0)))
    assert (np.isnan(x1(1)))

    # In bounds
    assert (x1(0, 0) == 0.0)
    assert (x1(-1, -1) == 0.0)

    # Out of bounds
    assert (np.isnan(x1(0, -2)))
    assert (np.isnan(x1(-2, 0)))
    assert (np.isnan(x1(-2, -2)))
    assert (np.isnan(x1(0, 1)))
    assert (np.isnan(x1(1, 0)))
    assert (np.isnan(x1(1, 1)))

    x3 = so(3)
    x3.set_vector([1.0, 2.0, 3.0])

    # In bounds
    assert (x3(0) == 1.0)
    assert (x3(1) == 2.0)
    assert (x3(2) == 3.0)
    assert (x3(-1) == 3.0)
    assert (x3(-2) == 2.0)
    assert (x3(-3) == 1.0)

    # Out of bounds
    assert (np.isnan(x3(-4)))
    assert (np.isnan(x3(3)))

    # In bounds
    assert (x3(0, 0) == 0.0)
    assert (x3(0, 1) == -3.0)
    assert (x3(0, 2) == 2.0)
    assert (x3(1, 0) == 3.0)
    assert (x3(1, 1) == 0.0)
    assert (x3(1, 2) == -1.0)
    assert (x3(2, 0) == -2.0)
    assert (x3(2, 1) == 1.0)
    assert (x3(2, 2) == 0.0)
    assert (x3(-1, -1) == 0.0)
    assert (x3(-1, -2) == 1.0)
    assert (x3(-1, -3) == -2.0)
    assert (x3(-2, -1) == -1.0)
    assert (x3(-2, -2) == 0.0)
    assert (x3(-2, -3) == 3.0)
    assert (x3(-3, -1) == 2.0)
    assert (x3(-3, -2) == -3.0)
    assert (x3(-3, -3) == 0.0)

    # Out of bounds
    assert (np.isnan(x3(0, -4)))
    assert (np.isnan(x3(-4, 0)))
    assert (np.isnan(x3(-4, -4)))
    assert (np.isnan(x3(0, 3)))
    assert (np.isnan(x3(3, 0)))
    assert (np.isnan(x3(3, 3)))

# TODO: Math ops int

def test_so_math_ops_double():
    from lielab.domain import so

    x1 = so(3)
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

def test_so_math_ops_so():
    from lielab.domain import so

    x1 = so(3)
    x2 = so(3)
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

def test_so_from_vector():
    """
    * Tests the from_vector operation.
    """

    from lielab.domain import so

    x0 = so.from_vector([])
    x0bar = x0.get_vector()

    assert (x0.get_shape() == 1)
    assert (x0bar.size == 0)

    x1 = so.from_vector([1.0])
    x1bar = x1.get_vector()

    assert (x1.get_shape() == 2)
    assert (x1bar.size == 1)
    assert (x1bar[0] == 1.0)

    x2 = so.from_vector([1.0, 2.0])
    x2bar = x2.get_vector()

    assert (x2.get_shape() == 3)
    assert (x2bar.size == 3)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 0.0)

    x3 = so.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.get_vector()

    assert (x3.get_shape() == 3)
    assert (x3bar.size == 3)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)

def test_so_project():
    from lielab.domain import so

    rand_2_2 = np.random.rand(2, 2)
    proj_2_2 = so.project(rand_2_2)

    assert (proj_2_2.shape[0] == 2)
    assert (proj_2_2.shape[1] == 2)
    assert (proj_2_2[0, 0] == 0.0)
    assert (proj_2_2[0, 1] == -proj_2_2[1, 0])
    assert (proj_2_2[1, 0] == -proj_2_2[0, 1])
    assert (proj_2_2[1, 1] == 0.0)

    rand_3_3 = np.random.rand(3, 3)
    proj_3_3 = so.project(rand_3_3)

    assert (proj_3_3.shape[0] == 3)
    assert (proj_3_3.shape[1] == 3)
    assert (proj_3_3[0, 0] == 0.0)
    assert (proj_3_3[0, 1] == -proj_3_3[1, 0])
    assert (proj_3_3[0, 2] == -proj_3_3[2, 0])
    assert (proj_3_3[1, 0] == -proj_3_3[0, 1])
    assert (proj_3_3[1, 1] == 0.0)
    assert (proj_3_3[1, 2] == -proj_3_3[2, 1])
    assert (proj_3_3[2, 0] == -proj_3_3[0, 2])
    assert (proj_3_3[2, 1] == -proj_3_3[1, 2])
    assert (proj_3_3[2, 2] == 0.0)

    rand_2_3 = np.random.rand(2, 3)
    proj_2_3 = so.project(rand_2_3)

    assert (proj_2_3.shape[0] == 2)
    assert (proj_2_3.shape[1] == 2)
    assert (proj_2_3[0, 0] == 0.0)
    assert (proj_2_3[0, 1] == -proj_2_3[1, 0])
    assert (proj_2_3[1, 0] == -proj_2_3[0, 1])
    assert (proj_2_3[1, 1] == 0.0)

    rand_3_2 = np.random.rand(3, 2)
    proj_3_2 = so.project(rand_3_2)

    assert (proj_3_2.shape[0] == 2)
    assert (proj_3_2.shape[1] == 2)
    assert (proj_3_2[0, 0] == 0.0)
    assert (proj_3_2[0, 1] == -proj_3_2[1, 0])
    assert (proj_3_2[1, 0] == -proj_3_2[0, 1])
    assert (proj_3_2[1, 1] == 0.0)

def test_so2():
    """
    Tests the algebra with so(2).
    
    Note that some of these test cases are trivial since so(2) is 1-dimensional.
    """

    from lielab.domain import so
    from lielab.functions import commutator

    x = so.basis(0,2)
    zero = x*0

    assert_domain(commutator(x, x), zero)

def test_so3():
    """
    * Tests the algebra with so(3).
    """

    from lielab.domain import so
    from lielab.functions import commutator

    x = so.basis(0,3)
    y = so.basis(1,3)
    z = so.basis(2,3)
    zero = x*0

    assert_domain(commutator(x, y), z)
    assert_domain(commutator(y, z), x)
    assert_domain(commutator(z, x), y)
    assert_domain(commutator(y, x), -z)
    assert_domain(commutator(z, y), -x)
    assert_domain(commutator(x, z), -y)
