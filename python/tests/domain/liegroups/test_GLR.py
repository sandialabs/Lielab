from lielab.testing import *
import pytest

def test_GLR_to_string():
    from lielab.domain import GLR

    x0 = GLR(0)
    assert (x0.to_string() == "GL(0, R)")
    x1 = GLR(1)
    assert (x1.to_string() == "GL(1, R)")
    x10 = GLR(10)
    assert (x10.to_string() == "GL(10, R)")

def test_GLR_main_initializer():
    from lielab.domain import GLR

    xblank = GLR()
    assert (xblank.get_dimension() == 0)

    x0 = GLR(0)
    assert (x0.get_dimension() == 0)
    x1 = GLR(1)
    assert (x1.get_dimension() == 1)
    x10 = GLR(10)
    assert (x10.get_dimension() == 100)

def test_GLR_matrix_initializer():
    from lielab.domain import GLR

    x0 = GLR(np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = GLR(np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = GLR(np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        GLR(np.random.rand(2, 3))
    
    with pytest.raises(RuntimeError):
        GLR(np.random.rand(3, 2))

def test_GLR_from_shape_initializer():
    from lielab.domain import GLR

    x0 = GLR.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = GLR.from_shape(1)
    assert (x1.get_dimension() == 1)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = GLR.from_shape(2)
    assert (x2.get_dimension() == 4)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_GLR_get_dimension():
    from lielab.domain import GLR

    zero = GLR(0)
    one = GLR(1)
    two = GLR(2)
    three = GLR(3)
    four = GLR(4)
    five = GLR(5)
    six = GLR(6)
    seven = GLR(7)
    eight = GLR(8)

    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 1)
    assert (two.get_dimension() == 4)
    assert (three.get_dimension() == 9)
    assert (four.get_dimension() == 16)
    assert (five.get_dimension() == 25)
    assert (six.get_dimension() == 36)
    assert (seven.get_dimension() == 49)
    assert (eight.get_dimension() == 64)

def test_GLR_serialize_unserialize():
    """
    * Tests the serialize/unserialize operation.
    """

    from lielab.domain import GLR

    x0 = GLR(0)
    x0.unserialize([])
    x0bar = x0.serialize()

    assert (x0bar.size == 0)

    x1 = GLR(1)
    x1.unserialize([1.0, 2.0])
    x1bar = x1.serialize()

    assert (x1bar.size == 1)
    assert (x1bar[0] == 1.0)

    x1.unserialize([3.0])
    x1bar = x1.serialize()

    assert (x1bar.size == 1)
    assert (x1bar[0] == 3.0)

    x1.unserialize([])
    x1bar = x1.serialize()

    assert (x1bar.size == 1)
    assert (x1bar[0] == 3.0)

    x2 = GLR(2)
    x2.unserialize([1.0, 2.0, 3.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 3.0)
    assert (x2bar[3] == 1.0)

    x2.unserialize([4.0, 5.0, 6.0, 7.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 4.0)
    assert (x2bar[1] == 5.0)
    assert (x2bar[2] == 6.0)
    assert (x2bar[3] == 7.0)

    x2.unserialize([8.0, 9.0, 10.0, 11.0, 12.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 8.0)
    assert (x2bar[1] == 9.0)
    assert (x2bar[2] == 10.0)
    assert (x2bar[3] == 11.0)

def test_GLR_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import GLR

    x0 = GLR(0)
    x0.unserialize([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = GLR(1)
    x1.unserialize([])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == 1.0)

    x1.unserialize([1.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == 1.0)

    x1.unserialize([2.0, 3.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == 2.0)

    x2 = GLR(2)
    x2.unserialize([1.0, 2.0, 3.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 1.0)
    assert (x2hat[0, 1] == 2.0)
    assert (x2hat[1, 0] == 3.0)
    assert (x2hat[1, 1] == 1.0)

    x2.unserialize([4.0, 5.0, 6.0, 7.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 4.0)
    assert (x2hat[0, 1] == 5.0)
    assert (x2hat[1, 0] == 6.0)
    assert (x2hat[1, 1] == 7.0)

    x2.unserialize([8.0, 9.0, 10.0, 11.0, 12.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == 8.0)
    assert (x2hat[0, 1] == 9.0)
    assert (x2hat[1, 0] == 10.0)
    assert (x2hat[1, 1] == 11.0)

def test_GLR_operator_parenthesis():
    from lielab.domain import GLR

    x0 = GLR.from_shape(0)
    x0.unserialize([])

    # Out of bounds
    # assert (np.isnan(x0(-1)))
    # assert (np.isnan(x0(0)))
    # assert (np.isnan(x0(1)))

    # Out of bounds
    assert (np.isnan(x0(0, -1)))
    assert (np.isnan(x0(-1, 0)))
    assert (np.isnan(x0(-1, -1)))
    assert (np.isnan(x0(0, 0)))
    assert (np.isnan(x0(0, 1)))
    assert (np.isnan(x0(1, 0)))
    assert (np.isnan(x0(1, 1)))

    x1 = GLR(1)
    x1.unserialize([1.0])

    # In bounds
    # assert (x1(0) == 1.0)

    # Out of bounds
    # assert (np.isnan(x1(-1)))
    # assert (np.isnan(x1(1)))

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

    x2 = GLR(2)
    x2.unserialize([1.0, 2.0, 3.0, 4.0])

    # In bounds
    # assert (x2(0) == 1.0)
    # assert (x2(1) == 2.0)
    # assert (x2(2) == 3.0)
    # assert (x2(3) == 4.0)

    # Out of bounds
    # assert (np.isnan(x2(-1)))
    # assert (np.isnan(x2(4)))

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

def test_GLR_math_ops_GLR():
    from lielab.domain import GLR

    x1 = GLR(2)
    x2 = GLR(2)
    x1.unserialize([1.0, 2.0, 3.0, 4.0])
    x2.unserialize([1.25, 2.5, 3.75, 1.0])

    x1_prod_x2 = x1*x2
    x1_prod_x2bar = x1_prod_x2.serialize()
    assert (x1_prod_x2bar.size == 4)
    assert (x1_prod_x2bar[0] == 8.75)
    assert (x1_prod_x2bar[1] == 4.5)
    assert (x1_prod_x2bar[2] == 18.75)
    assert (x1_prod_x2bar[3] == 11.5)

    x1 *= x2
    x1bar = x1.serialize()
    assert (x1bar.size == 4)
    assert (x1bar[0] == 8.75)
    assert (x1bar[1] == 4.5)
    assert (x1bar[2] == 18.75)
    assert (x1bar[3] == 11.5)

    x1.unserialize([1.0, 2.0, 3.0, 4.0])

    x1_inv = x1.inverse()
    x1_invbar = x1_inv.serialize()
    assert (x1_invbar.size == 4)
    assert (np.abs(x1_invbar[0] - -2.0) <= 1e-14)
    assert (np.abs(x1_invbar[1] - 1.0) <= 1e-14)
    assert (np.abs(x1_invbar[2] - 1.5) <= 1e-14)
    assert (np.abs(x1_invbar[3] - -0.5) <= 1e-14)

def test_GLR_project():
    from lielab.domain import GLR

    rand_2_2 = np.random.rand(2, 2)
    proj_2_2 = GLR.project(rand_2_2)

    assert (proj_2_2.shape[0] == 2)
    assert (proj_2_2.shape[1] == 2)
    assert (proj_2_2[0, 0] == rand_2_2[0, 0])
    assert (proj_2_2[0, 1] == rand_2_2[0, 1])
    assert (proj_2_2[1, 0] == rand_2_2[1, 0])
    assert (proj_2_2[1, 1] == rand_2_2[1, 1])

    rand_3_3 = np.random.rand(3, 3)
    proj_3_3 = GLR.project(rand_3_3)

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
    proj_2_3 = GLR.project(rand_2_3)

    assert (proj_2_3.shape[0] == 2)
    assert (proj_2_3.shape[1] == 2)
    assert (proj_2_3[0, 0] == rand_2_3[0, 0])
    assert (proj_2_3[0, 1] == rand_2_3[0, 1])
    assert (proj_2_3[1, 0] == rand_2_3[1, 0])
    assert (proj_2_3[1, 1] == rand_2_3[1, 1])

    rand_3_2 = np.random.rand(3, 2)
    proj_3_2 = GLR.project(rand_3_2)

    assert (proj_3_2.shape[0] == 2)
    assert (proj_3_2.shape[1] == 2)
    assert (proj_3_2[0, 0] == rand_3_2[0, 0])
    assert (proj_3_2[0, 1] == rand_3_2[0, 1])
    assert (proj_3_2[1, 0] == rand_3_2[1, 0])
    assert (proj_3_2[1, 1] == rand_3_2[1, 1])
