from lielab.testing import *
import pytest

def complex(a,b):
    return a + b*1j

def test_CN_to_string():
    from lielab.domain import CN

    xzero = CN.from_shape(0)
    assert (xzero.to_string() == "C^nan")
    x0 = CN(0)
    assert (x0.to_string() == "C^0")
    x1 = CN(1)
    assert (x1.to_string() == "C^1")
    x10 = CN(10)
    assert (x10.to_string() == "C^10")

def test_CN_main_initializer():
    from lielab.domain import CN

    xblank = CN()
    assert (xblank.get_dimension() == 0)

    x0 = CN(0)
    assert (x0.get_dimension() == 0)
    x1 = CN(1)
    assert (x1.get_dimension() == 2)
    x10 = CN(10)
    assert (x10.get_dimension() == 20)

def test_CN_matrix_initializer():
    from lielab.domain import CN

    x0 = CN(np.random.rand(0, 0) + 1j*np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = CN(np.random.rand(1, 1) + 1j*np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = CN(np.random.rand(2, 2) + 1j*np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        CN(np.random.rand(2, 3) + 1j*np.random.rand(2, 3))
    
    with pytest.raises(RuntimeError):
        CN(np.random.rand(3, 2) + 1j*np.random.rand(3, 2))

def test_CN_from_shape_initializer():
    from lielab.domain import CN

    x0 = CN.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = CN.from_shape(1)
    assert (x1.get_dimension() == 0)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = CN.from_shape(2)
    assert (x2.get_dimension() == 2)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_CN_get_dimension():
    from lielab.domain import CN

    veryzero = CN.from_shape(0)
    zero = CN(0)
    one = CN(1)
    two = CN(2)
    three = CN(3)
    four = CN(4)
    five = CN(5)
    six = CN(6)
    seven = CN(7)
    eight = CN(8)

    assert (veryzero.get_dimension() == 0)
    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 2)
    assert (two.get_dimension() == 4)
    assert (three.get_dimension() == 6)
    assert (four.get_dimension() == 8)
    assert (five.get_dimension() == 10)
    assert (six.get_dimension() == 12)
    assert (seven.get_dimension() == 14)
    assert (eight.get_dimension() == 16)

def test_CN_serialize_unserialize():
    """
    * Tests the serialize/unserialize operation.
    """

    from lielab.domain import CN

    xzero = CN.from_shape(0)
    xzero.unserialize([])
    xzerobar = xzero.serialize()

    assert (xzerobar.size == 0)

    x0 = CN(0)
    x0.unserialize([])
    x0bar = x0.serialize()

    assert (x0bar.size == 0)

    x1 = CN(1)
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

    x2 = CN(2)
    x2.unserialize([1.0, 2.0, 3.0])
    x2bar = x2.serialize()

    assert (x2bar.size == 4)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 3.0)
    assert (x2bar[3] == 0.0)

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

def test_CN_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import CN

    xzero = CN.from_shape(0)
    xzero.unserialize([])
    xzerohat = xzero.get_matrix()

    assert (xzerohat.shape[0] == 0)
    assert (xzerohat.shape[1] == 0)

    x0 = CN(0)
    x0.unserialize([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 1)
    assert (x0hat.shape[1] == 1)
    assert (x0hat[0, 0] == complex(1.0, 0.0))

    x1 = CN(1)
    x1.unserialize([1.0, 2.0, 3.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == complex(1.0, 2.0))

    assert (x1hat[0, 0] == complex(1.0, 0.0))
    assert (x1hat[1, 0] == complex(0.0, 0.0))
    assert (x1hat[1, 1] == complex(1.0, 0.0))

    x1.unserialize([4.0, 5.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == complex(4.0, 5.0))

    assert (x1hat[0, 0] == complex(1.0, 0.0))
    assert (x1hat[1, 0] == complex(0.0, 0.0))
    assert (x1hat[1, 1] == complex(1.0, 0.0))

    x1.unserialize([6.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == complex(6.0, 5.0))

    assert (x1hat[0, 0] == complex(1.0, 0.0))
    assert (x1hat[1, 0] == complex(0.0, 0.0))
    assert (x1hat[1, 1] == complex(1.0, 0.0))

    x2 = CN(2)
    x2.unserialize([1.0, 2.0, 3.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == complex(1.0, 2.0))
    assert (x2hat[1, 2] == complex(3.0, 0.0))

    assert (x2hat[0, 0] == complex(1.0, 0.0))
    assert (x2hat[0, 1] == complex(0.0, 0.0))
    assert (x2hat[1, 0] == complex(0.0, 0.0))
    assert (x2hat[1, 1] == complex(1.0, 0.0))
    assert (x2hat[2, 0] == complex(0.0, 0.0))
    assert (x2hat[2, 1] == complex(0.0, 0.0))
    assert (x2hat[2, 2] == complex(1.0, 0.0))

    x2.unserialize([4.0, 5.0, 6.0, 7.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == complex(4.0, 5.0))
    assert (x2hat[1, 2] == complex(6.0, 7.0))

    assert (x2hat[0, 0] == complex(1.0, 0.0))
    assert (x2hat[0, 1] == complex(0.0, 0.0))
    assert (x2hat[1, 0] == complex(0.0, 0.0))
    assert (x2hat[1, 1] == complex(1.0, 0.0))
    assert (x2hat[2, 0] == complex(0.0, 0.0))
    assert (x2hat[2, 1] == complex(0.0, 0.0))
    assert (x2hat[2, 2] == complex(1.0, 0.0))

    x2.unserialize([8.0, 9.0, 10.0, 11.0, 12.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == complex(8.0, 9.0))
    assert (x2hat[1, 2] == complex(10.0, 11.0))

    assert (x2hat[0, 0] == complex(1.0, 0.0))
    assert (x2hat[0, 1] == complex(0.0, 0.0))
    assert (x2hat[1, 0] == complex(0.0, 0.0))
    assert (x2hat[1, 1] == complex(1.0, 0.0))
    assert (x2hat[2, 0] == complex(0.0, 0.0))
    assert (x2hat[2, 1] == complex(0.0, 0.0))
    assert (x2hat[2, 2] == complex(1.0, 0.0))

def test_CN_operator_parenthesis():
    from lielab.domain import CN

    xzero = CN.from_shape(0)
    xzero.unserialize([])

    # Out of bounds
    assert (np.isnan(xzero(-1)))
    assert (np.isnan(xzero(0)))
    assert (np.isnan(xzero(1)))

    # Out of bounds
    assert (np.isnan(np.real(xzero(0, -1))))
    assert (np.isnan(np.imag(xzero(0, -1))))
    assert (np.isnan(np.real(xzero(-1, 0))))
    assert (np.isnan(np.imag(xzero(-1, 0))))
    assert (np.isnan(np.real(xzero(-1, -1))))
    assert (np.isnan(np.imag(xzero(-1, -1))))
    assert (np.isnan(np.real(xzero(0, 0))))
    assert (np.isnan(np.imag(xzero(0, 0))))
    assert (np.isnan(np.real(xzero(0, 1))))
    assert (np.isnan(np.imag(xzero(0, 1))))
    assert (np.isnan(np.real(xzero(1, 0))))
    assert (np.isnan(np.imag(xzero(1, 0))))
    assert (np.isnan(np.real(xzero(1, 1))))
    assert (np.isnan(np.imag(xzero(1, 1))))

    x1 = CN(1)
    x1.unserialize([1.0, 2.0])

    # In bounds
    assert (x1(0) == 1.0)
    assert (x1(1) == 2.0)
    assert (x1(-1) == 2.0)
    assert (x1(-2) == 1.0)

    # Out of bounds
    assert (np.isnan(x1(-3)))
    assert (np.isnan(x1(2)))

    # In bounds
    assert (x1(0, 0) == complex(1.0, 0.0))
    assert (x1(0, 1) == complex(1.0, 2.0))
    assert (x1(1, 0) == complex(0.0, 0.0))
    assert (x1(1, 1) == complex(1.0, 0.0))
    assert (x1(-1, -1) == complex(1.0, 0.0))
    assert (x1(-1, -2) == complex(0.0, 0.0))
    assert (x1(-2, -1) == complex(1.0, 2.0))
    assert (x1(-2, -2) == complex(1.0, 0.0))

    # Out of bounds
    assert (np.isnan(np.real(x1(0, -3))))
    assert (np.isnan(np.imag(x1(0, -3))))
    assert (np.isnan(np.real(x1(-3, 0))))
    assert (np.isnan(np.imag(x1(-3, 0))))
    assert (np.isnan(np.real(x1(-3, -3))))
    assert (np.isnan(np.imag(x1(-3, -3))))
    assert (np.isnan(np.real(x1(0, 2))))
    assert (np.isnan(np.imag(x1(0, 2))))
    assert (np.isnan(np.real(x1(2, 0))))
    assert (np.isnan(np.imag(x1(2, 0))))
    assert (np.isnan(np.real(x1(2, 2))))
    assert (np.isnan(np.imag(x1(2, 2))))

    x2 = CN(2)
    x2.unserialize([1.0, 2.0, 3.0, 4.0])

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
    assert (x2(0, 0) == complex(1.0, 0.0))
    assert (x2(0, 1) == complex(0.0, 0.0))
    assert (x2(0, 2) == complex(1.0, 2.0))
    assert (x2(1, 0) == complex(0.0, 0.0))
    assert (x2(1, 1) == complex(1.0, 0.0))
    assert (x2(1, 2) == complex(3.0, 4.0))
    assert (x2(2, 0) == complex(0.0, 0.0))
    assert (x2(2, 1) == complex(0.0, 0.0))
    assert (x2(2, 2) == complex(1.0, 0.0))
    assert (x2(-1, -1) == complex(1.0, 0.0))
    assert (x2(-1, -2) == complex(0.0, 0.0))
    assert (x2(-1, -3) == complex(0.0, 0.0))
    assert (x2(-2, -1) == complex(3.0, 4.0))
    assert (x2(-2, -2) == complex(1.0, 0.0))
    assert (x2(-2, -3) == complex(0.0, 0.0))
    assert (x2(-3, -1) == complex(1.0, 2.0))
    assert (x2(-3, -2) == complex(0.0, 0.0))
    assert (x2(-3, -3) == complex(1.0, 0.0))

    # Out of bounds
    assert (np.isnan(np.real(x2(0, -4))))
    assert (np.isnan(np.imag(x2(0, -4))))
    assert (np.isnan(np.real(x2(-4, 0))))
    assert (np.isnan(np.imag(x2(-4, 0))))
    assert (np.isnan(np.real(x2(-4, -4))))
    assert (np.isnan(np.imag(x2(-4, -4))))
    assert (np.isnan(np.real(x2(0, 3))))
    assert (np.isnan(np.imag(x2(0, 3))))
    assert (np.isnan(np.real(x2(3, 0))))
    assert (np.isnan(np.imag(x2(3, 0))))
    assert (np.isnan(np.real(x2(3, 3))))
    assert (np.isnan(np.imag(x2(3, 3))))

def test_CN_operator_bracket():
    from lielab.domain import CN

    xzero = CN.from_shape(0)
    xzero.unserialize([])

    # Out of bounds
    assert (np.isnan(np.real(xzero[-1])))
    assert (np.isnan(np.imag(xzero[-1])))
    assert (np.isnan(np.real(xzero[0])))
    assert (np.isnan(np.imag(xzero[0])))
    assert (np.isnan(np.real(xzero[1])))
    assert (np.isnan(np.imag(xzero[1])))

    x1 = CN(1)
    x1.unserialize([1.0, 2.0])

    # In bounds
    assert (x1[0] == complex(1.0, 2.0))
    assert (x1[-1] == complex(1.0, 2.0))

    # Out of bounds
    assert (np.isnan(np.real(x1[-2])))
    assert (np.isnan(np.imag(x1[-2])))
    assert (np.isnan(np.real(x1[1])))
    assert (np.isnan(np.imag(x1[1])))

    x2 = CN(2)
    x2.unserialize([1.0, 2.0, 3.0, 4.0])

    # In bounds
    assert (x2[0] == complex(1.0, 2.0))
    assert (x2[1] == complex(3.0, 4.0))
    assert (x2[-1] == complex(3.0, 4.0))
    assert (x2[-2] == complex(1.0, 2.0))

    # Out of bounds
    assert (np.isnan(np.real(x2[-3])))
    assert (np.isnan(np.imag(x2[-3])))
    assert (np.isnan(np.real(x2[2])))
    assert (np.isnan(np.imag(x2[2])))

def test_CN_math_ops_CN():
    from lielab.domain import CN

    x1 = CN(2)
    x2 = CN(2)
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

def test_CN_from_vector():
    """
    * Tests the from_vector operation.
    """

    from lielab.domain import CN

    x0 = CN.from_vector([])
    x0bar = x0.serialize()

    assert (x0.get_shape() == 1)
    assert (x0bar.size == 0)

    x1 = CN.from_vector([1.0])
    x1bar = x1.serialize()

    assert (x1.get_shape() == 2)
    assert (x1bar.size == 2)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 0.0)

    x2 = CN.from_vector([1.0, 2.0])
    x2bar = x2.serialize()

    assert (x2.get_shape() == 2)
    assert (x2bar.size == 2)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)

    x3 = CN.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.serialize()

    assert (x3.get_shape() == 3)
    assert (x3bar.size == 4)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)
    assert (x3bar[3] == 0.0)

def test_CN_from_to_complex_vector():
    """
    * Tests the from_to_complex_vector operations.
    """

    from lielab.domain import CN

    x0 = CN.from_complex_vector([])
    x0bar = x0.to_complex_vector()

    assert (x0.get_shape() == 1)
    assert (x0bar.size == 0)

    x1 = CN.from_complex_vector([1.0])
    x1bar = x1.to_complex_vector()

    assert (x1.get_shape() == 2)
    assert (x1bar.size == 1)
    assert (x1bar[0] == complex(1.0, 0.0))

    x2 = CN.from_complex_vector([1.0 + 2.0j])
    x2bar = x2.to_complex_vector()

    assert (x2.get_shape() == 2)
    assert (x2bar.size == 1)
    assert (x2bar[0] == complex(1.0, 2.0))

    x3 = CN.from_complex_vector([1.0 + 2.0j, 3.0])
    x3bar = x3.to_complex_vector()

    assert (x3.get_shape() == 3)
    assert (x3bar.size == 2)
    assert (x3bar[0] == complex(1.0, 2.0))
    assert (x3bar[1] == complex(3.0, 0.0))

def test_CN_project():
    from lielab.domain import CN

    rand_2_2 = np.random.rand(2, 2) + 1j*np.random.rand(2, 2)
    proj_2_2 = CN.project(rand_2_2)

    assert (proj_2_2.shape[0] == 2)
    assert (proj_2_2.shape[1] == 2)
    assert (proj_2_2[0, 0] == complex(1.0, 0.0))
    assert (proj_2_2[0, 1] == rand_2_2[0, 1])
    assert (proj_2_2[1, 0] == complex(0.0, 0.0))
    assert (proj_2_2[1, 1] == complex(1.0, 0.0))

    rand_3_3 = np.random.rand(3, 3) + 1j*np.random.rand(3, 3)
    proj_3_3 = CN.project(rand_3_3)

    assert (proj_3_3.shape[0] == 3)
    assert (proj_3_3.shape[1] == 3)
    assert (proj_3_3[0, 0] == complex(1.0, 0.0))
    assert (proj_3_3[0, 1] == complex(0.0, 0.0))
    assert (proj_3_3[0, 2] == rand_3_3[0, 2])
    assert (proj_3_3[1, 0] == complex(0.0, 0.0))
    assert (proj_3_3[1, 1] == complex(1.0, 0.0))
    assert (proj_3_3[1, 2] == rand_3_3[1, 2])
    assert (proj_3_3[2, 0] == complex(0.0, 0.0))
    assert (proj_3_3[2, 1] == complex(0.0, 0.0))
    assert (proj_3_3[2, 2] == complex(1.0, 0.0))

    rand_2_3 = np.random.rand(2, 3) + 1j*np.random.rand(2, 3)
    proj_2_3 = CN.project(rand_2_3)

    assert (proj_2_3.shape[0] == 2)
    assert (proj_2_3.shape[1] == 2)
    assert (proj_2_3[0, 0] == complex(1.0, 0.0))
    assert (proj_2_3[0, 1] == rand_2_3[0, 1])
    assert (proj_2_3[1, 0] == complex(0.0, 0.0))
    assert (proj_2_3[1, 1] == complex(1.0, 0.0))

    rand_3_2 = np.random.rand(3, 2) + np.random.rand(3, 2)
    proj_3_2 = CN.project(rand_3_2)

    assert (proj_3_2.shape[0] == 2)
    assert (proj_3_2.shape[1] == 2)
    assert (proj_3_2[0, 0] == complex(1.0, 0.0))
    assert (proj_3_2[0, 1] == rand_3_2[0, 1])
    assert (proj_3_2[1, 0] == complex(0.0, 0.0))
    assert (proj_3_2[1, 1] == complex(1.0, 0.0))
