from lielab.testing import *
import pytest

def complex(a,b):
    return a + b*1j

def test_glc_to_string():
    from lielab.domain import glc

    x0 = glc(0)
    assert (x0.to_string() == "gl(0, C)")
    x1 = glc(1)
    assert (x1.to_string() == "gl(1, C)")
    x10 = glc(10)
    assert (x10.to_string() == "gl(10, C)")


def test_glc_main_initializer():
    from lielab.domain import glc

    xblank = glc()
    assert (xblank.get_dimension() == 0)

    x0 = glc(0)
    assert (x0.get_dimension() == 0)
    x1 = glc(1)
    assert (x1.get_dimension() == 2)
    x10 = glc(10)
    assert (x10.get_dimension() == 200)

def test_glc_matrix_initializer():
    from lielab.domain import glc

    x0 = glc(np.random.rand(0, 0) + 1j*np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = glc(np.random.rand(1, 1) + 1j*np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = glc(np.random.rand(2, 2) + 1j*np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        glc(np.random.rand(2, 3) + 1j*np.random.rand(2, 3))
    
    with pytest.raises(RuntimeError):
        glc(np.random.rand(3, 2) + 1j*np.random.rand(3, 2))

def test_glc_basis_initializer():
    from lielab.domain import glc

    xm10 = glc.basis(-1, 0)
    assert (xm10.get_dimension() == 0)
    xm10bar = xm10.get_vector()
    assert (xm10bar.size == 0)

    x00 = glc.basis(0, 0)
    assert (x00.get_dimension() == 0)
    x00bar = x00.get_vector()
    assert (x00bar.size == 0)

    x10 = glc.basis(1, 0)
    assert (x10.get_dimension() == 0)
    x10bar = x10.get_vector()
    assert (x10bar.size == 0)

    x01 = glc.basis(0, 1)
    assert (x01.get_dimension() == 2)
    x01bar = x01.get_vector()
    assert (x01bar.size == 2)
    assert (x01bar[0] == 1.0)
    assert (x01bar[1] == 0.0)

    x11 = glc.basis(1, 1)
    assert (x11.get_dimension() == 2)
    x11bar = x11.get_vector()
    assert (x11bar.size == 2)
    assert (x11bar[0] == 0.0)
    assert (x11bar[1] == 1.0)

    x21 = glc.basis(2, 1)
    assert (x21.get_dimension() == 2)
    x21bar = x21.get_vector()
    assert (x21bar.size == 2)
    assert (x21bar[0] == 0.0)
    assert (x21bar[1] == 0.0)

    x02 = glc.basis(0, 2)
    assert (x02.get_dimension() == 8)
    x02bar = x02.get_vector()
    assert (x02bar.size == 8)
    assert (x02bar[0] == 1.0)
    assert (x02bar[1] == 0.0)
    assert (x02bar[2] == 0.0)
    assert (x02bar[3] == 0.0)
    assert (x02bar[4] == 0.0)
    assert (x02bar[5] == 0.0)
    assert (x02bar[6] == 0.0)
    assert (x02bar[7] == 0.0)

def test_glc_from_shape_initializer():
    from lielab.domain import glc

    x0 = glc.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = glc.from_shape(1)
    assert (x1.get_dimension() == 2)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = glc.from_shape(2)
    assert (x2.get_dimension() == 8)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_glc_get_dimension():
    from lielab.domain import glc

    zero = glc(0)
    one = glc(1)
    two = glc(2)
    three = glc(3)
    four = glc(4)
    five = glc(5)
    six = glc(6)
    seven = glc(7)
    eight = glc(8)

    # Dimensions
    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 2)
    assert (two.get_dimension() == 8)
    assert (three.get_dimension() == 18)
    assert (four.get_dimension() == 32)
    assert (five.get_dimension() == 50)
    assert (six.get_dimension() == 72)
    assert (seven.get_dimension() == 98)
    assert (eight.get_dimension() == 128)

def test_glc_set_get_vector():
    """
    * Tests the set/get_vector operation.
    """

    from lielab.domain import glc

    x0 = glc(0)
    x0.set_vector([])
    x0bar = x0.get_vector()

    assert (x0bar.size == 0)

    x1 = glc(1)
    x1.set_vector([1.0, 2.0, 3.0])
    x1bar = x1.get_vector()

    assert (x1bar.size == 2)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 2.0)

    x1.set_vector([4.0, 5.0])
    x1bar = x1.get_vector()

    assert (x1bar.size == 2)
    assert (x1bar[0] == 4.0)
    assert (x1bar[1] == 5.0)

    x1.set_vector([6.0])
    x1bar = x1.get_vector()

    assert (x1bar.size == 2)
    assert (x1bar[0] == 6.0)
    assert (x1bar[1] == 5.0)

    x2 = glc(2)
    x2.set_vector([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 8)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 3.0)
    assert (x2bar[3] == 4.0)
    assert (x2bar[4] == 5.0)
    assert (x2bar[5] == 6.0)
    assert (x2bar[6] == 7.0)
    assert (x2bar[7] == 0.0)

    x2.set_vector([8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 8)
    assert (x2bar[0] == 8.0)
    assert (x2bar[1] == 9.0)
    assert (x2bar[2] == 10.0)
    assert (x2bar[3] == 11.0)
    assert (x2bar[4] == 12.0)
    assert (x2bar[5] == 13.0)
    assert (x2bar[6] == 14.0)
    assert (x2bar[7] == 15.0)

    x2.set_vector([16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 8)
    assert (x2bar[0] == 16.0)
    assert (x2bar[1] == 17.0)
    assert (x2bar[2] == 18.0)
    assert (x2bar[3] == 19.0)
    assert (x2bar[4] == 20.0)
    assert (x2bar[5] == 21.0)
    assert (x2bar[6] == 22.0)
    assert (x2bar[7] == 23.0)

def test_glc_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import glc

    x0 = glc(0)
    x0.set_vector([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = glc(1)
    x1.set_vector([1.0, 2.0, 3.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == complex(1.0, 2.0))

    x1.set_vector([4.0, 5.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == complex(4.0, 5.0))

    x1.set_vector([6.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == complex(6.0, 5.0))

    x2 = glc(2)
    x2.set_vector([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == complex(1.0, 2.0))
    assert (x2hat[0, 1] == complex(3.0, 4.0))
    assert (x2hat[1, 0] == complex(5.0, 6.0))
    assert (x2hat[1, 1] == complex(7.0, 0.0))

    x2.set_vector([8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == complex(8.0, 9.0))
    assert (x2hat[0, 1] == complex(10.0, 11.0))
    assert (x2hat[1, 0] == complex(12.0, 13.0))
    assert (x2hat[1, 1] == complex(14.0, 15.0))

    x2.set_vector([16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == complex(16.0, 17.0))
    assert (x2hat[0, 1] == complex(18.0, 19.0))
    assert (x2hat[1, 0] == complex(20.0, 21.0))
    assert (x2hat[1, 1] == complex(22.0, 23.0))

def test_glc_operator_parenthesis():
    from lielab.domain import glc

    x0 = glc.from_shape(0)
    x0.set_vector([])

    # Out of bounds
    assert (np.isnan(x0(-1)))
    assert (np.isnan(x0(0)))
    assert (np.isnan(x0(1)))

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

    x1 = glc(1)
    x1.set_vector([1.0, 2.0])

    # In bounds
    assert (x1(0) == 1.0)
    assert (x1(1) == 2.0)
    assert (x1(-1) == 2.0)
    assert (x1(-2) == 1.0)

    # Out of bounds
    assert (np.isnan(x1(-3)))
    assert (np.isnan(x1(2)))

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

    x2 = glc(2)
    x2.set_vector([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])

    # In bounds
    assert (x2(0) == 1.0)
    assert (x2(1) == 2.0)
    assert (x2(2) == 3.0)
    assert (x2(3) == 4.0)
    assert (x2(4) == 5.0)
    assert (x2(5) == 6.0)
    assert (x2(6) == 7.0)
    assert (x2(7) == 8.0)
    assert (x2(-1) == 8.0)
    assert (x2(-2) == 7.0)
    assert (x2(-3) == 6.0)
    assert (x2(-4) == 5.0)
    assert (x2(-5) == 4.0)
    assert (x2(-6) == 3.0)
    assert (x2(-7) == 2.0)
    assert (x2(-8) == 1.0)

    # Out of bounds
    assert (np.isnan(x2(-9)))
    assert (np.isnan(x2(8)))

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

# TODO: mathops int here

def test_glc_math_ops_double():
    from lielab.domain import glc

    x1 = glc(2)
    x1.set_vector([1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0])

    x1_lm_2 = 2.0*x1
    assert (x1_lm_2(0, 0) == complex(2.5, 5.0))
    assert (x1_lm_2(0, 1) == complex(7.5, 2.0))
    assert (x1_lm_2(1, 0) == complex(2.0, 2.0))
    assert (x1_lm_2(1, 1) == complex(2.0, 2.0))

    x1_rm_2 = x1*2.0
    assert (x1_rm_2(0, 0) == complex(2.5, 5.0))
    assert (x1_rm_2(0, 1) == complex(7.5, 2.0))
    assert (x1_rm_2(1, 0) == complex(2.0, 2.0))
    assert (x1_rm_2(1, 1) == complex(2.0, 2.0))

    x1 *= 2.0
    assert (x1(0, 0) == complex(2.5, 5.0))
    assert (x1(0, 1) == complex(7.5, 2.0))
    assert (x1(1, 0) == complex(2.0, 2.0))
    assert (x1(1, 1) == complex(2.0, 2.0))

    x1.set_vector([1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0])

    x1_lm_2j = (2.0*1.0j)*x1
    assert (x1_lm_2j(0, 0) == complex(-5.0, 2.5))
    assert (x1_lm_2j(0, 1) == complex(-2.0, 7.5))
    assert (x1_lm_2j(1, 0) == complex(-2.0, 2.0))
    assert (x1_lm_2j(1, 1) == complex(-2.0, 2.0))

    x1_rm_2j = x1*(2.0*1.0j)
    assert (x1_rm_2j(0, 0) == complex(-5.0, 2.5))
    assert (x1_rm_2j(0, 1) == complex(-2.0, 7.5))
    assert (x1_rm_2j(1, 0) == complex(-2.0, 2.0))
    assert (x1_rm_2j(1, 1) == complex(-2.0, 2.0))

    x1 *= 2.0*1.0j
    assert (x1(0, 0) == complex(-5.0, 2.5))
    assert (x1(0, 1) == complex(-2.0, 7.5))
    assert (x1(1, 0) == complex(-2.0, 2.0))
    assert (x1(1, 1) == complex(-2.0, 2.0))

    x1.set_vector([1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0])

    x1_d_2 = x1/2.0
    assert (x1_d_2(0, 0) == complex(0.625, 1.25))
    assert (x1_d_2(0, 1) == complex(1.875, 0.5))
    assert (x1_d_2(1, 0) == complex(0.5, 0.5))
    assert (x1_d_2(1, 1) == complex(0.5, 0.5))

    x1 /= 2.0
    assert (x1(0, 0) == complex(0.625, 1.25))
    assert (x1(0, 1) == complex(1.875, 0.5))
    assert (x1(1, 0) == complex(0.5, 0.5))
    assert (x1(1, 1) == complex(0.5, 0.5))

    x1.set_vector([1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0])

    x1_d_2j = x1/(2.0*1.0j)
    assert (x1_d_2j(0, 0) == complex(1.25, -0.625))
    assert (x1_d_2j(0, 1) == complex(0.5, -1.875))
    assert (x1_d_2j(1, 0) == complex(0.5, -0.5))
    assert (x1_d_2j(1, 1) == complex(0.5, -0.5))

    x1 /= 2.0*1.0j
    assert (x1(0, 0) == complex(1.25, -0.625))
    assert (x1(0, 1) == complex(0.5, -1.875))
    assert (x1(1, 0) == complex(0.5, -0.5))
    assert (x1(1, 1) == complex(0.5, -0.5))

def test_glc_math_ops_cn():
    from lielab.domain import glc

    x1 = glc(2)
    x2 = glc(2)
    x1.set_vector([1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0])
    x2.set_vector([1.25, 2.5, 3.75, 1.0, 1.0, 1.0, 1.0, 1.0])

    x1_add_x2 = x1 + x2
    assert (x1_add_x2(0, 0) == complex(2.25, 4.5))
    assert (x1_add_x2(0, 1) == complex(6.75, 5.0))
    assert (x1_add_x2(1, 0) == complex(2.0, 2.0))
    assert (x1_add_x2(1, 1) == complex(2.0, 2.0))

    x1 += x2
    assert (x1(0, 0) == complex(2.25, 4.5))
    assert (x1(0, 1) == complex(6.75, 5.0))
    assert (x1(1, 0) == complex(2.0, 2.0))
    assert (x1(1, 1) == complex(2.0, 2.0))

    x1.set_vector([1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0])

    x1_sub_x2 = x1 - x2
    assert (x1_sub_x2(0, 0) == complex(-0.25, -0.5))
    assert (x1_sub_x2(0, 1) == complex(-0.75, 3.0))
    assert (x1_sub_x2(1, 0) == complex(0.0, 0.0))
    assert (x1_sub_x2(1, 1) == complex(0.0, 0.0))

    x1 -= x2
    assert (x1(0, 0) == complex(-0.25, -0.5))
    assert (x1(0, 1) == complex(-0.75, 3.0))
    assert (x1(1, 0) == complex(0.0, 0.0))
    assert (x1(1, 1) == complex(0.0, 0.0))

    x1.set_vector([1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 1.0])

    x1_unary_sub = (-x1)
    assert (x1_unary_sub(0, 0) == complex(-1.0, -2.0))
    assert (x1_unary_sub(0, 1) == complex(-3.0, -4.0))
    assert (x1_unary_sub(1, 0) == complex(-1.0, -1.0))
    assert (x1_unary_sub(1, 1) == complex(-1.0, -1.0))

def test_glc_from_vector():
    """
    * Tests the from_vector operation.
    """

    from lielab.domain import glc

    x0 = glc.from_vector([])
    x0bar = x0.get_vector()

    assert (x0.get_shape() == 0)
    assert (x0bar.size == 0)

    x1 = glc.from_vector([1.0])
    x1bar = x1.get_vector()

    assert (x1.get_shape() == 1)
    assert (x1bar.size == 2)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 0.0)

    x2 = glc.from_vector([1.0, 2.0])
    x2bar = x2.get_vector()

    assert (x2.get_shape() == 1)
    assert (x2bar.size == 2)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)

    x3 = glc.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.get_vector()

    assert (x3.get_shape() == 2)
    assert (x3bar.size == 8)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)
    assert (x3bar[3] == 0.0)
    assert (x3bar[4] == 0.0)
    assert (x3bar[5] == 0.0)
    assert (x3bar[6] == 0.0)
    assert (x3bar[7] == 0.0)

def test_glc_from_complex_vector():
    """
    * Tests the from_complex_vector operation.
    """

    from lielab.domain import glc

    x0 = glc.from_complex_vector([])
    x0bar = x0.get_vector()

    assert (x0.get_shape() == 0)
    assert (x0bar.size == 0)

    x1 = glc.from_complex_vector([1.0])
    x1bar = x1.get_vector()

    assert (x1.get_shape() == 1)
    assert (x1bar.size == 2)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 0.0)

    x2 = glc.from_complex_vector([1.0 + 2.0*1.0j])
    x2bar = x2.get_vector()

    assert (x2.get_shape() == 1)
    assert (x2bar.size == 2)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)

    x3 = glc.from_complex_vector([1.0 + 2.0*1.0j, 3.0])
    x3bar = x3.get_vector()

    assert (x3.get_shape() == 2)
    assert (x3bar.size == 8)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)
    assert (x3bar[3] == 0.0)
    assert (x3bar[4] == 0.0)
    assert (x3bar[5] == 0.0)
    assert (x3bar[6] == 0.0)
    assert (x3bar[7] == 0.0)

def test_glc_project():
    from lielab.domain import glc

    rand_2_2 = np.random.rand(2, 2) + 1.0j*np.random.rand(2, 2)
    proj_2_2 = glc.project(rand_2_2)

    assert (proj_2_2.shape[0] == 2)
    assert (proj_2_2.shape[1] == 2)
    assert (proj_2_2[0, 0] == rand_2_2[0, 0])
    assert (proj_2_2[0, 1] == rand_2_2[0, 1])
    assert (proj_2_2[1, 0] == rand_2_2[1, 0])
    assert (proj_2_2[1, 1] == rand_2_2[1, 1])

    rand_3_3 = np.random.rand(3, 3) + 1.0j*np.random.rand(3, 3)
    proj_3_3 = glc.project(rand_3_3)

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

    rand_2_3 = np.random.rand(2, 3) + 1.0j*np.random.rand(2, 3)
    proj_2_3 = glc.project(rand_2_3)

    assert (proj_2_3.shape[0] == 2)
    assert (proj_2_3.shape[1] == 2)
    assert (proj_2_3[0, 0] == rand_2_3[0, 0])
    assert (proj_2_3[0, 1] == rand_2_3[0, 1])
    assert (proj_2_3[1, 0] == rand_2_3[1, 0])
    assert (proj_2_3[1, 1] == rand_2_3[1, 1])

    rand_3_2 = np.random.rand(3, 2) + 1.0j*np.random.rand(3, 2)
    proj_3_2 = glc.project(rand_3_2)

    assert (proj_3_2.shape[0] == 2)
    assert (proj_3_2.shape[1] == 2)
    assert (proj_3_2[0, 0] == rand_3_2[0, 0])
    assert (proj_3_2[0, 1] == rand_3_2[0, 1])
    assert (proj_3_2[1, 0] == rand_3_2[1, 0])
    assert (proj_3_2[1, 1] == rand_3_2[1, 1])

def test_glc_get_from_vector():
    """
    * Tests the get/from vector operation.
    """

    from lielab.domain import glc

    x0 = glc.from_vector([])
    x0bar = x0.get_vector()

    assert (x0bar.size == 0)
    assert (x0.get_dimension() == 0)
    
    x1 = glc.from_vector([1.0])
    x1bar = x1.get_vector()

    assert (x1bar.size == 2)
    assert (x1.get_dimension() == 2)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 0.0)

    x2 = glc.from_vector([1.0, 2.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 2)
    assert (x2.get_dimension() == 2)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)

    x3 = glc.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 8)
    assert (x3.get_dimension() == 8)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)
    assert (x3bar[3] == 0.0)
    assert (x3bar[4] == 0.0)
    assert (x3bar[5] == 0.0)
    assert (x3bar[6] == 0.0)
    assert (x3bar[7] == 0.0)

    x4 = glc.from_vector([1.0, 2.0, 3.0, 4.0])
    x4bar = x4.get_vector()

    assert (x4bar.size == 8)
    assert (x4.get_dimension() == 8)
    assert (x4bar[0] == 1.0)
    assert (x4bar[1] == 2.0)
    assert (x4bar[2] == 3.0)
    assert (x4bar[3] == 4.0)
    assert (x4bar[4] == 0.0)
    assert (x4bar[5] == 0.0)
    assert (x4bar[6] == 0.0)
    assert (x4bar[7] == 0.0)

    x5 = glc.from_vector([1.0, 2.0, 3.0, 4.0, 5.0])
    x5bar = x5.get_vector()

    assert (x5bar.size == 8)
    assert (x5.get_dimension() == 8)
    assert (x5bar[0] == 1.0)
    assert (x5bar[1] == 2.0)
    assert (x5bar[2] == 3.0)
    assert (x5bar[3] == 4.0)
    assert (x5bar[4] == 5.0)
    assert (x5bar[5] == 0.0)
    assert (x5bar[6] == 0.0)
    assert (x5bar[7] == 0.0)
