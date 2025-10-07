from lielab.testing import *
import pytest

def complex(a,b):
    return a + b*1j

def test_cn_to_string():
    from lielab.domain import cn

    xzero = cn.from_shape(0)
    assert (xzero.to_string() == "c^nan")
    x0 = cn(0)
    assert (x0.to_string() == "c^0")
    x1 = cn(1)
    assert (x1.to_string() == "c^1")
    x10 = cn(10)
    assert (x10.to_string() == "c^10")

def test_cn_main_initializer():
    from lielab.domain import cn

    xblank = cn()
    assert (xblank.get_dimension() == 0)

    x0 = cn(0)
    assert (x0.get_dimension() == 0)
    x1 = cn(1)
    assert (x1.get_dimension() == 2)
    x10 = cn(10)
    assert (x10.get_dimension() == 20)

def test_cn_matrix_initializer():
    from lielab.domain import cn

    x0 = cn(np.random.rand(0, 0) + 1j*np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = cn(np.random.rand(1, 1) + 1j*np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = cn(np.random.rand(2, 2) + 1j*np.random.rand(1, 1))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        cn(np.random.rand(2, 3) + 1j*np.random.rand(2, 3))

    with pytest.raises(RuntimeError):
        cn(np.random.rand(3, 2) + 1j*np.random.rand(3, 2))

def test_cn_basis_initializer():

    from lielab.domain import cn

    xm10 = cn.basis(-1, 0)
    assert (xm10.get_dimension() == 0)
    xm10bar = xm10.get_vector()
    assert (xm10bar.size == 0)

    x00 = cn.basis(0, 0)
    assert (x00.get_dimension() == 0)
    x00bar = x00.get_vector()
    assert (x00bar.size == 0)

    x10 = cn.basis(1, 0)
    assert (x10.get_dimension() == 0)
    x10bar = x10.get_vector()
    assert (x10bar.size == 0)

    x01 = cn.basis(0, 1)
    assert (x01.get_dimension() == 2)
    x01bar = x01.get_vector()
    assert (x01bar.size == 2)
    assert (x01bar[0] == 1.0)
    assert (x01bar[1] == 0.0)

    x11 = cn.basis(1, 1)
    assert (x11.get_dimension() == 2)
    x11bar = x11.get_vector()
    assert (x11bar.size == 2)
    assert (x11bar[0] == 0.0)
    assert (x11bar[1] == 1.0)

    x21 = cn.basis(2, 1)
    assert (x21.get_dimension() == 2)
    x21bar = x21.get_vector()
    assert (x21bar.size == 2)
    assert (x21bar[0] == 0.0)
    assert (x21bar[1] == 0.0)

    x02 = cn.basis(0, 2)
    assert (x02.get_dimension() == 4)
    x02bar = x02.get_vector()
    assert (x02bar.size == 4)
    assert (x02bar[0] == 1.0)
    assert (x02bar[1] == 0.0)
    assert (x02bar[2] == 0.0)
    assert (x02bar[3] == 0.0)


def test_cn_from_shape_initializer():

    from lielab.domain import cn

    x0 = cn.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = cn.from_shape(1)
    assert (x1.get_dimension() == 0)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = cn.from_shape(2)
    assert (x2.get_dimension() == 2)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)


def test_cn_get_dimension():

    from lielab.domain import cn

    veryzero = cn.from_shape(0)
    zero = cn(0)
    one = cn(1)
    two = cn(2)
    three = cn(3)
    four = cn(4)
    five = cn(5)
    six = cn(6)
    seven = cn(7)
    eight = cn(8)

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


def test_cn_set_get_vector():
    """
    Tests the set/get_vector operation.
    """

    from lielab.domain import cn

    xzero = cn.from_shape(0)
    xzero.set_vector([])
    xzerobar = xzero.get_vector()

    assert (xzerobar.size == 0)

    x0 = cn(0)
    x0.set_vector([])
    x0bar = x0.get_vector()

    assert (x0bar.size == 0)

    x1 = cn(1)
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

    x2 = cn(2)
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


def test_cn_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import cn

    xzero = cn.from_shape(0)
    xzero.set_vector([])
    xzerohat = xzero.get_matrix()

    assert (xzerohat.shape[0] == 0)
    assert (xzerohat.shape[1] == 0)

    x0 = cn(0)
    x0.set_vector([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 1)
    assert (x0hat.shape[1] == 1)
    assert (x0hat[0, 0] == complex(0.0, 0.0))

    x1 = cn(1)
    x1.set_vector([1.0, 2.0, 3.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == complex(1.0, 2.0))

    assert (x1hat[0, 0] == complex(0.0, 0.0))
    assert (x1hat[1, 0] == complex(0.0, 0.0))
    assert (x1hat[1, 1] == complex(0.0, 0.0))

    x1.set_vector([4.0, 5.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == complex(4.0, 5.0))

    assert (x1hat[0, 0] == complex(0.0, 0.0))
    assert (x1hat[1, 0] == complex(0.0, 0.0))
    assert (x1hat[1, 1] == complex(0.0, 0.0))

    x1.set_vector([6.0])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 2)
    assert (x1hat.shape[1] == 2)
    assert (x1hat[0, 1] == complex(6.0, 5.0))

    assert (x1hat[0, 0] == complex(0.0, 0.0))
    assert (x1hat[1, 0] == complex(0.0, 0.0))
    assert (x1hat[1, 1] == complex(0.0, 0.0))

    x2 = cn(2)
    x2.set_vector([1.0, 2.0, 3.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == complex(1.0, 2.0))
    assert (x2hat[1, 2] == complex(3.0, 0.0))

    assert (x2hat[0, 0] == complex(0.0, 0.0))
    assert (x2hat[0, 1] == complex(0.0, 0.0))
    assert (x2hat[1, 0] == complex(0.0, 0.0))
    assert (x2hat[1, 1] == complex(0.0, 0.0))
    assert (x2hat[2, 0] == complex(0.0, 0.0))
    assert (x2hat[2, 1] == complex(0.0, 0.0))
    assert (x2hat[2, 2] == complex(0.0, 0.0))

    x2.set_vector([4.0, 5.0, 6.0, 7.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == complex(4.0, 5.0))
    assert (x2hat[1, 2] == complex(6.0, 7.0))

    assert (x2hat[0, 0] == complex(0.0, 0.0))
    assert (x2hat[0, 1] == complex(0.0, 0.0))
    assert (x2hat[1, 0] == complex(0.0, 0.0))
    assert (x2hat[1, 1] == complex(0.0, 0.0))
    assert (x2hat[2, 0] == complex(0.0, 0.0))
    assert (x2hat[2, 1] == complex(0.0, 0.0))
    assert (x2hat[2, 2] == complex(0.0, 0.0))

    x2.set_vector([8.0, 9.0, 10.0, 11.0, 12.0])
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 3)
    assert (x2hat.shape[1] == 3)
    assert (x2hat[0, 2] == complex(8.0, 9.0))
    assert (x2hat[1, 2] == complex(10.0, 11.0))

    assert (x2hat[0, 0] == complex(0.0, 0.0))
    assert (x2hat[0, 1] == complex(0.0, 0.0))
    assert (x2hat[1, 0] == complex(0.0, 0.0))
    assert (x2hat[1, 1] == complex(0.0, 0.0))
    assert (x2hat[2, 0] == complex(0.0, 0.0))
    assert (x2hat[2, 1] == complex(0.0, 0.0))
    assert (x2hat[2, 2] == complex(0.0, 0.0))


def test_cn_operator_parenthesis():
    from lielab.domain import cn

    xzero = cn.from_shape(0)
    xzero.set_vector([])

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

    x1 = cn(1)
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
    assert (x1(0, 0) == complex(0.0, 0.0))
    assert (x1(0, 1) == complex(1.0, 2.0))
    assert (x1(1, 0) == complex(0.0, 0.0))
    assert (x1(1, 1) == complex(0.0, 0.0))
    assert (x1(-1, -1) == complex(0.0, 0.0))
    assert (x1(-1, -2) == complex(0.0, 0.0))
    assert (x1(-2, -1) == complex(1.0, 2.0))
    assert (x1(-2, -2) == complex(0.0, 0.0))

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

    x2 = cn(2)
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
    assert (x2(0, 0) == complex(0.0, 0.0))
    assert (x2(0, 1) == complex(0.0, 0.0))
    assert (x2(0, 2) == complex(1.0, 2.0))
    assert (x2(1, 0) == complex(0.0, 0.0))
    assert (x2(1, 1) == complex(0.0, 0.0))
    assert (x2(1, 2) == complex(3.0, 4.0))
    assert (x2(2, 0) == complex(0.0, 0.0))
    assert (x2(2, 1) == complex(0.0, 0.0))
    assert (x2(2, 2) == complex(0.0, 0.0))
    assert (x2(-1, -1) == complex(0.0, 0.0))
    assert (x2(-1, -2) == complex(0.0, 0.0))
    assert (x2(-1, -3) == complex(0.0, 0.0))
    assert (x2(-2, -1) == complex(3.0, 4.0))
    assert (x2(-2, -2) == complex(0.0, 0.0))
    assert (x2(-2, -3) == complex(0.0, 0.0))
    assert (x2(-3, -1) == complex(1.0, 2.0))
    assert (x2(-3, -2) == complex(0.0, 0.0))
    assert (x2(-3, -3) == complex(0.0, 0.0))

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


def test_cn_operator_bracket():
    from lielab.domain import cn

    xzero = cn.from_shape(0)
    xzero.set_vector([])

    # Out of bounds
    assert (np.isnan(np.real(xzero[-1])))
    assert (np.isnan(np.imag(xzero[-1])))
    assert (np.isnan(np.real(xzero[0])))
    assert (np.isnan(np.imag(xzero[0])))
    assert (np.isnan(np.real(xzero[1])))
    assert (np.isnan(np.imag(xzero[1])))

    x1 = cn(1)
    x1.set_vector([1.0, 2.0])

    # In bounds
    assert (x1[0] == complex(1.0, 2.0))
    assert (x1[-1] == complex(1.0, 2.0))

    # Out of bounds
    assert (np.isnan(np.real(x1[-2])))
    assert (np.isnan(np.imag(x1[-2])))
    assert (np.isnan(np.real(x1[1])))
    assert (np.isnan(np.imag(x1[1])))

    x2 = cn(2)
    x2.set_vector([1.0, 2.0, 3.0, 4.0])

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


def test_cn_math_ops_int():
    from lielab.domain import cn

    x1 = cn(2)
    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_lm_2 = 2*x1
    assert (x1_lm_2(0) == 2.5)
    assert (x1_lm_2(1) == 5.0)
    assert (x1_lm_2(2) == 7.5)
    assert (x1_lm_2(3) == 2.0)

    x1_rm_2 = x1*2
    assert (x1_rm_2(0) == 2.5)
    assert (x1_rm_2(1) == 5.0)
    assert (x1_rm_2(2) == 7.5)
    assert (x1_rm_2(3) == 2.0)

    x1 *= 2
    assert (x1(0) == 2.5)
    assert (x1(1) == 5.0)
    assert (x1(2) == 7.5)
    assert (x1(3) == 2.0)

    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_lm_2j = (2*1j)*x1
    assert (x1_lm_2j(0) == -5.0)
    assert (x1_lm_2j(1) == 2.5)
    assert (x1_lm_2j(2) == -2.0)
    assert (x1_lm_2j(3) == 7.5)

    x1_rm_2j = x1*(2*1j)
    assert (x1_rm_2j(0) == -5.0)
    assert (x1_rm_2j(1) == 2.5)
    assert (x1_rm_2j(2) == -2.0)
    assert (x1_rm_2j(3) == 7.5)

    x1 *= 2*1j
    assert (x1(0) == -5.0)
    assert (x1(1) == 2.5)
    assert (x1(2) == -2.0)
    assert (x1(3) == 7.5)

    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_d_2 = x1/2
    assert (x1_d_2(0) == 0.625)
    assert (x1_d_2(1) == 1.25)
    assert (x1_d_2(2) == 1.875)
    assert (x1_d_2(3) == 0.5)

    x1 /= 2
    assert (x1(0) == 0.625)
    assert (x1(1) == 1.25)
    assert (x1(2) == 1.875)
    assert (x1(3) == 0.5)

    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_d_2j = x1/(2*1j)
    assert (x1_d_2j(0) == 1.25)
    assert (x1_d_2j(1) == -0.625)
    assert (x1_d_2j(2) == 0.5)
    assert (x1_d_2j(3) == -1.875)

    x1 /= 2*1j
    assert (x1(0) == 1.25)
    assert (x1(1) == -0.625)
    assert (x1(2) == 0.5)
    assert (x1(3) == -1.875)


def test_cn_math_ops_double():
    from lielab.domain import cn

    x1 = cn(2)
    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_lm_2 = 2.0*x1
    assert (x1_lm_2(0) == 2.5)
    assert (x1_lm_2(1) == 5.0)
    assert (x1_lm_2(2) == 7.5)
    assert (x1_lm_2(3) == 2.0)

    x1_rm_2 = x1*2.0
    assert (x1_rm_2(0) == 2.5)
    assert (x1_rm_2(1) == 5.0)
    assert (x1_rm_2(2) == 7.5)
    assert (x1_rm_2(3) == 2.0)

    x1 *= 2.0
    assert (x1(0) == 2.5)
    assert (x1(1) == 5.0)
    assert (x1(2) == 7.5)
    assert (x1(3) == 2.0)

    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_lm_2j = (2.0*1j)*x1
    assert (x1_lm_2j(0) == -5.0)
    assert (x1_lm_2j(1) == 2.5)
    assert (x1_lm_2j(2) == -2.0)
    assert (x1_lm_2j(3) == 7.5)

    x1_rm_2j = x1*(2.0*1j)
    assert (x1_rm_2j(0) == -5.0)
    assert (x1_rm_2j(1) == 2.5)
    assert (x1_rm_2j(2) == -2.0)
    assert (x1_rm_2j(3) == 7.5)

    x1 *= 2.0*1j
    assert (x1(0) == -5.0)
    assert (x1(1) == 2.5)
    assert (x1(2) == -2.0)
    assert (x1(3) == 7.5)

    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_d_2 = x1/2.0
    assert (x1_d_2(0) == 0.625)
    assert (x1_d_2(1) == 1.25)
    assert (x1_d_2(2) == 1.875)
    assert (x1_d_2(3) == 0.5)

    x1 /= 2.0
    assert (x1(0) == 0.625)
    assert (x1(1) == 1.25)
    assert (x1(2) == 1.875)
    assert (x1(3) == 0.5)

    x1.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_d_2j = x1/(2.0*1j)
    assert (x1_d_2j(0) == 1.25)
    assert (x1_d_2j(1) == -0.625)
    assert (x1_d_2j(2) == 0.5)
    assert (x1_d_2j(3) == -1.875)

    x1 /= 2.0*1j
    assert (x1(0) == 1.25)
    assert (x1(1) == -0.625)
    assert (x1(2) == 0.5)
    assert (x1(3) == -1.875)


def test_cn_math_ops_cn():
    from lielab.domain import cn

    x1 = cn(2)
    x2 = cn(2)
    x1.set_vector([1.0, 2.0, 3.0, 4.0])
    x2.set_vector([1.25, 2.5, 3.75, 1.0])

    x1_add_x2 = x1 + x2
    assert (x1_add_x2(0) == 2.25)
    assert (x1_add_x2(1) == 4.5)
    assert (x1_add_x2(2) == 6.75)
    assert (x1_add_x2(3) == 5.0)

    x1 += x2
    assert (x1(0) == 2.25)
    assert (x1(1) == 4.5)
    assert (x1(2) == 6.75)
    assert (x1(3) == 5.0)

    x1.set_vector([1.0, 2.0, 3.0, 4.0])

    x1_sub_x2 = x1 - x2
    assert (x1_sub_x2(0) == -0.25)
    assert (x1_sub_x2(1) == -0.5)
    assert (x1_sub_x2(2) == -0.75)
    assert (x1_sub_x2(3) == 3.0)

    x1 -= x2
    assert (x1(0) == -0.25)
    assert (x1(1) == -0.5)
    assert (x1(2) == -0.75)
    assert (x1(3) == 3.0)

    x1.set_vector([1.0, 2.0, 3.0, 4.0])

    x1_unary_sub = (-x1)
    assert (x1_unary_sub(0) == -1.0)
    assert (x1_unary_sub(1) == -2.0)
    assert (x1_unary_sub(2) == -3.0)
    assert (x1_unary_sub(3) == -4.0)


def test_cn_from_vector():
    """
    Tests the from_vector operation.
    """

    from lielab.domain import cn

    x0 = cn.from_vector([])
    x0bar = x0.get_vector()

    assert (x0.get_shape() == 1)
    assert (x0bar.size == 0)

    x1 = cn.from_vector([1.0])
    x1bar = x1.get_vector()

    assert (x1.get_shape() == 2)
    assert (x1bar.size == 2)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 0.0)

    x2 = cn.from_vector([1.0, 2.0])
    x2bar = x2.get_vector()

    assert (x2.get_shape() == 2)
    assert (x2bar.size == 2)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)

    x3 = cn.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.get_vector()

    assert (x3.get_shape() == 3)
    assert (x3bar.size == 4)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)
    assert (x3bar[3] == 0.0)


def test_cn_from_to_complex_vector():
    """
    Tests the from_to_complex_vector operation.
    """

    from lielab.domain import cn

    x0 = cn.from_complex_vector([])
    x0bar = x0.to_complex_vector()

    assert (x0.get_shape() == 1)
    assert (x0bar.size == 0)

    x1 = cn.from_complex_vector([1.0])
    x1bar = x1.to_complex_vector()

    assert (x1.get_shape() == 2)
    assert (x1bar.size == 1)
    assert (x1bar[0] == complex(1.0, 0.0))

    x2 = cn.from_complex_vector([1.0 + 2.0*1j])
    x2bar = x2.to_complex_vector()

    assert (x2.get_shape() == 2)
    assert (x2bar.size == 1)
    assert (x2bar[0] == complex(1.0, 2.0))

    x3 = cn.from_complex_vector([1.0 + 2.0*1j, 3.0])
    x3bar = x3.to_complex_vector()

    assert (x3.get_shape() == 3)
    assert (x3bar.size == 2)
    assert (x3bar[0] == complex(1.0, 2.0))
    assert (x3bar[1] == complex(3.0, 0.0))


def test_cn_project():
    from lielab.domain import cn

    rand_2_2 = np.random.rand(2, 2) + 1j*np.random.rand(2, 2)
    proj_2_2 = cn.project(rand_2_2)

    assert (proj_2_2.shape[0] == 2)
    assert (proj_2_2.shape[1] == 2)
    assert (proj_2_2[0, 0] == complex(0.0, 0.0))
    assert (proj_2_2[0, 1] == rand_2_2[0, 1])
    assert (proj_2_2[1, 0] == complex(0.0, 0.0))
    assert (proj_2_2[1, 1] == complex(0.0, 0.0))

    rand_3_3 = np.random.rand(3, 3) + 1j*np.random.rand(3, 3)
    proj_3_3 = cn.project(rand_3_3)

    assert (proj_3_3.shape[0] == 3)
    assert (proj_3_3.shape[1] == 3)
    assert (proj_3_3[0, 0] == complex(0.0, 0.0))
    assert (proj_3_3[0, 1] == complex(0.0, 0.0))
    assert (proj_3_3[0, 2] == rand_3_3[0, 2])
    assert (proj_3_3[1, 0] == complex(0.0, 0.0))
    assert (proj_3_3[1, 1] == complex(0.0, 0.0))
    assert (proj_3_3[1, 2] == rand_3_3[1, 2])
    assert (proj_3_3[2, 0] == complex(0.0, 0.0))
    assert (proj_3_3[2, 1] == complex(0.0, 0.0))
    assert (proj_3_3[2, 2] == complex(0.0, 0.0))

    rand_2_3 = np.random.rand(2, 3) + 1j*np.random.rand(2, 3)
    proj_2_3 = cn.project(rand_2_3)

    assert (proj_2_3.shape[0] == 2)
    assert (proj_2_3.shape[1] == 2)
    assert (proj_2_3[0, 0] == complex(0.0, 0.0))
    assert (proj_2_3[0, 1] == rand_2_3[0, 1])
    assert (proj_2_3[1, 0] == complex(0.0, 0.0))
    assert (proj_2_3[1, 1] == complex(0.0, 0.0))

    rand_3_2 = np.random.rand(3, 2) + 1j*np.random.rand(3, 2)
    proj_3_2 = cn.project(rand_3_2)

    assert (proj_3_2.shape[0] == 2)
    assert (proj_3_2.shape[1] == 2)
    assert (proj_3_2[0, 0] == complex(0.0, 0.0))
    assert (proj_3_2[0, 1] == rand_3_2[0, 1])
    assert (proj_3_2[1, 0] == complex(0.0, 0.0))
    assert (proj_3_2[1, 1] == complex(0.0, 0.0))
