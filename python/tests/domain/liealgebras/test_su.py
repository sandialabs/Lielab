from lielab.testing import *
import pytest

def complex(a,b):
    return a + b*1j

def test_su_to_string():
    from lielab.domain import su

    x0 = su(0)
    assert (x0.to_string() == "su(0)")
    x1 = su(1)
    assert (x1.to_string() == "su(1)")
    x10 = su(10)
    assert (x10.to_string() == "su(10)")

def test_su_main_initializer():
    from lielab.domain import su

    xblank = su()
    assert (xblank.get_dimension() == 0)

    x0 = su(0)
    assert (x0.get_dimension() == 0)
    x1 = su(1)
    assert (x1.get_dimension() == 0)
    x10 = su(2)
    assert (x10.get_dimension() == 3)

def test_su_matrix_initializer():
    from lielab.domain import su

    x0 = su(np.random.rand(0, 0) + 1j*np.random.rand(0, 0))
    assert (x0.get_shape() == 0)

    x1 = su(np.random.rand(1, 1) + 1j*np.random.rand(1, 1))
    assert (x1.get_shape() == 1)

    x2 = su(np.random.rand(2, 2) + 1j*np.random.rand(2, 2))
    assert (x2.get_shape() == 2)

    with pytest.raises(RuntimeError):
        su(np.random.rand(2, 3) + 1j*np.random.rand(2, 3))
    
    with pytest.raises(RuntimeError):
        su(np.random.rand(3, 2) + 1j*np.random.rand(3, 2))

def test_su_basis_initializer():
    from lielab.domain import su

    xm10 = su.basis(-1, 0)
    assert (xm10.get_dimension() == 0)
    xm10bar = xm10.get_vector()
    assert (xm10bar.size == 0)

    x00 = su.basis(0, 0)
    assert (x00.get_dimension() == 0)
    x00bar = x00.get_vector()
    assert (x00bar.size == 0)

    x10 = su.basis(1, 0)
    assert (x10.get_dimension() == 0)
    x10bar = x10.get_vector()
    assert (x10bar.size == 0)

    x02 = su.basis(0, 2)
    assert (x02.get_dimension() == 3)
    x02bar = x02.get_vector()
    assert (x02bar.size == 3)
    assert (x02bar[0] == 1.0)
    assert (x02bar[1] == 0.0)
    assert (x02bar[2] == 0.0)

    x12 = su.basis(1, 2)
    assert (x12.get_dimension() == 3)
    x12bar = x12.get_vector()
    assert (x12bar.size == 3)
    assert (x12bar[0] == 0.0)
    assert (x12bar[1] == 1.0)
    assert (x12bar[2] == 0.0)

    x22 = su.basis(2, 2)
    assert (x22.get_dimension() == 3)
    x22bar = x22.get_vector()
    assert (x22bar.size == 3)
    assert (x22bar[0] == 0.0)
    assert (x22bar[1] == 0.0)
    assert (x22bar[2] == 1.0)

    x32 = su.basis(3, 2)
    assert (x32.get_dimension() == 3)
    x32bar = x32.get_vector()
    assert (x32bar.size == 3)
    assert (x32bar[0] == 0.0)
    assert (x32bar[1] == 0.0)
    assert (x32bar[2] == 0.0)

    x03 = su.basis(0, 3)
    assert (x03.get_dimension() == 8)
    x03bar = x03.get_vector()
    assert (x03bar.size == 8)
    assert (x03bar[0] == 1.0)
    assert (x03bar[1] == 0.0)
    assert (x03bar[2] == 0.0)
    assert (x03bar[3] == 0.0)
    assert (x03bar[4] == 0.0)
    assert (x03bar[5] == 0.0)
    assert (x03bar[6] == 0.0)
    assert (x03bar[7] == 0.0)

def test_su_from_shape_initializer():
    from lielab.domain import su

    x0 = su.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = su.from_shape(1)
    assert (x1.get_dimension() == 0)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = su.from_shape(2)
    assert (x2.get_dimension() == 3)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_su_get_dimension():
    from lielab.domain import su

    zero = su(0)
    one = su(1)
    two = su(2)
    three = su(3)
    four = su(4)
    five = su(5)
    six = su(6)
    seven = su(7)
    eight = su(8)

    assert (zero.get_dimension() == 0)
    assert (one.get_dimension() == 0)
    assert (two.get_dimension() == 3)
    assert (three.get_dimension() == 8)
    assert (four.get_dimension() == 15)
    assert (five.get_dimension() == 24)
    assert (six.get_dimension() == 35)
    assert (seven.get_dimension() == 48)
    assert (eight.get_dimension() == 63)

def test_su_set_get_vector():
    """
    * Tests the set/get_vector operation.
    """

    from lielab.domain import su

    x0 = su(0)
    x0.set_vector([])
    x0bar = x0.get_vector()

    assert (x0bar.size == 0)

    x2 = su(2)
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

    x2.set_vector([8.0, 9.0])
    x2bar = x2.get_vector()

    assert (x2bar.size == 3)
    assert (x2bar[0] == 8.0)
    assert (x2bar[1] == 9.0)
    assert (x2bar[2] == 7.0)

    x3 = su(3)
    x3.set_vector([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 8)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)
    assert (x3bar[3] == 4.0)
    assert (x3bar[4] == 5.0)
    assert (x3bar[5] == 6.0)
    assert (x3bar[6] == 7.0)
    assert (x3bar[7] == 0.0)

    x3.set_vector([8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 8)
    assert (x3bar[0] == 8.0)
    assert (x3bar[1] == 9.0)
    assert (x3bar[2] == 10.0)
    assert (x3bar[3] == 11.0)
    assert (x3bar[4] == 12.0)
    assert (x3bar[5] == 13.0)
    assert (x3bar[6] == 14.0)
    assert (x3bar[7] == 15.0)

    x3.set_vector([16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0])
    x3bar = x3.get_vector()

    assert (x3bar.size == 8)
    assert (x3bar[0] == 16.0)
    assert (x3bar[1] == 17.0)
    assert (x3bar[2] == 18.0)
    assert (x3bar[3] == 19.0)
    assert (x3bar[4] == 20.0)
    assert (x3bar[5] == 21.0)
    assert (x3bar[6] == 22.0)
    assert (x3bar[7] == 23.0)

def test_su_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import su

    x0 = su(0)
    x0.set_vector([])
    x0hat = x0.get_matrix()

    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = su(1)
    x1.set_vector([])
    x1hat = x1.get_matrix()

    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)
    assert (x1hat[0, 0] == complex(0.0, 0.0))

    x2 = su(2)
    x2.set_vector([1.0, 2.0, 3.0])
    x2hat = x2.get_matrix()

    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == complex(0.0, 3.0))
    assert (x2hat[0, 1] == complex(-2.0, 1.0))
    assert (x2hat[1, 0] == complex(2.0, 1.0))
    assert (x2hat[1, 1] == complex(0.0, -3.0))

    x2.set_vector([4.0, 5.0])
    x2hat = x2.get_matrix()

    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == complex(0.0, 3.0))
    assert (x2hat[0, 1] == complex(-5.0, 4.0))
    assert (x2hat[1, 0] == complex(5.0, 4.0))
    assert (x2hat[1, 1] == complex(0.0, -3.0))

    x2.set_vector([6.0])
    x2hat = x2.get_matrix()

    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)
    assert (x2hat[0, 0] == complex(0.0, 3.0))
    assert (x2hat[0, 1] == complex(-5.0, 6.0))
    assert (x2hat[1, 0] == complex(5.0, 6.0))
    assert (x2hat[1, 1] == complex(0.0, -3.0))

def test_su_operator_parenthesis():
    from lielab.domain import su

    x0 = su(0)
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

    x1 = su(1)
    x1.set_vector([])

    # Out of bounds
    assert (np.isnan(x1(-1)))
    assert (np.isnan(x1(0)))
    assert (np.isnan(x1(1)))

    # In bounds
    assert (x1(0, 0) == complex(0.0, 0.0))
    assert (x1(-1, -1) == complex(0.0, 0.0))

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

    x2 = su(2)
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
    assert (x2(0, 0) == complex(0.0, 3.0))
    assert (x2(0, 1) == complex(-2.0, 1.0))
    assert (x2(1, 0) == complex(2.0, 1.0))
    assert (x2(1, 1) == complex(0.0, -3.0))
    assert (x2(-1, -1) == complex(0.0, -3.0))
    assert (x2(-1, -2) == complex(2.0, 1.0))
    assert (x2(-2, -1) == complex(-2.0, 1.0))
    assert (x2(-2, -2) == complex(0.0, 3.0))

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

# TODO: Math ops int

def test_su_math_ops_double():
    from lielab.domain import su

    x1 = su(2)
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
    # TODO: Imaginary ops move elements out of the algebra. Dunno how this should be handled.
    x1_lm_2j = (2.0*1.0j)*x1
    assert (x1_lm_2j(0) == 0.0)
    assert (x1_lm_2j(1) == 0.0)
    assert (x1_lm_2j(2) == 0.0)

    x1_rm_2j = x1*(2.0*1.0j)
    assert (x1_rm_2j(0) == 0.0)
    assert (x1_rm_2j(1) == 0.0)
    assert (x1_rm_2j(2) == 0.0)

    x1 *= 2.0*1.0j
    assert (x1(0) == 0.0)
    assert (x1(1) == 0.0)
    assert (x1(2) == 0.0)

    x1.set_vector([1.25, 2.5, 3.75])

    x1_d_2 = x1/2.0
    assert (x1_d_2(0) == 0.625)
    assert (x1_d_2(1) == 1.25)
    assert (x1_d_2(2) == 1.875)

    x1 /= 2.0
    assert (x1(0) == 0.625)
    assert (x1(1) == 1.25)
    assert (x1(2) == 1.875)

    x1.set_vector([1.25, 2.5, 3.75])

    x1_d_2j = x1/(2.0*1.0j)
    assert (x1_d_2j(0) == 0.0)
    assert (x1_d_2j(1) == 0.0)
    assert (x1_d_2j(2) == 0.0)

    x1 /= 2.0*1.0j
    assert (x1(0) == 0.0)
    assert (x1(1) == 0.0)
    assert (x1(2) == 0.0)

def test_su_math_ops_cn():
    from lielab.domain import su

    x1 = su(2)
    x2 = su(2)
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

def test_su_from_vector():
    """
    * Tests the from_vector operation.
    """

    from lielab.domain import su

    x0 = su.from_vector([])
    x0bar = x0.get_vector()

    assert (x0.get_shape() == 1)
    assert (x0bar.size == 0)

    x1 = su.from_vector([1.0])
    x1bar = x1.get_vector()

    assert (x1.get_shape() == 2)
    assert (x1bar.size == 3)
    assert (x1bar[0] == 1.0)
    assert (x1bar[1] == 0.0)
    assert (x1bar[2] == 0.0)

    x2 = su.from_vector([1.0, 2.0])
    x2bar = x2.get_vector()

    assert (x2.get_shape() == 2)
    assert (x2bar.size == 3)
    assert (x2bar[0] == 1.0)
    assert (x2bar[1] == 2.0)
    assert (x2bar[2] == 0.0)

    x3 = su.from_vector([1.0, 2.0, 3.0])
    x3bar = x3.get_vector()

    assert (x3.get_shape() == 2)
    assert (x3bar.size == 3)
    assert (x3bar[0] == 1.0)
    assert (x3bar[1] == 2.0)
    assert (x3bar[2] == 3.0)

# TODO: Project

def test_su2():
    """
    * Tests the algebra with su(2).
    """

    from lielab.domain import su
    from lielab.functions import commutator

    D = su.basis(0, 2).get_dimension()

    # Construct the su2 basis
    basis = []
    for ii in range(D):
        basis.append(su.basis(ii, 2))

    # su2 specific identities
    assert_domain(commutator(basis[0], basis[1]),  2*basis[2])
    assert_domain(commutator(basis[1], basis[2]),  2*basis[0])
    assert_domain(commutator(basis[2], basis[0]),  2*basis[1])
    assert_domain(commutator(basis[1], basis[0]), -2*basis[2])
    assert_domain(commutator(basis[2], basis[1]), -2*basis[0])
    assert_domain(commutator(basis[0], basis[2]), -2*basis[1])

    # Hamilton's identities
    # Note that i^2 = j^2 = k^2 = -1^2 isn't checked since this isn't true for the algebra

    # ij = -ji = k
    # assert_domain(b[0]*b[1], -b[1]*b[0])
    # assert_domain(b[0]*b[1],  b[2])

    # # jk = -kj = i
    # assert_domain(b[1]*b[2], -b[2]*b[1])
    # assert_domain(b[1]*b[2],  b[0])

    # # ki = -ik = j
    # assert_domain(b[2]*b[0], -b[0]*b[2])
    # assert_domain(b[2]*b[0],  b[1])

def test_su3():
    """
    * Tests the algebra with su(3).
    """

    from lielab.domain import su

    D = su.basis(0, 3).get_dimension()

    # Construct the su3 basis
    basis = []
    for ii in range(D):
        basis.append(su.basis(ii, 3))

    # TODO: Implement these. It's tough to work out what these should be with GGM
    # su3 specific identities
    # assert_domain( commutator(t1, t2), _i*t3)
    # assert_domain( commutator(t1, t4), _i*t7 / 2.0)
    # assert_domain(-commutator(t1, t5), _i*t6 / 2.0)
    # assert_domain( commutator(t2, t4), _i*t6 / 2.0)
    # assert_domain( commutator(t2, t5), _i*t7 / 2.0)
    # assert_domain( commutator(t3, t4), _i*t5 / 2.0)
    # assert_domain(-commutator(t3, t6), _i*t7 / 2.0)
    # TODO: Why aren't these evaluating true?
    # assert_domain( commutator(t4, t5), _i*t8 * std::sqrt(3.0) / 2.0)
    # assert_domain( commutator(t6, t7), _i*t8 * std::sqrt(3.0) / 2.0)
