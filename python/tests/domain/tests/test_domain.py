"""
Test cases for the domain Python bindings.
"""

import pytest

from lielab.testing import *

def test_cn():
    """
    Checks the basic cn class
    """

    from lielab.domain import cn
    from lielab.functions import commutator

    one = cn(1)
    two = cn(2)
    three = cn(3)
    four = cn(4)
    five = cn(5)
    six = cn(6)
    seven = cn(7)
    eight = cn(8)

    # Dimensions
    assert one.get_dimension() == 0
    assert two.get_dimension() == 2
    assert three.get_dimension() == 4
    assert four.get_dimension() == 6
    assert five.get_dimension() == 8
    assert six.get_dimension() == 10
    assert seven.get_dimension() == 12
    assert eight.get_dimension() == 14

    x = cn.basis(0,4)
    y = cn.basis(1,4)
    z = cn.basis(2,4)
    zero = x*0

    # Unary subtraction
    -x

    # Vector Addition
    x + y
    assert_domain(x + y, y + x)

    # Vector Subtraction
    x - y
    assert_domain(x - y, -(y - x))

    # Scalar Multiplication
    an_int*x
    a_double*x
    x*an_int
    x*a_double
    assert_domain(an_int*x, x*an_int)
    assert_domain(a_double*x, x*a_double)

    # Scalar division
    x/an_int
    x/a_double

    # Vector multiplication
    x * y
    assert_domain(commutator(x,y), -commutator(y,x))
    assert_domain(commutator(x,y), commutator(y,x)) # Abelian check


def test_gl():
    """
    Checks the basic gl class
    """

    from lielab.domain import gl
    from lielab.functions import commutator

    one = gl(1)
    two = gl(2)
    three = gl(3)
    four = gl(4)
    five = gl(5)
    six = gl(6)
    seven = gl(7)
    eight = gl(8)

    # Dimensions
    assert one.get_dimension() == 1
    assert two.get_dimension() == 4
    assert three.get_dimension() == 9
    assert four.get_dimension() == 16
    assert five.get_dimension() == 25
    assert six.get_dimension() == 36
    assert seven.get_dimension() == 49
    assert eight.get_dimension() == 64

    x = gl.basis(0,4)
    y = gl.basis(1,4)
    z = gl.basis(2,4)
    zero = x*0

    # Unary subtraction
    -x

    # Vector Addition
    x + y
    assert_domain(x + y, y + x)

    # Vector Subtraction
    x - y
    assert_domain(x - y, -(y - x))

    # Scalar Multiplication
    an_int*x
    a_double*x
    x*an_int
    x*a_double
    assert_domain(an_int*x, x*an_int)
    assert_domain(a_double*x, x*a_double)

    # Scalar division
    x/an_int
    x/a_double

    # Vector multiplication
    x * y
    assert_domain(commutator(x,y), -commutator(y,x))


def test_glc():
    """
    Checks the basic glc class
    """

    from lielab.domain import glc
    from lielab.functions import commutator

    one = glc(1)
    two = glc(2)
    three = glc(3)
    four = glc(4)
    five = glc(5)
    six = glc(6)
    seven = glc(7)
    eight = glc(8)

    # Dimensions
    assert one.get_dimension() == 2
    assert two.get_dimension() == 8
    assert three.get_dimension() == 18
    assert four.get_dimension() == 32
    assert five.get_dimension() == 50
    assert six.get_dimension() == 72
    assert seven.get_dimension() == 98
    assert eight.get_dimension() == 128

    x = glc.basis(0,4)
    y = glc.basis(1,4)
    z = glc.basis(2,4)
    zero = x*0

    # Unary subtraction
    -x

    # Vector Addition
    x + y
    assert_domain(x + y, y + x)

    # Vector Subtraction
    x - y
    assert_domain(x - y, -(y - x))

    # Scalar Multiplication
    an_int*x
    a_double*x
    x*an_int
    x*a_double
    assert_domain(an_int*x, x*an_int)
    assert_domain(a_double*x, x*a_double)

    # Scalar division
    x/an_int
    x/a_double

    # Vector multiplication
    x * y
    assert_domain(commutator(x,y), -commutator(y,x))


def test_rn():
    """
    Checks the basic rn class
    """

    from lielab.domain import rn
    from lielab.functions import commutator

    one = rn(1)
    two = rn(2)
    three = rn(3)
    four = rn(4)
    five = rn(5)
    six = rn(6)
    seven = rn(7)
    eight = rn(8)

    # Dimensions
    assert one.get_dimension() == 0
    assert two.get_dimension() == 1
    assert three.get_dimension() == 2
    assert four.get_dimension() == 3
    assert five.get_dimension() == 4
    assert six.get_dimension() == 5
    assert seven.get_dimension() == 6
    assert eight.get_dimension() == 7

    x = rn.basis(0,4)
    y = rn.basis(1,4)
    z = rn.basis(2,4)
    zero = x*0

    # Unary subtraction
    -x

    # Vector Addition
    x + y
    assert_domain(x + y, y + x)

    # Vector Subtraction
    x - y
    assert_domain(x - y, -(y - x))

    # Scalar Multiplication
    an_int*x
    a_double*x
    x*an_int
    x*a_double
    assert_domain(an_int*x, x*an_int)
    assert_domain(a_double*x, x*a_double)

    # Scalar division
    x/an_int
    x/a_double

    # Vector multiplication
    x * y
    assert_domain(commutator(x,y), -commutator(y,x))
    assert_domain(commutator(x,y), commutator(y,x)) # Abelian check


def test_se():
    """
    Checks the basic se class
    """

    from lielab.domain import se
    from lielab.functions import commutator

    two = se(2)
    three = se(3)
    four = se(4)
    five = se(5)
    six = se(6)
    seven = se(7)

    # Dimensions
    assert two.get_dimension() == 1
    assert three.get_dimension() == 3
    assert four.get_dimension() == 6
    assert five.get_dimension() == 10
    assert six.get_dimension() == 15
    assert seven.get_dimension() == 21

    x = se.basis(0, 4)
    y = se.basis(1, 4)
    z = se.basis(2, 4)
    u = se.basis(3, 4)
    v = se.basis(4, 4)
    w = se.basis(5, 4)
    zero = x*0

    # Unary subtraction
    -x

    # Vector Addition
    x + y
    assert_domain(x + y, y + x)

    # Vector Subtraction
    x - y
    assert_domain(x - y, -(y - x))

    # Scalar Multiplication
    an_int*x
    a_double*x
    x*an_int
    x*a_double
    assert_domain(an_int*x, x*an_int)
    assert_domain(a_double*x, x*a_double)

    # Scalar division
    x/an_int
    x/a_double

    # Vector multiplication
    x * y
    assert_domain(commutator(u,v), -commutator(v,u))


def test_se4():
    """
    Tests the se algebra with se(4)
    """

    from lielab.domain import se
    from lielab.functions import commutator

    x = se.basis(0, 4)
    y = se.basis(1, 4)
    z = se.basis(2, 4)
    u = se.basis(3, 4)
    v = se.basis(4, 4)
    w = se.basis(5, 4)
    zero = se(4)

    # Vector addition
    x + y
    x + z
    x + u
    x + v
    x + w
    y + z
    y + u
    y + v
    y + w
    z + u
    z + v
    z + w
    u + v
    u + w
    v + w

    # Vector subtraction
    x - y
    x - z
    x - u
    x - v
    x - w
    y - z

    # Right distributive
    assert_domain((x + y) * z, x * z + y * z)

    # Left distributive
    assert_domain(x * (y + z), x * y + x * z)

    # Scalar multiplication
    assert_domain((a_double * x) * (a_double * y), (a_double * a_double) * (x * y))

    # Bilinearity
    assert_domain(commutator(a_double * x + a_double * y, z), a_double * commutator(x, z) + a_double * commutator(y, z))

    # Alternating
    assert_domain(commutator(x, x), zero)
    assert_domain(commutator(y, y), zero)
    assert_domain(commutator(z, z), zero)
    assert_domain(commutator(u, u), zero)
    assert_domain(commutator(v, v), zero)
    assert_domain(commutator(w, w), zero)

    # Jacobi identity
    assert_domain(commutator(x, commutator(y, z)) + commutator(z, commutator(x, y)) + commutator(y, commutator(z, x)), zero)

    # Anticommutivity
    assert_domain(commutator(x, y), -commutator(y, x))
    assert_domain(commutator(x, z), -commutator(z, x))
    assert_domain(commutator(y, z), -commutator(z, y))
    assert_domain(commutator(u, v), -commutator(v, u))
    assert_domain(commutator(u, w), -commutator(w, u))
    assert_domain(commutator(v, w), -commutator(w, v))

    # se3 specific identities
    assert_domain(commutator(x, y), zero)
    assert_domain(commutator(y, z), zero)
    assert_domain(commutator(z, x), zero)
    assert_domain(commutator(y, x), zero)
    assert_domain(commutator(z, y), zero)
    assert_domain(commutator(x, z), zero)
    assert_domain(commutator(x, v), z)
    assert_domain(commutator(y, w), x)
    assert_domain(commutator(z, u), y)
    assert_domain(commutator(y, u), -z)
    assert_domain(commutator(z, v), -x)
    assert_domain(commutator(x, w), -y)
    assert_domain(commutator(u, v), w)
    assert_domain(commutator(v, w), u)
    assert_domain(commutator(w, u), v)
    assert_domain(commutator(v, u), -w)
    assert_domain(commutator(w, v), -u)
    assert_domain(commutator(u, w), -v)

    # se3 specific vectors
    assert abs(u(0,0) - 0.0) < TOL_FINE
    assert abs(u(0,1) - 0.0) < TOL_FINE
    assert abs(u(0,2) - 0.0) < TOL_FINE
    assert abs(u(1,0) - 0.0) < TOL_FINE
    assert abs(u(1,1) - 0.0) < TOL_FINE
    assert abs(u(1,2) + 1.0) < TOL_FINE
    assert abs(u(2,0) - 0.0) < TOL_FINE
    assert abs(u(2,1) - 1.0) < TOL_FINE
    assert abs(u(2,2) - 0.0) < TOL_FINE


def test_so():
    """
    Checks the basic so class
    """

    from lielab.domain import so
    from lielab.functions import commutator

    one = so(1)
    four = so(4)
    five = so(5)
    six = so(6)
    seven = so(7)
    eight = so(8)

    # Dimensions
    assert one.get_dimension() == 0
    # so2 is checked independendly
    # so3 is checked independently
    assert four.get_dimension() == 6
    assert five.get_dimension() == 10
    assert six.get_dimension() == 15
    assert seven.get_dimension() == 21
    assert eight.get_dimension() == 28

    x = so.basis(0,3)
    y = so.basis(1,3)
    z = so.basis(2,3)
    zero = x*0

    # Unary subtraction
    -x

    # Vector Addition
    x + y
    assert_domain(x + y, y + x)

    # Vector Subtraction
    x - y
    assert_domain(x - y, -(y - x))

    # Scalar Multiplication
    an_int*x
    a_double*x
    x*an_int
    x*a_double
    assert_domain(an_int*x, x*an_int)
    assert_domain(a_double*x, x*a_double)

    # Scalar division
    x/an_int
    x/a_double

    # Vector multiplication
    x * y
    assert_domain(commutator(x,y), -commutator(y,x))


def test_so2():
    """
    Tests the so algebra with so(2).
    
    Note that some of these test cases are trivial since so(2) is 1-dimensional.
    """

    from lielab.domain import so
    from lielab.functions import commutator

    x = so.basis(0,2)
    zero = x*0

    # Check dimension
    assert x.get_dimension() == 1

    # Vector addition
    x + x
    assert_domain(x + x, 2*x)

    # Vector subtraction
    x - x
    assert_domain(x - x, zero)

    # Right distributive
    assert_domain((x + x) * x, x * x + x * x)

    # Left distributive
    assert_domain(x * (x + x), x * x + x * x)

    # Scalar multiplication
    b = 3.0
    assert_domain((a_double * x) * (b * x), (a_double * b) * (x * x))

    # Bilinearity
    assert_domain(commutator(a_double * x + b * x, x), a_double * commutator(x, x) + b * commutator(x, x))

    # Alternating
    assert_domain(commutator(x, x), zero)

    # Jacobi Identity
    assert_domain(commutator(x, commutator(x, x)) + commutator(x, commutator(x, x)) + commutator(x, commutator(x, x)), zero)

    # Anticommutivity
    assert_domain(commutator(x, x), -commutator(x, x))


def test_so3():
    """
    Tests the so algebra with so(3)
    """

    from lielab.domain import so
    from lielab.functions import commutator

    zero = so(3)
    x = so.basis(0,3)
    y = so.basis(1,3)
    z = so.basis(2,3)

    zero.set_vector([0,0,0])
    x.set_vector([1,0,0])
    y.set_vector([0,1,0])
    z.set_vector([0,0,1])

    # Vector addition
    x + y
    x + z
    y + z

    # Vector subtraction
    x - y
    x - z
    y - z

    # Right distributive
    assert_domain((x + y) * z, x * z + y * z)

    # Left distributive
    assert_domain(x * (y + z), x * y + x * z)

    # Scalar multiplication
    assert_domain((a_double * x) * (a_double * y), (a_double * a_double) * (x * y))

    # Bilinearity
    assert_domain(commutator(a_double * x + a_double * y, z), a_double * commutator(x, z) + a_double * commutator(y, z))

    # Alternating
    assert_domain(commutator(x, x), zero)
    assert_domain(commutator(y, y), zero)
    assert_domain(commutator(z, z), zero)

    # Jacobi identity
    assert_domain(commutator(x, commutator(y, z)) + commutator(z, commutator(x, y)) + commutator(y, commutator(z, x)), zero)

    # Anticommutivity
    assert_domain(commutator(x, y), -commutator(y, x))
    assert_domain(commutator(x, z), -commutator(z, x))
    assert_domain(commutator(y, z), -commutator(z, y))

    # so3 specific identities
    assert_domain(commutator(x, y), z)
    assert_domain(commutator(y, z), x)
    assert_domain(commutator(z, x), y)
    assert_domain(commutator(y, x), -z)
    assert_domain(commutator(z, y), -x)
    assert_domain(commutator(x, z), -y)

    # so3 specific vectors
    assert abs(x(0,0) - 0.0) < TOL_FINE
    assert abs(x(0,1) - 0.0) < TOL_FINE
    assert abs(x(0,2) - 0.0) < TOL_FINE
    assert abs(x(1,0) - 0.0) < TOL_FINE
    assert abs(x(1,1) - 0.0) < TOL_FINE
    assert abs(x(1,2) + 1.0) < TOL_FINE
    assert abs(x(2,0) - 0.0) < TOL_FINE
    assert abs(x(2,1) - 1.0) < TOL_FINE
    assert abs(x(2,2) - 0.0) < TOL_FINE

def test_sp():
    """
    Tests the sp algebra.
    """

    from lielab.domain import sp
    from lielab.functions import commutator
   
    two = sp(2)
    four = sp(4)
    six = sp(6)
    eight = sp(8)

    # Dimensions
    assert two.get_dimension() == 3
    assert four.get_dimension() == 10
    assert six.get_dimension() == 21
    assert eight.get_dimension() == 36

    x = sp.basis(0,2)
    y = sp.basis(1,2)
    z = sp.basis(2,2)
    zero = x*0

    # Unary subtraction
    -x

    # Vector Addition
    x + y
    assert_domain(x + y, y + x)

    # Vector Subtraction
    x - y
    assert_domain(x - y, -(y - x))

    # Scalar Multiplication
    an_int*x
    a_double*x
    x*an_int
    x*a_double
    assert_domain(an_int*x, x*an_int)
    assert_domain(a_double*x, x*a_double)

    # Scalar division
    x/an_int
    x/a_double

    # Vector multiplication
    x * y
    assert_domain(commutator(x,y), -commutator(y,x))


def test_su():
    """
    Tests the su class
    """

    from lielab.domain import su
    from lielab.functions import commutator

    u = su(2)
    v = su(2)
    w = su(2)

    zero = u*0

    # Check dimension
    assert u.get_dimension() == 3
    assert v.get_dimension() == 3
    assert w.get_dimension() == 3

    u.set_vector([1,0,0])
    v.set_vector([0,1,0])
    w.set_vector([0,0,1])

    # Unary subtraction
    -u

    # Vector Addition
    u + v
    assert_domain(u + v, v + u)

    # Vector Subtraction
    u - v
    assert_domain(u - v, -(v - u))

    # Scalar Multiplication
    an_int*u
    a_double*u
    u*an_int
    u*a_double
    assert_domain(an_int*u, u*an_int)
    assert_domain(a_double*u, u*a_double)

    # Vector multiplication
    u * v
    assert_domain(commutator(u,v), -commutator(v,u))


def test_su2():
    """
    Tests the su algebra with su(2).
    """

    from lielab.domain import su
    from lielab.functions import commutator

    zero = su(2)
    u = su(2)
    v = su(2)
    w = su(2)

    assert u.get_dimension() == 3
    assert v.get_dimension() == 3
    assert w.get_dimension() == 3

    u.set_vector([1,0,0])
    v.set_vector([0,1,0])
    w.set_vector([0,0,1])

    # Right distributive
    assert_domain((u + v) * w, u * w + v * w)

    # Left distributive
    assert_domain(u * (v + w), u * v + u * w)

    # Scalar multiplication
    assert_domain((a_double * u) * (a_double * v), (a_double * a_double) * (u * v))

    # Bilinearity
    assert_domain(commutator(a_double * u + a_double * v, w), a_double * commutator(u, w) + a_double * commutator(v, w))

    # Alternating
    assert_domain(commutator(u, u), zero)
    assert_domain(commutator(v, v), zero)
    assert_domain(commutator(w, w), zero)

    # Jacobi identity
    assert_domain(commutator(u, commutator(v, w)) + commutator(w, commutator(u, v)) + commutator(v, commutator(w, u)), zero)

    # Anticommutivity
    assert_domain(commutator(u, v), -commutator(v, u))

    # su2 specific identities
    assert_domain(commutator(u, v), 2*w)
    assert_domain(commutator(v, w), 2*u)
    assert_domain(commutator(w, u), 2*v)
    assert_domain(commutator(v, u), -2*w)
    assert_domain(commutator(w, v), -2*u)
    assert_domain(commutator(u, w), -2*v)

    # Hamilton's identities
    # Note that i^2 = j^2 = k^2 = -1^2 isn't checked since there is no identity

    # ij = -ji = k
    assert_domain(u*v, -v*u)
    assert_domain(u*v, w)

    # jk = -kj = i
    assert_domain(v*w, -w*v)
    assert_domain(v*w, u)

    # ki = -ik = j
    assert_domain(w*u, -u*w)
    assert_domain(w*u, v)


def test_su3():
    """
    Tests the su algebra with su(3).
    """

    from lielab.domain import su
    from lielab.functions import commutator

    zero = su(3)
    t1 = su(3)
    t2 = su(3)
    t3 = su(3)
    t4 = su(3)
    t5 = su(3)
    t6 = su(3)
    t7 = su(3)
    t8 = su(3)

    assert t1.get_dimension() == 8
    assert t2.get_dimension() == 8
    assert t3.get_dimension() == 8
    assert t4.get_dimension() == 8
    assert t5.get_dimension() == 8
    assert t6.get_dimension() == 8
    assert t7.get_dimension() == 8
    assert t8.get_dimension() == 8

    t1.set_vector([1, 0, 0, 0, 0, 0, 0, 0])
    t2.set_vector([0, 1, 0, 0, 0, 0, 0, 0])
    t3.set_vector([0, 0, 1, 0, 0, 0, 0, 0])
    t4.set_vector([0, 0, 0, 1, 0, 0, 0, 0])
    t5.set_vector([0, 0, 0, 0, 1, 0, 0, 0])
    t6.set_vector([0, 0, 0, 0, 0, 1, 0, 0])
    t7.set_vector([0, 0, 0, 0, 0, 0, 1, 0])
    t8.set_vector([0, 0, 0, 0, 0, 0, 0, 1])

    # Right distributive
    assert_domain((t1 + t2) * t3, t1 * t3 + t2 * t3)

    # Left distributive
    assert_domain(t1 * (t2 + t3), t1 * t2 + t1 * t3)

    # Scalar multiplication
    assert_domain((a_double * t1) * (a_double * t2), (a_double * a_double) * (t1 * t2))

    # Bilinearity
    assert_domain(commutator(a_double * t1 + a_double * t2, t3), a_double * commutator(t1, t3) + a_double * commutator(t2, t3))

    # Alternating
    assert_domain(commutator(t1, t1), zero)
    assert_domain(commutator(t2, t2), zero)
    assert_domain(commutator(t2, t2), zero)

    # Jacobi identity
    assert_domain(commutator(t1, commutator(t2, t3)) + commutator(t3, commutator(t1, t2)) + commutator(t2, commutator(t3, t1)), zero)

    # Anticommutivity
    assert_domain(commutator(t1, t2), -commutator(t2, t1))

    # su3 specific identities
    # assert_domain( commutator(t1, t2), 1j*t3)
    # assert_domain( commutator(t1, t4), 1j*t7 / 2.0)
    # assert_domain(-commutator(t1, t5), 1j*t6 / 2.0)
    # assert_domain( commutator(t2, t4), 1j*t6 / 2.0)
    # assert_domain( commutator(t2, t5), 1j*t7 / 2.0)
    # assert_domain( commutator(t3, t4), 1j*t5 / 2.0)
    # assert_domain(-commutator(t3, t6), 1j*t7 / 2.0)


def test_CN():
    """
    Tests CN against well-known identities.
    """

    from lielab.domain import CN

    x = CN([1,0,0])
    y = CN([0,1,0])
    z = CN([0,0,1])

    # Group multiplication
    x*y
    assert_domain(x*y, (y.inverse()*x.inverse()).inverse())
    assert_domain(x*y, y*x) # Group is abelian

    # Group inverse
    x.inverse()
    assert_domain(x.inverse(), CN(-x._data))


def test_GL():
    """
    Tests GL against well-known identities.
    """

    from lielab.domain import gl
    from lielab.functions import exp

    x = exp(gl.basis(0,2))
    y = exp(gl.basis(1,2))
    z = exp(gl.basis(2,2))

    # Group multiplication
    x*y
    assert_domain(x*y, (y.inverse()*x.inverse()).inverse())

    # Group inverse
    x.inverse()


def test_GLC():
    """
    Tests GLC against well-known identities.
    """

    from lielab.domain import glc
    from lielab.functions import exp

    x = exp(glc.basis(0,2))
    y = exp(glc.basis(1,2))
    z = exp(glc.basis(2,2))

    # Group multiplication
    x*y
    assert_domain(x*y, (y.inverse()*x.inverse()).inverse())

    # Group inverse
    x.inverse()


def test_RN():
    """
    Tests RN against well-known identities.
    """

    from lielab.domain import RN

    x = RN([1,0,0])
    y = RN([0,1,0])
    z = RN([0,0,1])

    # Group multiplication
    x*y
    assert_domain(x*y, (y.inverse()*x.inverse()).inverse())
    assert_domain(x*y, y*x) # Group is abelian

    # Group inverse
    x.inverse()
    assert_domain(x.inverse(), RN(-x._data))


def test_SO():
    """
    Tests SO against well-known identities.
    """

    from lielab.domain import so
    from lielab.functions import exp

    x = exp(so.basis(0,3))
    y = exp(so.basis(1,3))
    z = exp(so.basis(2,3))

    # Group product
    x*y
    assert_domain(x*y, (y.inverse()*x.inverse()).inverse())

    # Group inverse
    x.inverse()


def test_SP():
    """
    Tests SP against well-known identities.
    """

    from lielab.domain import sp
    from lielab.functions import exp

    x = exp(sp.basis(0,2))
    y = exp(sp.basis(1,2))
    z = exp(sp.basis(2,2))

    # Group product
    x*y
    assert_domain(x*y, (y.inverse()*x.inverse()).inverse())

    # Group inverse
    x.inverse()


def test_SU():
    """
    Tests SU against well-known identities.
    """

    from lielab.domain import SU

    u = SU(2)
    v = SU(2)
    w = SU(2)

    # Group product
    u*v
    assert_domain(u*v, (v.inverse()*u.inverse()).inverse())

    # Group inverse
    u.inverse()


def test_Quaternion():
    """
    Tests quaternions against well-known identities.
    """

    from lielab.domain import SU

    qm1 = SU.Quaternion(-1, 0, 0, 0)
    qi = SU.Quaternion(0, 1, 0, 0)
    qj = SU.Quaternion(0, 0, 1, 0)
    qk = SU.Quaternion(0, 0, 0, 1)

    # Hamilton's identities
    # i^2 = j^2 = k^2 = -1
    assert_matrix(qi*qi, qm1)
    assert_matrix(qj*qj, qm1)
    assert_matrix(qk*qk, qm1)

    # ij = -ji = -k
    assert_matrix(qi*qj, (qj*qi).inverse())
    assert_matrix(qi*qj, qk.inverse())

    # jk = -kj = -i
    assert_matrix(qj*qk, (qk*qj).inverse())
    assert_matrix(qj*qk, qi.inverse())

    # ki = -ik = -j
    assert_matrix(qk*qi, (qi*qk).inverse())
    assert_matrix(qk*qi, qj.inverse())
