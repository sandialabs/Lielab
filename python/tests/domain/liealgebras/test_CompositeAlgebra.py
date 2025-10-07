import lielab
import numpy as np

def complex(a,b):
    return a + b*1j

ycn1 = lielab.domain.cn.from_vector([1.0, 2.0, 3.0, 4.0])
yglc1 = lielab.domain.glc.from_vector([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
yglr1 = lielab.domain.glr.from_vector([13.0, 14.0, 15.0, 16.0])
yrn1 = lielab.domain.rn.from_vector([17.0, 18.0])
yse1 = lielab.domain.se.from_vector([19.0, 20.0, 21.0])
yso1 = lielab.domain.so.from_vector([22.0, 23.0, 24.0])
ysp1 = lielab.domain.sp.from_vector([25.0, 26.0, 27.0])
ysu1 = lielab.domain.su.from_vector([28.0, 29.0, 30.0])

def test_CompositeAlgebra_to_string():
    from lielab.domain import CompositeAlgebra

    y1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])
    elements = []
    elements.append(ycn1)
    elements.append(yglc1)
    elements.append(yglr1)
    elements.append(yrn1)
    elements.append(yse1)
    elements.append(yso1)
    elements.append(ysp1)
    elements.append(ysu1)
    y2 = CompositeAlgebra(elements)

def test_CompositeAlgebra_main_initializer():
    from lielab.domain import CompositeAlgebra

    xblank = CompositeAlgebra()
    assert (xblank.get_dimension() == 0)

    x0 = CompositeAlgebra(0)
    assert (x0.get_dimension() == 0)
    x1 = CompositeAlgebra(1)
    assert (x1.get_dimension() == 2)
    x10 = CompositeAlgebra(10)
    assert (x10.get_dimension() == 200)

def CompositeAlgebra_list_initializer():
    from lielab.domain import CompositeAlgebra

    y1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1]);
    elements = []
    elements.append(ycn1)
    elements.append(yglc1)
    elements.append(yglr1)
    elements.append(yrn1)
    elements.append(yse1)
    elements.append(yso1)
    elements.append(ysp1)
    elements.append(ysu1)
    y2 = CompositeAlgebra(elements)

def test_CompositeAlgebra_basis_initializer():
    from lielab.domain import CompositeAlgebra

    xm10 = CompositeAlgebra.basis(-1, 0)
    assert (xm10.get_dimension() == 0)
    xm10bar = xm10.get_vector()
    assert (xm10bar.size == 0)

    x00 = CompositeAlgebra.basis(0, 0)
    assert (x00.get_dimension() == 0)
    x00bar = x00.get_vector()
    assert (x00bar.size == 0)

    x10 = CompositeAlgebra.basis(1, 0)
    assert (x10.get_dimension() == 0)
    x10bar = x10.get_vector()
    assert (x10bar.size == 0)

    x01 = CompositeAlgebra.basis(0, 1)
    assert (x01.get_dimension() == 2)
    x01bar = x01.get_vector()
    assert (x01bar.size == 2)
    assert (x01bar[0] == 1.0)
    assert (x01bar[1] == 0.0)

    x11 = CompositeAlgebra.basis(1, 1)
    assert (x11.get_dimension() == 2)
    x11bar = x11.get_vector()
    assert (x11bar.size == 2)
    assert (x11bar[0] == 0.0)
    assert (x11bar[1] == 1.0)

    x21 = CompositeAlgebra.basis(2, 1)
    assert (x21.get_dimension() == 2)
    x21bar = x21.get_vector()
    assert (x21bar.size == 2)
    assert (x21bar[0] == 0.0)
    assert (x21bar[1] == 0.0)

    x02 = CompositeAlgebra.basis(0, 2)
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

def test_CompositeAlgebra_from_shape_initializer():
    from lielab.domain import CompositeAlgebra

    x0 = CompositeAlgebra.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = CompositeAlgebra.from_shape(1)
    assert (x1.get_dimension() == 2)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = CompositeAlgebra.from_shape(2)
    assert (x2.get_dimension() == 8)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_CompositeAlgebra_get_dimension():
    from lielab.domain import CompositeAlgebra

    zero = CompositeAlgebra(0)
    one = CompositeAlgebra(1)
    two = CompositeAlgebra(2)
    three = CompositeAlgebra(3)
    four = CompositeAlgebra(4)
    five = CompositeAlgebra(5)
    six = CompositeAlgebra(6)
    seven = CompositeAlgebra(7)
    eight = CompositeAlgebra(8)

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

def test_CompositeAlgebra_set_get_vector():
    """
    * Tests the set/get_vector operation.
    """

    from lielab.domain import CompositeAlgebra

    y1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])
    y1.set_vector([30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0,
                   20.0, 19.0, 18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0,
                   10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
    
    y1bar = y1.get_vector()
    assert (y1bar.size == 30)
    assert (y1bar[0] == 30.0)
    assert (y1bar[1] == 29.0)
    assert (y1bar[2] == 28.0)
    assert (y1bar[3] == 27.0)
    assert (y1bar[4] == 26.0)
    assert (y1bar[5] == 25.0)
    assert (y1bar[6] == 24.0)
    assert (y1bar[7] == 23.0)
    assert (y1bar[8] == 22.0)
    assert (y1bar[9] == 21.0)
    assert (y1bar[10] == 20.0)
    assert (y1bar[11] == 19.0)
    assert (y1bar[12] == 18.0)
    assert (y1bar[13] == 17.0)
    assert (y1bar[14] == 16.0)
    assert (y1bar[15] == 15.0)
    assert (y1bar[16] == 14.0)
    assert (y1bar[17] == 13.0)
    assert (y1bar[18] == 12.0)
    assert (y1bar[19] == 11.0)
    assert (y1bar[20] == 10.0)
    assert (y1bar[21] == 9.0)
    assert (y1bar[22] == 8.0)
    assert (y1bar[23] == 7.0)
    assert (y1bar[24] == 6.0)
    assert (y1bar[25] == 5.0)
    assert (y1bar[26] == 4.0)
    assert (y1bar[27] == 3.0)
    assert (y1bar[28] == 2.0)
    assert (y1bar[29] == 1.0)

def test_CompositeAlgebra_get_matrix():
    """
    * Tests the get_matrix operation.
    """

    from lielab.domain import CompositeAlgebra

    y1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])
    y1.set_vector([30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0,
                   20.0, 19.0, 18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0,
                   10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])

    y1hat = y1.get_matrix()

    assert (y1hat.shape[0] == 20)
    assert (y1hat.shape[1] == 20)
    
    # TODO: Check the 0's
    # TODO: Check each submatrix


def test_CompositeAlgebra_operator_parenthesis():
    from lielab.domain import CompositeAlgebra, cn, glc

    y1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])

    # In bounds cn component
    assert (y1(0) == 1.0)
    assert (y1(1) == 2.0)
    assert (y1(2) == 3.0)
    assert (y1(3) == 4.0)
    assert (y1(-30) == 1.0)
    assert (y1(-29) == 2.0)
    assert (y1(-28) == 3.0)
    assert (y1(-27) == 4.0)

    # In bounds glc component
    assert (y1(4) == 5.0)
    assert (y1(5) == 6.0)
    assert (y1(6) == 7.0)
    assert (y1(7) == 8.0)
    assert (y1(8) == 9.0)
    assert (y1(9) == 10.0)
    assert (y1(10) == 11.0)
    assert (y1(11) == 12.0)
    assert (y1(-26) == 5.0)
    assert (y1(-25) == 6.0)
    assert (y1(-24) == 7.0)
    assert (y1(-23) == 8.0)
    assert (y1(-22) == 9.0)
    assert (y1(-21) == 10.0)
    assert (y1(-20) == 11.0)
    assert (y1(-19) == 12.0)

    # In bounds glr component
    assert (y1(12) == 13.0)
    assert (y1(13) == 14.0)
    assert (y1(14) == 15.0)
    assert (y1(15) == 16.0)
    assert (y1(-18) == 13.0)
    assert (y1(-17) == 14.0)
    assert (y1(-16) == 15.0)
    assert (y1(-15) == 16.0)

    # In bounds rn component
    assert (y1(16) == 17.0)
    assert (y1(17) == 18.0)
    assert (y1(-14) == 17.0)
    assert (y1(-13) == 18.0)

    # In bounds se component
    assert (y1(18) == 19.0)
    assert (y1(19) == 20.0)
    assert (y1(20) == 21.0)
    assert (y1(-12) == 19.0)
    assert (y1(-11) == 20.0)
    assert (y1(-10) == 21.0)

    # In bounds so component
    assert (y1(21) == 22.0)
    assert (y1(22) == 23.0)
    assert (y1(23) == 24.0)
    assert (y1(-9) == 22.0)
    assert (y1(-8) == 23.0)
    assert (y1(-7) == 24.0)

    # In bounds sp component
    assert (y1(24) == 25.0)
    assert (y1(25) == 26.0)
    assert (y1(26) == 27.0)
    assert (y1(-6) == 25.0)
    assert (y1(-5) == 26.0)
    assert (y1(-4) == 27.0)

    # In bounds su component
    assert (y1(27) == 28.0)
    assert (y1(28) == 29.0)
    assert (y1(29) == 30.0)
    assert (y1(-3) == 28.0)
    assert (y1(-2) == 29.0)
    assert (y1(-1) == 30.0)

    # Out of bounds
    assert (np.isnan(y1(-31)))
    assert (np.isnan(y1(30)))

    # In bounds cn component
    assert (y1(0, 0) == complex(0.0, 0.0))
    assert (y1(0, 1) == complex(0.0, 0.0))
    assert (y1(0, 2) == complex(1.0, 2.0))
    assert (y1(1, 0) == complex(0.0, 0.0))
    assert (y1(1, 1) == complex(0.0, 0.0))
    assert (y1(1, 2) == complex(3.0, 4.0))
    assert (y1(2, 0) == complex(0.0, 0.0))
    assert (y1(2, 1) == complex(0.0, 0.0))
    assert (y1(2, 2) == complex(0.0, 0.0))
    assert (y1(-20, -20) == complex(0.0, 0.0))
    assert (y1(-20, -19) == complex(0.0, 0.0))
    assert (y1(-20, -18) == complex(1.0, 2.0))
    assert (y1(-19, -20) == complex(0.0, 0.0))
    assert (y1(-19, -19) == complex(0.0, 0.0))
    assert (y1(-19, -18) == complex(3.0, 4.0))
    assert (y1(-18, -20) == complex(0.0, 0.0))
    assert (y1(-19, -19) == complex(0.0, 0.0))
    assert (y1(-18, -18) == complex(0.0, 0.0))

    # In bounds glc component
    assert (y1(3, 3) == complex(5.0, 6.0))
    assert (y1(3, 4) == complex(7.0, 8.0))
    assert (y1(4, 3) == complex(9.0, 10.0))
    assert (y1(4, 4) == complex(11.0, 12.0))
    assert (y1(-17, -17) == complex(5.0, 6.0))
    assert (y1(-17, -16) == complex(7.0, 8.0))
    assert (y1(-16, -17) == complex(9.0, 10.0))
    assert (y1(-16, -16) == complex(11.0, 12.0))

    # In bounds glr component
    assert (y1(5, 5) == complex(13.0, 0.0))
    assert (y1(5, 6) == complex(14.0, 0.0))
    assert (y1(6, 5) == complex(15.0, 0.0))
    assert (y1(6, 6) == complex(16.0, 0.0))
    assert (y1(-15, -15) == complex(13.0, 0.0))
    assert (y1(-15, -14) == complex(14.0, 0.0))
    assert (y1(-14, -15) == complex(15.0, 0.0))
    assert (y1(-14, -14) == complex(16.0, 0.0))

    # In bounds rn component
    assert (y1(7, 7) == complex(0.0, 0.0))
    assert (y1(7, 8) == complex(0.0, 0.0))
    assert (y1(7, 9) == complex(17.0, 0.0))
    assert (y1(8, 7) == complex(0.0, 0.0))
    assert (y1(8, 8) == complex(0.0, 0.0))
    assert (y1(8, 9) == complex(18.0, 0.0))
    assert (y1(9, 7) == complex(0.0, 0.0))
    assert (y1(9, 8) == complex(0.0, 0.0))
    assert (y1(9, 9) == complex(0.0, 0.0))
    assert (y1(-13, -13) == complex(0.0, 0.0))
    assert (y1(-13, -12) == complex(0.0, 0.0))
    assert (y1(-13, -11) == complex(17.0, 0.0))
    assert (y1(-12, -13) == complex(0.0, 0.0))
    assert (y1(-12, -12) == complex(0.0, 0.0))
    assert (y1(-12, -11) == complex(18.0, 0.0))
    assert (y1(-11, -13) == complex(0.0, 0.0))
    assert (y1(-11, -12) == complex(0.0, 0.0))
    assert (y1(-11, -11) == complex(0.0, 0.0))

    # In bounds se component
    assert (y1(10, 10) == complex(0.0, 0.0))
    assert (y1(10, 11) == complex(-21.0, 0.0))
    assert (y1(10, 12) == complex(19.0, 0.0))
    assert (y1(11, 10) == complex(21.0, 0.0))
    assert (y1(11, 11) == complex(0.0, 0.0))
    assert (y1(11, 12) == complex(20.0, 0.0))
    assert (y1(12, 10) == complex(0.0, 0.0))
    assert (y1(12, 11) == complex(0.0, 0.0))
    assert (y1(12, 12) == complex(0.0, 0.0))
    assert (y1(-10, -10) == complex(0.0, 0.0))
    assert (y1(-10, -9) == complex(-21.0, 0.0))
    assert (y1(-10, -8) == complex(19.0, 0.0))
    assert (y1(-9, -10) == complex(21.0, 0.0))
    assert (y1(-9, -9) == complex(0.0, 0.0))
    assert (y1(-9, -8) == complex(20.0, 0.0))
    assert (y1(-8, -10) == complex(0.0, 0.0))
    assert (y1(-8, -9) == complex(0.0, 0.0))
    assert (y1(-8, -8) == complex(0.0, 0.0))

    # In bounds so component
    assert (y1(13, 13) == complex(0.0, 0.0))
    assert (y1(13, 14) == complex(-24.0, 0.0))
    assert (y1(13, 15) == complex(23.0, 0.0))
    assert (y1(14, 13) == complex(24.0, 0.0))
    assert (y1(14, 14) == complex(0.0, 0.0))
    assert (y1(14, 15) == complex(-22.0, 0.0))
    assert (y1(15, 13) == complex(-23.0, 0.0))
    assert (y1(15, 14) == complex(22.0, 0.0))
    assert (y1(15, 15) == complex(0.0, 0.0))
    assert (y1(-7, -7) == complex(0.0, 0.0))
    assert (y1(-7, -6) == complex(-24.0, 0.0))
    assert (y1(-7, -5) == complex(23.0, 0.0))
    assert (y1(-6, -7) == complex(24.0, 0.0))
    assert (y1(-6, -6) == complex(0.0, 0.0))
    assert (y1(-6, -5) == complex(-22.0, 0.0))
    assert (y1(-5, -7) == complex(-23.0, 0.0))
    assert (y1(-5, -6) == complex(22.0, 0.0))
    assert (y1(-5, -5) == complex(0.0, 0.0))

    # In bounds sp component
    assert (y1(16, 16) == complex(25.0, 0.0))
    assert (y1(16, 17) == complex(26.0, 0.0))
    assert (y1(17, 16) == complex(27.0, 0.0))
    assert (y1(17, 17) == complex(-25.0, 0.0))
    assert (y1(-4, -4) == complex(25.0, 0.0))
    assert (y1(-4, -3) == complex(26.0, 0.0))
    assert (y1(-3, -4) == complex(27.0, 0.0))
    assert (y1(-3, -3) == complex(-25.0, 0.0))

    # In bounds su component
    assert (y1(18, 18) == complex(0.0, 30.0))
    assert (y1(18, 19) == complex(-29.0, 28.0))
    assert (y1(19, 18) == complex(29.0, 28.0))
    assert (y1(19, 19) == complex(0.0, -30.0))
    assert (y1(-2, -2) == complex(0.0, 30.0))
    assert (y1(-2, -1) == complex(-29.0, 28.0))
    assert (y1(-1, -2) == complex(29.0, 28.0))
    assert (y1(-1, -1) == complex(0.0, -30.0))

    # In bounds sparse components
    assert (y1(0, 3) == complex(0.0, 0.0))
    assert (y1(0, 4) == complex(0.0, 0.0))
    assert (y1(1, 3) == complex(0.0, 0.0))
    assert (y1(1, 4) == complex(0.0, 0.0))
    assert (y1(2, 3) == complex(0.0, 0.0))
    assert (y1(2, 4) == complex(0.0, 0.0))
    assert (y1(3, 0) == complex(0.0, 0.0))
    assert (y1(4, 0) == complex(0.0, 0.0))
    assert (y1(3, 1) == complex(0.0, 0.0))
    assert (y1(4, 1) == complex(0.0, 0.0))
    assert (y1(3, 2) == complex(0.0, 0.0))
    assert (y1(4, 2) == complex(0.0, 0.0))
    # Probably don't need to check the rest

    # Out of bounds
    assert (np.isnan(np.real(y1(0, -21))))
    assert (np.isnan(np.imag(y1(0, -21))))
    assert (np.isnan(np.real(y1(-21, 0))))
    assert (np.isnan(np.imag(y1(-21, 0))))
    assert (np.isnan(np.real(y1(-21, -21))))
    assert (np.isnan(np.imag(y1(-21, -21))))
    assert (np.isnan(np.real(y1(0, 20))))
    assert (np.isnan(np.imag(y1(0, 20))))
    assert (np.isnan(np.real(y1(20, 0))))
    assert (np.isnan(np.imag(y1(20, 0))))
    assert (np.isnan(np.real(y1(20, 20))))
    assert (np.isnan(np.imag(y1(20, 20))))

# TODO: mathops int here

def test_CompositeAlgebra_math_ops_double():
    from lielab.domain import CompositeAlgebra

    x1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])

    x1_lm_2 = 2.0*x1
    x1_lm_2bar = x1_lm_2.get_vector()
    assert (x1_lm_2bar.size == 30)
    assert (x1_lm_2bar[0] == 2.0)
    assert (x1_lm_2bar[1] == 4.0)
    assert (x1_lm_2bar[2] == 6.0)
    assert (x1_lm_2bar[3] == 8.0)
    assert (x1_lm_2bar[4] == 10.0)
    assert (x1_lm_2bar[5] == 12.0)
    assert (x1_lm_2bar[6] == 14.0)
    assert (x1_lm_2bar[7] == 16.0)
    assert (x1_lm_2bar[8] == 18.0)
    assert (x1_lm_2bar[9] == 20.0)
    assert (x1_lm_2bar[10] == 22.0)
    assert (x1_lm_2bar[11] == 24.0)
    assert (x1_lm_2bar[12] == 26.0)
    assert (x1_lm_2bar[13] == 28.0)
    assert (x1_lm_2bar[14] == 30.0)
    assert (x1_lm_2bar[15] == 32.0)
    assert (x1_lm_2bar[16] == 34.0)
    assert (x1_lm_2bar[17] == 36.0)
    assert (x1_lm_2bar[18] == 38.0)
    assert (x1_lm_2bar[19] == 40.0)
    assert (x1_lm_2bar[20] == 42.0)
    assert (x1_lm_2bar[21] == 44.0)
    assert (x1_lm_2bar[22] == 46.0)
    assert (x1_lm_2bar[23] == 48.0)
    assert (x1_lm_2bar[24] == 50.0)
    assert (x1_lm_2bar[25] == 52.0)
    assert (x1_lm_2bar[26] == 54.0)
    assert (x1_lm_2bar[27] == 56.0)
    assert (x1_lm_2bar[28] == 58.0)
    assert (x1_lm_2bar[29] == 60.0)

    x1_rm_2 = x1*2.0
    x1_rm_2bar = x1_rm_2.get_vector()
    assert (x1_rm_2bar.size == 30)
    assert (x1_rm_2bar[0] == 2.0)
    assert (x1_rm_2bar[1] == 4.0)
    assert (x1_rm_2bar[2] == 6.0)
    assert (x1_rm_2bar[3] == 8.0)
    assert (x1_rm_2bar[4] == 10.0)
    assert (x1_rm_2bar[5] == 12.0)
    assert (x1_rm_2bar[6] == 14.0)
    assert (x1_rm_2bar[7] == 16.0)
    assert (x1_rm_2bar[8] == 18.0)
    assert (x1_rm_2bar[9] == 20.0)
    assert (x1_rm_2bar[10] == 22.0)
    assert (x1_rm_2bar[11] == 24.0)
    assert (x1_rm_2bar[12] == 26.0)
    assert (x1_rm_2bar[13] == 28.0)
    assert (x1_rm_2bar[14] == 30.0)
    assert (x1_rm_2bar[15] == 32.0)
    assert (x1_rm_2bar[16] == 34.0)
    assert (x1_rm_2bar[17] == 36.0)
    assert (x1_rm_2bar[18] == 38.0)
    assert (x1_rm_2bar[19] == 40.0)
    assert (x1_rm_2bar[20] == 42.0)
    assert (x1_rm_2bar[21] == 44.0)
    assert (x1_rm_2bar[22] == 46.0)
    assert (x1_rm_2bar[23] == 48.0)
    assert (x1_rm_2bar[24] == 50.0)
    assert (x1_rm_2bar[25] == 52.0)
    assert (x1_rm_2bar[26] == 54.0)
    assert (x1_rm_2bar[27] == 56.0)
    assert (x1_rm_2bar[28] == 58.0)
    assert (x1_rm_2bar[29] == 60.0)

    x1 *= 2.0
    x1_im_2bar = x1.get_vector()
    assert (x1_im_2bar.size == 30)
    assert (x1_im_2bar[0] == 2.0)
    assert (x1_im_2bar[1] == 4.0)
    assert (x1_im_2bar[2] == 6.0)
    assert (x1_im_2bar[3] == 8.0)
    assert (x1_im_2bar[4] == 10.0)
    assert (x1_im_2bar[5] == 12.0)
    assert (x1_im_2bar[6] == 14.0)
    assert (x1_im_2bar[7] == 16.0)
    assert (x1_im_2bar[8] == 18.0)
    assert (x1_im_2bar[9] == 20.0)
    assert (x1_im_2bar[10] == 22.0)
    assert (x1_im_2bar[11] == 24.0)
    assert (x1_im_2bar[12] == 26.0)
    assert (x1_im_2bar[13] == 28.0)
    assert (x1_im_2bar[14] == 30.0)
    assert (x1_im_2bar[15] == 32.0)
    assert (x1_im_2bar[16] == 34.0)
    assert (x1_im_2bar[17] == 36.0)
    assert (x1_im_2bar[18] == 38.0)
    assert (x1_im_2bar[19] == 40.0)
    assert (x1_im_2bar[20] == 42.0)
    assert (x1_im_2bar[21] == 44.0)
    assert (x1_im_2bar[22] == 46.0)
    assert (x1_im_2bar[23] == 48.0)
    assert (x1_im_2bar[24] == 50.0)
    assert (x1_im_2bar[25] == 52.0)
    assert (x1_im_2bar[26] == 54.0)
    assert (x1_im_2bar[27] == 56.0)
    assert (x1_im_2bar[28] == 58.0)
    assert (x1_im_2bar[29] == 60.0)

    x1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])

    x1_d_2 = x1/2.0
    x1_d_2bar = x1_d_2.get_vector()
    assert (x1_d_2bar.size == 30)
    assert (x1_d_2bar[0] == 0.5)
    assert (x1_d_2bar[1] == 1.0)
    assert (x1_d_2bar[2] == 1.5)
    assert (x1_d_2bar[3] == 2.0)
    assert (x1_d_2bar[4] == 2.5)
    assert (x1_d_2bar[5] == 3.0)
    assert (x1_d_2bar[6] == 3.5)
    assert (x1_d_2bar[7] == 4.0)
    assert (x1_d_2bar[8] == 4.5)
    assert (x1_d_2bar[9] == 5.0)
    assert (x1_d_2bar[10] == 5.5)
    assert (x1_d_2bar[11] == 6.0)
    assert (x1_d_2bar[12] == 6.5)
    assert (x1_d_2bar[13] == 7.0)
    assert (x1_d_2bar[14] == 7.5)
    assert (x1_d_2bar[15] == 8.0)
    assert (x1_d_2bar[16] == 8.5)
    assert (x1_d_2bar[17] == 9.0)
    assert (x1_d_2bar[18] == 9.5)
    assert (x1_d_2bar[19] == 10.0)
    assert (x1_d_2bar[20] == 10.5)
    assert (x1_d_2bar[21] == 11.0)
    assert (x1_d_2bar[22] == 11.5)
    assert (x1_d_2bar[23] == 12.0)
    assert (x1_d_2bar[24] == 12.5)
    assert (x1_d_2bar[25] == 13.0)
    assert (x1_d_2bar[26] == 13.5)
    assert (x1_d_2bar[27] == 14.0)
    assert (x1_d_2bar[28] == 14.5)
    assert (x1_d_2bar[29] == 15.0)

    x1 /= 2.0
    x1_id_2bar = x1.get_vector()
    assert (x1_id_2bar.size == 30)
    assert (x1_id_2bar[0] == 0.5)
    assert (x1_id_2bar[1] == 1.0)
    assert (x1_id_2bar[2] == 1.5)
    assert (x1_id_2bar[3] == 2.0)
    assert (x1_id_2bar[4] == 2.5)
    assert (x1_id_2bar[5] == 3.0)
    assert (x1_id_2bar[6] == 3.5)
    assert (x1_id_2bar[7] == 4.0)
    assert (x1_id_2bar[8] == 4.5)
    assert (x1_id_2bar[9] == 5.0)
    assert (x1_id_2bar[10] == 5.5)
    assert (x1_id_2bar[11] == 6.0)
    assert (x1_id_2bar[12] == 6.5)
    assert (x1_id_2bar[13] == 7.0)
    assert (x1_id_2bar[14] == 7.5)
    assert (x1_id_2bar[15] == 8.0)
    assert (x1_id_2bar[16] == 8.5)
    assert (x1_id_2bar[17] == 9.0)
    assert (x1_id_2bar[18] == 9.5)
    assert (x1_id_2bar[19] == 10.0)
    assert (x1_id_2bar[20] == 10.5)
    assert (x1_id_2bar[21] == 11.0)
    assert (x1_id_2bar[22] == 11.5)
    assert (x1_id_2bar[23] == 12.0)
    assert (x1_id_2bar[24] == 12.5)
    assert (x1_id_2bar[25] == 13.0)
    assert (x1_id_2bar[26] == 13.5)
    assert (x1_id_2bar[27] == 14.0)
    assert (x1_id_2bar[28] == 14.5)
    assert (x1_id_2bar[29] == 15.0)

def test_CompositeAlgebra_math_ops_CompositeAlgebra():
    from lielab.domain import CompositeAlgebra

    x1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])
    x2 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])

    x1_add_x2 = x1 + x2
    x1_add_x2bar = x1_add_x2.get_vector()
    assert (x1_add_x2bar.size == 30)
    assert (x1_add_x2bar[0] == 2.0)
    assert (x1_add_x2bar[1] == 4.0)
    assert (x1_add_x2bar[2] == 6.0)
    assert (x1_add_x2bar[3] == 8.0)
    assert (x1_add_x2bar[4] == 10.0)
    assert (x1_add_x2bar[5] == 12.0)
    assert (x1_add_x2bar[6] == 14.0)
    assert (x1_add_x2bar[7] == 16.0)
    assert (x1_add_x2bar[8] == 18.0)
    assert (x1_add_x2bar[9] == 20.0)
    assert (x1_add_x2bar[10] == 22.0)
    assert (x1_add_x2bar[11] == 24.0)
    assert (x1_add_x2bar[12] == 26.0)
    assert (x1_add_x2bar[13] == 28.0)
    assert (x1_add_x2bar[14] == 30.0)
    assert (x1_add_x2bar[15] == 32.0)
    assert (x1_add_x2bar[16] == 34.0)
    assert (x1_add_x2bar[17] == 36.0)
    assert (x1_add_x2bar[18] == 38.0)
    assert (x1_add_x2bar[19] == 40.0)
    assert (x1_add_x2bar[20] == 42.0)
    assert (x1_add_x2bar[21] == 44.0)
    assert (x1_add_x2bar[22] == 46.0)
    assert (x1_add_x2bar[23] == 48.0)
    assert (x1_add_x2bar[24] == 50.0)
    assert (x1_add_x2bar[25] == 52.0)
    assert (x1_add_x2bar[26] == 54.0)
    assert (x1_add_x2bar[27] == 56.0)
    assert (x1_add_x2bar[28] == 58.0)
    assert (x1_add_x2bar[29] == 60.0)

    x1 += x2
    x1_iadd_x2bar = x1.get_vector()
    assert (x1_iadd_x2bar.size == 30)
    assert (x1_iadd_x2bar[0] == 2.0)
    assert (x1_iadd_x2bar[1] == 4.0)
    assert (x1_iadd_x2bar[2] == 6.0)
    assert (x1_iadd_x2bar[3] == 8.0)
    assert (x1_iadd_x2bar[4] == 10.0)
    assert (x1_iadd_x2bar[5] == 12.0)
    assert (x1_iadd_x2bar[6] == 14.0)
    assert (x1_iadd_x2bar[7] == 16.0)
    assert (x1_iadd_x2bar[8] == 18.0)
    assert (x1_iadd_x2bar[9] == 20.0)
    assert (x1_iadd_x2bar[10] == 22.0)
    assert (x1_iadd_x2bar[11] == 24.0)
    assert (x1_iadd_x2bar[12] == 26.0)
    assert (x1_iadd_x2bar[13] == 28.0)
    assert (x1_iadd_x2bar[14] == 30.0)
    assert (x1_iadd_x2bar[15] == 32.0)
    assert (x1_iadd_x2bar[16] == 34.0)
    assert (x1_iadd_x2bar[17] == 36.0)
    assert (x1_iadd_x2bar[18] == 38.0)
    assert (x1_iadd_x2bar[19] == 40.0)
    assert (x1_iadd_x2bar[20] == 42.0)
    assert (x1_iadd_x2bar[21] == 44.0)
    assert (x1_iadd_x2bar[22] == 46.0)
    assert (x1_iadd_x2bar[23] == 48.0)
    assert (x1_iadd_x2bar[24] == 50.0)
    assert (x1_iadd_x2bar[25] == 52.0)
    assert (x1_iadd_x2bar[26] == 54.0)
    assert (x1_iadd_x2bar[27] == 56.0)
    assert (x1_iadd_x2bar[28] == 58.0)
    assert (x1_iadd_x2bar[29] == 60.0)

    x1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])

    x1_sub_x2 = x1 - x2
    x1_sub_x2bar = x1_sub_x2.get_vector()
    assert (x1_sub_x2bar.size == 30)
    assert (x1_sub_x2bar[0] == 0.0)
    assert (x1_sub_x2bar[1] == 0.0)
    assert (x1_sub_x2bar[2] == 0.0)
    assert (x1_sub_x2bar[3] == 0.0)
    assert (x1_sub_x2bar[4] == 0.0)
    assert (x1_sub_x2bar[5] == 0.0)
    assert (x1_sub_x2bar[6] == 0.0)
    assert (x1_sub_x2bar[7] == 0.0)
    assert (x1_sub_x2bar[8] == 0.0)
    assert (x1_sub_x2bar[9] == 0.0)
    assert (x1_sub_x2bar[10] == 0.0)
    assert (x1_sub_x2bar[11] == 0.0)
    assert (x1_sub_x2bar[12] == 0.0)
    assert (x1_sub_x2bar[13] == 0.0)
    assert (x1_sub_x2bar[14] == 0.0)
    assert (x1_sub_x2bar[15] == 0.0)
    assert (x1_sub_x2bar[16] == 0.0)
    assert (x1_sub_x2bar[17] == 0.0)
    assert (x1_sub_x2bar[18] == 0.0)
    assert (x1_sub_x2bar[19] == 0.0)
    assert (x1_sub_x2bar[20] == 0.0)
    assert (x1_sub_x2bar[21] == 0.0)
    assert (x1_sub_x2bar[22] == 0.0)
    assert (x1_sub_x2bar[23] == 0.0)
    assert (x1_sub_x2bar[24] == 0.0)
    assert (x1_sub_x2bar[25] == 0.0)
    assert (x1_sub_x2bar[26] == 0.0)
    assert (x1_sub_x2bar[27] == 0.0)
    assert (x1_sub_x2bar[28] == 0.0)
    assert (x1_sub_x2bar[29] == 0.0)

    x1 -= x2
    x1_isub_x2bar = x1.get_vector()
    assert (x1_isub_x2bar.size == 30)
    assert (x1_isub_x2bar[0] == 0.0)
    assert (x1_isub_x2bar[1] == 0.0)
    assert (x1_isub_x2bar[2] == 0.0)
    assert (x1_isub_x2bar[3] == 0.0)
    assert (x1_isub_x2bar[4] == 0.0)
    assert (x1_isub_x2bar[5] == 0.0)
    assert (x1_isub_x2bar[6] == 0.0)
    assert (x1_isub_x2bar[7] == 0.0)
    assert (x1_isub_x2bar[8] == 0.0)
    assert (x1_isub_x2bar[9] == 0.0)
    assert (x1_isub_x2bar[10] == 0.0)
    assert (x1_isub_x2bar[11] == 0.0)
    assert (x1_isub_x2bar[12] == 0.0)
    assert (x1_isub_x2bar[13] == 0.0)
    assert (x1_isub_x2bar[14] == 0.0)
    assert (x1_isub_x2bar[15] == 0.0)
    assert (x1_isub_x2bar[16] == 0.0)
    assert (x1_isub_x2bar[17] == 0.0)
    assert (x1_isub_x2bar[18] == 0.0)
    assert (x1_isub_x2bar[19] == 0.0)
    assert (x1_isub_x2bar[20] == 0.0)
    assert (x1_isub_x2bar[21] == 0.0)
    assert (x1_isub_x2bar[22] == 0.0)
    assert (x1_isub_x2bar[23] == 0.0)
    assert (x1_isub_x2bar[24] == 0.0)
    assert (x1_isub_x2bar[25] == 0.0)
    assert (x1_isub_x2bar[26] == 0.0)
    assert (x1_isub_x2bar[27] == 0.0)
    assert (x1_isub_x2bar[28] == 0.0)
    assert (x1_isub_x2bar[29] == 0.0)

    x1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])

    x1_unary_sub = (-x1)
    x1_usub_bar = x1_unary_sub.get_vector()
    assert (x1_usub_bar.size == 30)
    assert (x1_usub_bar[0] == -1.0)
    assert (x1_usub_bar[1] == -2.0)
    assert (x1_usub_bar[2] == -3.0)
    assert (x1_usub_bar[3] == -4.0)
    assert (x1_usub_bar[4] == -5.0)
    assert (x1_usub_bar[5] == -6.0)
    assert (x1_usub_bar[6] == -7.0)
    assert (x1_usub_bar[7] == -8.0)
    assert (x1_usub_bar[8] == -9.0)
    assert (x1_usub_bar[9] == -10.0)
    assert (x1_usub_bar[10] == -11.0)
    assert (x1_usub_bar[11] == -12.0)
    assert (x1_usub_bar[12] == -13.0)
    assert (x1_usub_bar[13] == -14.0)
    assert (x1_usub_bar[14] == -15.0)
    assert (x1_usub_bar[15] == -16.0)
    assert (x1_usub_bar[16] == -17.0)
    assert (x1_usub_bar[17] == -18.0)
    assert (x1_usub_bar[18] == -19.0)
    assert (x1_usub_bar[19] == -20.0)
    assert (x1_usub_bar[20] == -21.0)
    assert (x1_usub_bar[21] == -22.0)
    assert (x1_usub_bar[22] == -23.0)
    assert (x1_usub_bar[23] == -24.0)
    assert (x1_usub_bar[24] == -25.0)
    assert (x1_usub_bar[25] == -26.0)
    assert (x1_usub_bar[26] == -27.0)
    assert (x1_usub_bar[27] == -28.0)
    assert (x1_usub_bar[28] == -29.0)
    assert (x1_usub_bar[29] == -30.0)

def test_CompositeAlgebra_operator_bracket():
    from lielab.domain import CompositeAlgebra

    x1 = CompositeAlgebra([ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])
    
    x10 = x1[0]
    x11 = x1[1]
    x12 = x1[2]
    x13 = x1[3]
    x14 = x1[4]
    x15 = x1[5]
    x16 = x1[6]
    x17 = x1[7]
    x1m8 = x1[-8]
    x1m7 = x1[-7]
    x1m6 = x1[-6]
    x1m5 = x1[-5]
    x1m4 = x1[-4]
    x1m3 = x1[-3]
    x1m2 = x1[-2]
    x1m1 = x1[-1]

    assert (x10.to_string() == "c^2")
    assert (x11.to_string() == "gl(2, C)")
    assert (x12.to_string() == "gl(2, R)")
    assert (x13.to_string() == "r^2")
    assert (x14.to_string() == "se(2)")
    assert (x15.to_string() == "so(3)")
    assert (x16.to_string() == "sp(2, R)")
    assert (x17.to_string() == "su(2)")
    assert (x1m8.to_string() == "c^2")
    assert (x1m7.to_string() == "gl(2, C)")
    assert (x1m6.to_string() == "gl(2, R)")
    assert (x1m5.to_string() == "r^2")
    assert (x1m4.to_string() == "se(2)")
    assert (x1m3.to_string() == "so(3)")
    assert (x1m2.to_string() == "sp(2, R)")
    assert (x1m1.to_string() == "su(2)")

    # Out of bounds
    x18 = x1[8]
    x1m9 = x1[-9]

    assert (x18.to_string() == "gl(0, C)")
    assert (x1m9.to_string() == "gl(0, C)")
