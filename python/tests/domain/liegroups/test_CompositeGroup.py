import numpy as np

def complex(a,b):
    return a + b*1j

def _make_cgroup():
    from lielab.domain import CompositeGroup
    import lielab

    yCN1 = lielab.domain.CN(2)
    yCN1.unserialize([1.0, 2.0, 3.0, 4.0])
    yGLC1 = lielab.domain.GLC(2)
    yGLC1.unserialize([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    yGLR1 = lielab.domain.GLR(2)
    yGLR1.unserialize([13.0, 14.0, 15.0, 16.0])
    yRN1 = lielab.domain.RN(2)
    yRN1.unserialize([17.0, 18.0])
    ySE1 = lielab.domain.SE(2)
    ySE1.unserialize([19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0])
    ySO1 = lielab.domain.SO(2)
    ySO1.unserialize([28.0, 29.0, 30.0, 31.0])
    ySP1 = lielab.domain.SP(2)
    ySP1.unserialize([32.0, 33.0, 34.0, 35.0])
    ySU1 = lielab.domain.SU(2)
    ySU1.unserialize([36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0])

    return CompositeGroup([yCN1, yGLC1, yGLR1, yRN1, ySE1, ySO1, ySP1, ySU1])

def test_CompositeGroup_to_string():
    from lielab.domain import CompositeGroup

    y1 = _make_cgroup()

    assert (y1.to_string() == "C^2 x GL(2, C) x GL(2, R) x R^2 x SE(2) x SO(2) x SP(2, R) x SU(2)")

def test_CompositeGroup_main_initializer():
    from lielab.domain import CompositeGroup

    xblank = CompositeGroup()
    assert (xblank.get_dimension() == 0)

    x0 = CompositeGroup(0)
    assert (x0.get_dimension() == 0)
    x1 = CompositeGroup(1)
    assert (x1.get_dimension() == 2)
    x10 = CompositeGroup(10)
    assert (x10.get_dimension() == 200)

def test_CompositeGroup_from_shape_initializer():
    from lielab.domain import CompositeGroup

    x0 = CompositeGroup.from_shape(0)
    assert (x0.get_dimension() == 0)
    x0hat = x0.get_matrix()
    assert (x0hat.shape[0] == 0)
    assert (x0hat.shape[1] == 0)

    x1 = CompositeGroup.from_shape(1)
    assert (x1.get_dimension() == 2)
    x1hat = x1.get_matrix()
    assert (x1hat.shape[0] == 1)
    assert (x1hat.shape[1] == 1)

    x2 = CompositeGroup.from_shape(2)
    assert (x2.get_dimension() == 8)
    x2hat = x2.get_matrix()
    assert (x2hat.shape[0] == 2)
    assert (x2hat.shape[1] == 2)

def test_CompositeGroup_get_dimension():
    from lielab.domain import CompositeGroup

    zero = CompositeGroup(0)
    one = CompositeGroup(1)
    two = CompositeGroup(2)
    three = CompositeGroup(3)
    four = CompositeGroup(4)
    five = CompositeGroup(5)
    six = CompositeGroup(6)
    seven = CompositeGroup(7)
    eight = CompositeGroup(8)

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

    y1 = _make_cgroup()
    
    dims = y1.get_dimensions()
    assert (len(dims) == 8)
    assert (dims[0] == 4)
    assert (dims[1] == 8)
    assert (dims[2] == 4)
    assert (dims[3] == 2)
    assert (dims[4] == 3)
    assert (dims[5] == 1)
    assert (dims[6] == 3)
    assert (dims[7] == 3)

    assert (y1.get_dimension() == 28)

def test_CompositeGroup_get_size():
    from lielab.domain import CompositeGroup

    zero = CompositeGroup(0)
    one = CompositeGroup(1)
    two = CompositeGroup(2)
    three = CompositeGroup(3)
    four = CompositeGroup(4)
    five = CompositeGroup(5)
    six = CompositeGroup(6)
    seven = CompositeGroup(7)
    eight = CompositeGroup(8)

    # Dimensions
    assert (zero.get_size() == 0)
    assert (one.get_size() == 2)
    assert (two.get_size() == 8)
    assert (three.get_size() == 18)
    assert (four.get_size() == 32)
    assert (five.get_size() == 50)
    assert (six.get_size() == 72)
    assert (seven.get_size() == 98)
    assert (eight.get_dimension() == 128)

    y1 = _make_cgroup()
    
    sizes = y1.get_sizes()
    assert (len(sizes) == 8)
    assert (sizes[0] == 4)
    assert (sizes[1] == 8)
    assert (sizes[2] == 4)
    assert (sizes[3] == 2)
    assert (sizes[4] == 9)
    assert (sizes[5] == 4)
    assert (sizes[6] == 4)
    assert (sizes[7] == 8)

    assert (y1.get_size() == 43)

def test_CompositeGroup_serialize_unserialize():
    """
    Tests the serialize/unserialize operation.
    """

    from lielab.domain import CompositeGroup

    y1 = _make_cgroup()
    y1.unserialize([43.0, 42.0, 41.0, 40.0, 39.0, 38.0, 37.0, 36.0, 35.0, 34.0, 33.0, 32.0, 31.0, 30.0,
                    29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0, 20.0, 19.0, 18.0, 17.0, 16.0,
                    15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
    
    y1bar = y1.serialize()
    assert (y1bar.size == 43)
    assert (y1bar[0] == 43.0)
    assert (y1bar[1] == 42.0)
    assert (y1bar[2] == 41.0)
    assert (y1bar[3] == 40.0)
    assert (y1bar[4] == 39.0)
    assert (y1bar[5] == 38.0)
    assert (y1bar[6] == 37.0)
    assert (y1bar[7] == 36.0)
    assert (y1bar[8] == 35.0)
    assert (y1bar[9] == 34.0)
    assert (y1bar[10] == 33.0)
    assert (y1bar[11] == 32.0)
    assert (y1bar[12] == 31.0)
    assert (y1bar[13] == 30.0)
    assert (y1bar[14] == 29.0)
    assert (y1bar[15] == 28.0)
    assert (y1bar[16] == 27.0)
    assert (y1bar[17] == 26.0)
    assert (y1bar[18] == 25.0)
    assert (y1bar[19] == 24.0)
    assert (y1bar[20] == 23.0)
    assert (y1bar[21] == 22.0)
    assert (y1bar[22] == 21.0)
    assert (y1bar[23] == 20.0)
    assert (y1bar[24] == 19.0)
    assert (y1bar[25] == 18.0)
    assert (y1bar[26] == 17.0)
    assert (y1bar[27] == 16.0)
    assert (y1bar[28] == 15.0)
    assert (y1bar[29] == 14.0)
    assert (y1bar[30] == 13.0)
    assert (y1bar[31] == 12.0)
    assert (y1bar[32] == 11.0)
    assert (y1bar[33] == 10.0)
    assert (y1bar[34] == 9.0)
    assert (y1bar[35] == 8.0)
    assert (y1bar[36] == 7.0)
    assert (y1bar[37] == 6.0)
    assert (y1bar[38] == 5.0)
    assert (y1bar[39] == 4.0)
    assert (y1bar[40] == 3.0)
    assert (y1bar[41] == 2.0)
    assert (y1bar[42] == 1.0)

def test_CompositeGroup_get_matrix():
    """
    Tests the get_matrix operation.
    """

    from lielab.domain import CompositeGroup

    y1 = _make_cgroup()

    y1hat = y1.get_matrix()

    assert (y1hat.shape[0] == 19)
    assert (y1hat.shape[1] == 19)
    
    # TODO: Check the 0's
    # TODO: Check each submatrix

def test_CompositeGroup_operator_parenthesis():
    from lielab.domain import CompositeGroup

    y1 = _make_cgroup()

    # In bounds CN component
    assert (y1(0, 0) == complex(1.0, 0.0))
    assert (y1(0, 1) == complex(0.0, 0.0))
    assert (y1(0, 2) == complex(1.0, 2.0))
    assert (y1(1, 0) == complex(0.0, 0.0))
    assert (y1(1, 1) == complex(1.0, 0.0))
    assert (y1(1, 2) == complex(3.0, 4.0))
    assert (y1(2, 0) == complex(0.0, 0.0))
    assert (y1(2, 1) == complex(0.0, 0.0))
    assert (y1(2, 2) == complex(1.0, 0.0))
    assert (y1(-19, -19) == complex(1.0, 0.0))
    assert (y1(-19, -18) == complex(0.0, 0.0))
    assert (y1(-19, -17) == complex(1.0, 2.0))
    assert (y1(-18, -19) == complex(0.0, 0.0))
    assert (y1(-18, -18) == complex(1.0, 0.0))
    assert (y1(-18, -17) == complex(3.0, 4.0))
    assert (y1(-17, -19) == complex(0.0, 0.0))
    assert (y1(-17, -18) == complex(0.0, 0.0))
    assert (y1(-17, -17) == complex(1.0, 0.0))

    # In bounds GLC component
    assert (y1(3, 3) == complex(5.0, 6.0))
    assert (y1(3, 4) == complex(7.0, 8.0))
    assert (y1(4, 3) == complex(9.0, 10.0))
    assert (y1(4, 4) == complex(11.0, 12.0))
    assert (y1(-16, -16) == complex(5.0, 6.0))
    assert (y1(-16, -15) == complex(7.0, 8.0))
    assert (y1(-15, -16) == complex(9.0, 10.0))
    assert (y1(-15, -15) == complex(11.0, 12.0))

    # In bounds GLR component
    assert (y1(5, 5) == complex(13.0, 0.0))
    assert (y1(5, 6) == complex(14.0, 0.0))
    assert (y1(6, 5) == complex(15.0, 0.0))
    assert (y1(6, 6) == complex(16.0, 0.0))
    assert (y1(-14, -14) == complex(13.0, 0.0))
    assert (y1(-14, -13) == complex(14.0, 0.0))
    assert (y1(-13, -14) == complex(15.0, 0.0))
    assert (y1(-13, -13) == complex(16.0, 0.0))

    # In bounds RN component
    assert (y1(7, 7) == complex(1.0, 0.0))
    assert (y1(7, 8) == complex(0.0, 0.0))
    assert (y1(7, 9) == complex(17.0, 0.0))
    assert (y1(8, 7) == complex(0.0, 0.0))
    assert (y1(8, 8) == complex(1.0, 0.0))
    assert (y1(8, 9) == complex(18.0, 0.0))
    assert (y1(9, 7) == complex(0.0, 0.0))
    assert (y1(9, 8) == complex(0.0, 0.0))
    assert (y1(9, 9) == complex(1.0, 0.0))
    assert (y1(-12, -12) == complex(1.0, 0.0))
    assert (y1(-12, -11) == complex(0.0, 0.0))
    assert (y1(-12, -10) == complex(17.0, 0.0))
    assert (y1(-11, -12) == complex(0.0, 0.0))
    assert (y1(-11, -11) == complex(1.0, 0.0))
    assert (y1(-11, -10) == complex(18.0, 0.0))
    assert (y1(-10, -12) == complex(0.0, 0.0))
    assert (y1(-10, -11) == complex(0.0, 0.0))
    assert (y1(-10, -10) == complex(1.0, 0.0))

    # In bounds SE component
    assert (y1(10, 10) == complex(19.0, 0.0))
    assert (y1(10, 11) == complex(20.0, 0.0))
    assert (y1(10, 12) == complex(21.0, 0.0))
    assert (y1(11, 10) == complex(22.0, 0.0))
    assert (y1(11, 11) == complex(23.0, 0.0))
    assert (y1(11, 12) == complex(24.0, 0.0))
    assert (y1(12, 10) == complex(25.0, 0.0))
    assert (y1(12, 11) == complex(26.0, 0.0))
    assert (y1(12, 12) == complex(27.0, 0.0))
    assert (y1(-9, -9) == complex(19.0, 0.0))
    assert (y1(-9, -8) == complex(20.0, 0.0))
    assert (y1(-9, -7) == complex(21.0, 0.0))
    assert (y1(-8, -9) == complex(22.0, 0.0))
    assert (y1(-8, -8) == complex(23.0, 0.0))
    assert (y1(-8, -7) == complex(24.0, 0.0))
    assert (y1(-7, -9) == complex(25.0, 0.0))
    assert (y1(-7, -8) == complex(26.0, 0.0))
    assert (y1(-7, -7) == complex(27.0, 0.0))

    # In bounds SO component
    assert (y1(13, 13) == complex(28.0, 0.0))
    assert (y1(13, 14) == complex(29.0, 0.0))
    assert (y1(14, 13) == complex(30.0, 0.0))
    assert (y1(14, 14) == complex(31.0, 0.0))
    assert (y1(-6, -6) == complex(28.0, 0.0))
    assert (y1(-6, -5) == complex(29.0, 0.0))
    assert (y1(-5, -6) == complex(30.0, 0.0))
    assert (y1(-5, -5) == complex(31.0, 0.0))

    # In bounds SP component
    assert (y1(15, 15) == complex(32.0, 0.0))
    assert (y1(15, 16) == complex(33.0, 0.0))
    assert (y1(16, 15) == complex(34.0, 0.0))
    assert (y1(16, 16) == complex(35.0, 0.0))
    assert (y1(-4, -4) == complex(32.0, 0.0))
    assert (y1(-4, -3) == complex(33.0, 0.0))
    assert (y1(-3, -4) == complex(34.0, 0.0))
    assert (y1(-3, -3) == complex(35.0, 0.0))

    # In bounds SU component
    assert (y1(17, 17) == complex(36.0, 37.0))
    assert (y1(17, 18) == complex(38.0, 39.0))
    assert (y1(18, 17) == complex(40.0, 41.0))
    assert (y1(18, 18) == complex(42.0, 43.0))
    assert (y1(-2, -2) == complex(36.0, 37.0))
    assert (y1(-2, -1) == complex(38.0, 39.0))
    assert (y1(-1, -2) == complex(40.0, 41.0))
    assert (y1(-1, -1) == complex(42.0, 43.0))

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
    assert (np.isnan(np.real(y1(0, -20))))
    assert (np.isnan(np.imag(y1(0, -20))))
    assert (np.isnan(np.real(y1(-20, 0))))
    assert (np.isnan(np.imag(y1(-20, 0))))
    assert (np.isnan(np.real(y1(-20, -20))))
    assert (np.isnan(np.imag(y1(-20, -20))))
    assert (np.isnan(np.real(y1(0, 19))))
    assert (np.isnan(np.imag(y1(0, 19))))
    assert (np.isnan(np.real(y1(19, 0))))
    assert (np.isnan(np.imag(y1(19, 0))))
    assert (np.isnan(np.real(y1(19, 19))))
    assert (np.isnan(np.imag(y1(19, 19))))

def test_CompositeGroup_math_ops_CompositeGroup():
    from lielab.domain import CompositeGroup

    y1 = _make_cgroup()
    y2 = _make_cgroup()

    y1_prod_y2 = y1*y2
    assert (len(y1_prod_y2.space) == len(y1.space))

    y1_prod_y2_0bar = y1_prod_y2[0].serialize()
    assert (y1_prod_y2_0bar.size == 4)
    assert (y1_prod_y2_0bar[0] == 2.0)
    assert (y1_prod_y2_0bar[1] == 4.0)
    assert (y1_prod_y2_0bar[2] == 6.0)
    assert (y1_prod_y2_0bar[3] == 8.0)

    y1_prod_y2_1bar = y1_prod_y2[1].serialize()
    assert (y1_prod_y2_1bar.size == 8)
    assert (y1_prod_y2_1bar[0] == -28.0)
    assert (y1_prod_y2_1bar[1] == 202.0)
    assert (y1_prod_y2_1bar[2] == -32.0)
    assert (y1_prod_y2_1bar[3] == 254.0)
    assert (y1_prod_y2_1bar[4] == -36.0)
    assert (y1_prod_y2_1bar[5] == 322.0)
    assert (y1_prod_y2_1bar[6] == -40.0)
    assert (y1_prod_y2_1bar[7] == 406.0)

    y1_prod_y2_2bar = y1_prod_y2[2].serialize()
    assert (y1_prod_y2_2bar.size == 4)
    assert (y1_prod_y2_2bar[0] == 379.0)
    assert (y1_prod_y2_2bar[1] == 406.0)
    assert (y1_prod_y2_2bar[2] == 435.0)
    assert (y1_prod_y2_2bar[3] == 466.0)

    y1_prod_y2_3bar = y1_prod_y2[3].serialize()
    assert (y1_prod_y2_3bar.size == 2)
    assert (y1_prod_y2_3bar[0] == 34.0)
    assert (y1_prod_y2_3bar[1] == 36.0)

    y1_prod_y2_4bar = y1_prod_y2[4].serialize()
    assert (y1_prod_y2_4bar.size == 9)
    assert (y1_prod_y2_4bar[0] == 1326.0)
    assert (y1_prod_y2_4bar[1] == 1386.0)
    assert (y1_prod_y2_4bar[2] == 1446.0)
    assert (y1_prod_y2_4bar[3] == 1524.0)
    assert (y1_prod_y2_4bar[4] == 1593.0)
    assert (y1_prod_y2_4bar[5] == 1662.0)
    assert (y1_prod_y2_4bar[6] == 1722.0)
    assert (y1_prod_y2_4bar[7] == 1800.0)
    assert (y1_prod_y2_4bar[8] == 1878.0)

    y1_prod_y2_5bar = y1_prod_y2[5].serialize()
    assert (y1_prod_y2_5bar.size == 4)
    assert (y1_prod_y2_5bar[0] == 1654.0)
    assert (y1_prod_y2_5bar[1] == 1711.0)
    assert (y1_prod_y2_5bar[2] == 1770.0)
    assert (y1_prod_y2_5bar[3] == 1831.0)

    y1_prod_y2_6bar = y1_prod_y2[6].serialize()
    assert (y1_prod_y2_6bar.size == 4)
    assert (y1_prod_y2_6bar[0] == 2146.0)
    assert (y1_prod_y2_6bar[1] == 2211.0)
    assert (y1_prod_y2_6bar[2] == 2278.0)
    assert (y1_prod_y2_6bar[3] == 2347.0)

    y1_prod_y2_7bar = y1_prod_y2[7].serialize()
    assert (y1_prod_y2_7bar.size == 8)
    assert (y1_prod_y2_7bar[0] == -152.0)
    assert (y1_prod_y2_7bar[1] == 5782.0)
    assert (y1_prod_y2_7bar[2] == -156.0)
    assert (y1_prod_y2_7bar[3] == 6082.0)
    assert (y1_prod_y2_7bar[4] == -160.0)
    assert (y1_prod_y2_7bar[5] == 6398.0)
    assert (y1_prod_y2_7bar[6] == -164.0)
    assert (y1_prod_y2_7bar[7] == 6730.0)

    y1 *= y2
    assert (len(y1.space) == len(y2.space))

    y1_iprod_y2_0bar = y1[0].serialize()
    assert (y1_iprod_y2_0bar.size == 4)
    assert (y1_iprod_y2_0bar[0] == 2.0)
    assert (y1_iprod_y2_0bar[1] == 4.0)
    assert (y1_iprod_y2_0bar[2] == 6.0)
    assert (y1_iprod_y2_0bar[3] == 8.0)

    y1_iprod_y2_1bar = y1[1].serialize()
    assert (y1_iprod_y2_1bar.size == 8)
    assert (y1_iprod_y2_1bar[0] == -28.0)
    assert (y1_iprod_y2_1bar[1] == 202.0)
    assert (y1_iprod_y2_1bar[2] == -32.0)
    assert (y1_iprod_y2_1bar[3] == 254.0)
    assert (y1_iprod_y2_1bar[4] == -36.0)
    assert (y1_iprod_y2_1bar[5] == 322.0)
    assert (y1_iprod_y2_1bar[6] == -40.0)
    assert (y1_iprod_y2_1bar[7] == 406.0)

    y1_iprod_y2_2bar = y1[2].serialize()
    assert (y1_iprod_y2_2bar.size == 4)
    assert (y1_iprod_y2_2bar[0] == 379.0)
    assert (y1_iprod_y2_2bar[1] == 406.0)
    assert (y1_iprod_y2_2bar[2] == 435.0)
    assert (y1_iprod_y2_2bar[3] == 466.0)

    y1_iprod_y2_3bar = y1[3].serialize()
    assert (y1_iprod_y2_3bar.size == 2)
    assert (y1_iprod_y2_3bar[0] == 34.0)
    assert (y1_iprod_y2_3bar[1] == 36.0)

    y1_iprod_y2_4bar = y1[4].serialize()
    assert (y1_iprod_y2_4bar.size == 9)
    assert (y1_iprod_y2_4bar[0] == 1326.0)
    assert (y1_iprod_y2_4bar[1] == 1386.0)
    assert (y1_iprod_y2_4bar[2] == 1446.0)
    assert (y1_iprod_y2_4bar[3] == 1524.0)
    assert (y1_iprod_y2_4bar[4] == 1593.0)
    assert (y1_iprod_y2_4bar[5] == 1662.0)
    assert (y1_iprod_y2_4bar[6] == 1722.0)
    assert (y1_iprod_y2_4bar[7] == 1800.0)
    assert (y1_iprod_y2_4bar[8] == 1878.0)

    y1_iprod_y2_5bar = y1[5].serialize()
    assert (y1_iprod_y2_5bar.size == 4)
    assert (y1_iprod_y2_5bar[0] == 1654.0)
    assert (y1_iprod_y2_5bar[1] == 1711.0)
    assert (y1_iprod_y2_5bar[2] == 1770.0)
    assert (y1_iprod_y2_5bar[3] == 1831.0)

    y1_iprod_y2_6bar = y1[6].serialize()
    assert (y1_iprod_y2_6bar.size == 4)
    assert (y1_iprod_y2_6bar[0] == 2146.0)
    assert (y1_iprod_y2_6bar[1] == 2211.0)
    assert (y1_iprod_y2_6bar[2] == 2278.0)
    assert (y1_iprod_y2_6bar[3] == 2347.0)

    y1_iprod_y2_7bar = y1[7].serialize()
    assert (y1_iprod_y2_7bar.size == 8)
    assert (y1_iprod_y2_7bar[0] == -152.0)
    assert (y1_iprod_y2_7bar[1] == 5782.0)
    assert (y1_iprod_y2_7bar[2] == -156.0)
    assert (y1_iprod_y2_7bar[3] == 6082.0)
    assert (y1_iprod_y2_7bar[4] == -160.0)
    assert (y1_iprod_y2_7bar[5] == 6398.0)
    assert (y1_iprod_y2_7bar[6] == -164.0)
    assert (y1_iprod_y2_7bar[7] == 6730.0)

    y1 = _make_cgroup()

    y1_inv = y1.inverse()
    assert (len(y1_inv.space) == len(y1.space))

    y1_inv_0bar = y1_inv[0].serialize()
    assert (y1_inv_0bar.size == 4)
    assert (y1_inv_0bar[0] == -1.0)
    assert (y1_inv_0bar[1] == -2.0)
    assert (y1_inv_0bar[2] == -3.0)
    assert (y1_inv_0bar[3] == -4.0)

    y1_inv_1bar = y1_inv[1].serialize()
    assert (y1_inv_1bar.size == 8)
    assert (np.abs(y1_inv_1bar[0] - -0.75) < 1e-14)
    assert (np.abs(y1_inv_1bar[1] - 0.6875) < 1e-14)
    assert (np.abs(y1_inv_1bar[2] - 0.5) < 1e-14)
    assert (np.abs(y1_inv_1bar[3] - -0.4375) < 1e-14)
    assert (np.abs(y1_inv_1bar[4] - 0.625) < 1e-14)
    assert (np.abs(y1_inv_1bar[5] - -0.5625) < 1e-14)
    assert (np.abs(y1_inv_1bar[6] - -0.375) < 1e-14)
    assert (np.abs(y1_inv_1bar[7] - 0.3125) < 1e-14)

    y1_inv_2bar = y1_inv[2].serialize()
    assert (y1_inv_2bar.size == 4)
    assert (np.abs(y1_inv_2bar[0] - -8.0) < 1e-13)
    assert (np.abs(y1_inv_2bar[1] - 7.0) < 1e-13)
    assert (np.abs(y1_inv_2bar[2] - 7.5) < 1e-13)
    assert (np.abs(y1_inv_2bar[3] - -6.5) < 1e-13)

    y1_inv_3bar = y1_inv[3].serialize()
    assert (y1_inv_3bar.size == 2)
    assert (y1_inv_3bar[0] == -17.0)
    assert (y1_inv_3bar[1] == -18.0)

    y1_inv_4bar = y1_inv[4].serialize()
    assert (y1_inv_4bar.size == 9)
    # TODO: Change default values to something more reasonable
    # assert (y1_inv_4bar[0] == 1326.0)
    # assert (y1_inv_4bar[1] == 1386.0)
    # assert (y1_inv_4bar[2] == 1446.0)
    # assert (y1_inv_4bar[3] == 1524.0)
    # assert (y1_inv_4bar[4] == 1593.0)
    # assert (y1_inv_4bar[5] == 1662.0)
    # assert (y1_inv_4bar[6] == 1722.0)
    # assert (y1_inv_4bar[7] == 1800.0)
    # assert (y1_inv_4bar[8] == 1878.0)

    y1_inv_5bar = y1_inv[5].serialize()
    assert (y1_inv_5bar.size == 4)
    assert (y1_inv_5bar[0] == 28.0)
    assert (y1_inv_5bar[1] == 30.0)
    assert (y1_inv_5bar[2] == 29.0)
    assert (y1_inv_5bar[3] == 31.0)
    
    y1_inv_6bar = y1_inv[6].serialize()
    assert (y1_inv_6bar.size == 4)
    assert (np.abs(y1_inv_6bar[0] - -17.5) < 1e-12)
    assert (np.abs(y1_inv_6bar[1] - 16.5) < 1e-12)
    assert (np.abs(y1_inv_6bar[2] - 17.0) < 1e-12)
    assert (np.abs(y1_inv_6bar[3] - -16.0) < 1e-12)

    y1_inv_7bar = y1_inv[7].serialize()
    assert (y1_inv_7bar.size == 8)
    assert (np.abs(y1_inv_7bar[0] - -2.6875) < 1e-13)
    assert (np.abs(y1_inv_7bar[1] - 2.625) < 1e-13)
    assert (np.abs(y1_inv_7bar[2] - 2.4375) < 1e-13)
    assert (np.abs(y1_inv_7bar[3] - -2.375) < 1e-13)
    assert (np.abs(y1_inv_7bar[4] - 2.5625) < 1e-13)
    assert (np.abs(y1_inv_7bar[5] - -2.5) < 1e-13)
    assert (np.abs(y1_inv_7bar[6] - -2.3125) < 1e-13)
    assert (np.abs(y1_inv_7bar[7] - 2.25) < 1e-13)

def test_CompositeGroup_operator_bracket():
    from lielab.domain import CompositeGroup

    x1 = _make_cgroup()
    
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

    assert (x10.to_string() == "C^2")
    assert (x11.to_string() == "GL(2, C)")
    assert (x12.to_string() == "GL(2, R)")
    assert (x13.to_string() == "R^2")
    assert (x14.to_string() == "SE(2)")
    assert (x15.to_string() == "SO(2)")
    assert (x16.to_string() == "SP(2, R)")
    assert (x17.to_string() == "SU(2)")
    assert (x1m8.to_string() == "C^2")
    assert (x1m7.to_string() == "GL(2, C)")
    assert (x1m6.to_string() == "GL(2, R)")
    assert (x1m5.to_string() == "R^2")
    assert (x1m4.to_string() == "SE(2)")
    assert (x1m3.to_string() == "SO(2)")
    assert (x1m2.to_string() == "SP(2, R)")
    assert (x1m1.to_string() == "SU(2)")

    # Out of bounds
    x18 = x1[8]
    x1m9 = x1[-9]

    assert (x18.to_string() == "GL(0, C)")
    assert (x1m9.to_string() == "GL(0, C)")
