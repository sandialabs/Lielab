import numpy as np
import pytest

def complex(a,b):
    return a + b*1j

def _make_cmanifold():
    from lielab.domain import CompositeManifold
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
    ySE1.unserialize([19.0, 20.0, 21.0, 22.0, 23.0, 24.0])
    ySO1 = lielab.domain.SO(2)
    ySO1.unserialize([28.0, 29.0, 30.0, 31.0])
    ySP1 = lielab.domain.SP(2)
    ySP1.unserialize([32.0, 33.0, 34.0, 35.0])
    ySU1 = lielab.domain.SU(2)
    ySU1.unserialize([36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0])

    ycn1 = lielab.domain.cn.from_vector([44.0, 45.0, 46.0, 47.0])
    yglc1 = lielab.domain.glc.from_vector([48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0])
    yglr1 = lielab.domain.glr.from_vector([56.0, 57.0, 58.0, 59.0])
    yrn1 = lielab.domain.rn.from_vector([60.0, 61.0])
    yse1 = lielab.domain.se.from_vector([62.0, 63.0, 64.0])
    yso1 = lielab.domain.so.from_vector([65.0, 66.0, 67.0])
    ysp1 = lielab.domain.sp.from_vector([68.0, 69.0, 70.0])
    ysu1 = lielab.domain.su.from_vector([71.0, 72.0, 73.0])

    return CompositeManifold([yCN1, yGLC1, yGLR1, yRN1, ySE1, ySO1, ySP1, ySU1, ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1])

def test_CompositeManifold_to_string():
    from lielab.domain import CompositeManifold

    y1 = _make_cmanifold()

    assert (y1.to_string() == "C^2 x GL(2, C) x GL(2, R) x R^2 x SE(2) x SO(2) x SP(2, R) x SU(2) x c^2 x gl(2, C) x gl(2, R) x r^2 x se(2) x so(3) x sp(2, R) x su(2)")

def test_CompositeManifold_main_initializer():
    from lielab.domain import CompositeManifold

    xblank = CompositeManifold()
    assert (xblank.get_dimension() == 0)

    x0 = CompositeManifold(0)
    assert (len(x0) == 0)
    x1 = CompositeManifold(1)
    assert (len(x1) == 1)
    x10 = CompositeManifold(10)
    assert (len(x10) == 10)

def test_CompositeManifold_get_dimension():
    from lielab.domain import CompositeManifold

    y1 = _make_cmanifold()
    
    dims = y1.get_dimensions()
    assert (len(dims) == 16)
    assert (dims[0] == 4)
    assert (dims[1] == 8)
    assert (dims[2] == 4)
    assert (dims[3] == 2)
    assert (dims[4] == 3)
    assert (dims[5] == 1)
    assert (dims[6] == 3)
    assert (dims[7] == 3)
    assert (dims[8] == 4)
    assert (dims[9] == 8)
    assert (dims[10] == 4)
    assert (dims[11] == 2)
    assert (dims[12] == 3)
    assert (dims[13] == 3)
    assert (dims[14] == 3)
    assert (dims[15] == 3)

    assert (y1.get_dimension() == 58)

def test_CompositeManifold_get_size():
    from lielab.domain import CompositeManifold

    y1 = _make_cmanifold()
    
    sizes = y1.get_sizes()
    assert (len(sizes) == 16)
    assert (sizes[0] == 4)
    assert (sizes[1] == 8)
    assert (sizes[2] == 4)
    assert (sizes[3] == 2)
    assert (sizes[4] == 6)
    assert (sizes[5] == 4)
    assert (sizes[6] == 4)
    assert (sizes[7] == 8)
    assert (sizes[8] == 4)
    assert (sizes[9] == 8)
    assert (sizes[10] == 4)
    assert (sizes[11] == 2)
    assert (sizes[12] == 3)
    assert (sizes[13] == 3)
    assert (sizes[14] == 3)
    assert (sizes[15] == 3)

    assert (y1.get_size() == 70)

def test_CompositeManifold_serialize_unserialize():
    """
    Tests the serialize/unserialize operation.
    """

    from lielab.domain import CompositeManifold

    y1 = _make_cmanifold()
    y1.unserialize([73.0, 72.0, 71.0, 70.0, 69.0, 68.0, 67.0, 66.0, 65.0, 64.0, 63.0, 62.0, 61.0, 60.0,
                    59.0, 58.0, 57.0, 56.0, 55.0, 54.0, 53.0, 52.0, 51.0, 50.0, 49.0, 48.0, 47.0, 46.0,
                    45.0, 44.0, 43.0, 42.0, 41.0, 40.0, 39.0, 38.0, 37.0, 36.0, 35.0, 34.0, 33.0, 32.0,
                    31.0, 30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0, 20.0, 19.0, 18.0,
                    17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0])
    
    y1bar = y1.serialize()
    assert (y1bar.size == 70)
    assert (y1bar[0] == 73.0)
    assert (y1bar[1] == 72.0)
    assert (y1bar[2] == 71.0)
    assert (y1bar[3] == 70.0)
    assert (y1bar[4] == 69.0)
    assert (y1bar[5] == 68.0)
    assert (y1bar[6] == 67.0)
    assert (y1bar[7] == 66.0)
    assert (y1bar[8] == 65.0)
    assert (y1bar[9] == 64.0)
    assert (y1bar[10] == 63.0)
    assert (y1bar[11] == 62.0)
    assert (y1bar[12] == 61.0)
    assert (y1bar[13] == 60.0)
    assert (y1bar[14] == 59.0)
    assert (y1bar[15] == 58.0)
    assert (y1bar[16] == 57.0)
    assert (y1bar[17] == 56.0)
    assert (y1bar[18] == 55.0)
    assert (y1bar[19] == 54.0)
    assert (y1bar[20] == 53.0)
    assert (y1bar[21] == 52.0)
    assert (y1bar[22] == 51.0)
    assert (y1bar[23] == 50.0)
    assert (y1bar[24] == 49.0)
    assert (y1bar[25] == 48.0)
    assert (y1bar[26] == 47.0)
    assert (y1bar[27] == 46.0)
    assert (y1bar[28] == 45.0)
    assert (y1bar[29] == 44.0)
    assert (y1bar[30] == 43.0)
    assert (y1bar[31] == 42.0)
    assert (y1bar[32] == 41.0)
    assert (y1bar[33] == 40.0)
    assert (y1bar[34] == 39.0)
    assert (y1bar[35] == 38.0)
    assert (y1bar[36] == 37.0)
    assert (y1bar[37] == 36.0)
    assert (y1bar[38] == 35.0)
    assert (y1bar[39] == 34.0)
    assert (y1bar[40] == 33.0)
    assert (y1bar[41] == 32.0)
    assert (y1bar[42] == 31.0)
    assert (y1bar[43] == 30.0)
    assert (y1bar[44] == 29.0)
    assert (y1bar[45] == 28.0)
    assert (y1bar[46] == 27.0)
    assert (y1bar[47] == 26.0)
    assert (y1bar[48] == 25.0)
    assert (y1bar[49] == 24.0)
    assert (y1bar[50] == 23.0)
    assert (y1bar[51] == 22.0)
    assert (y1bar[52] == 21.0)
    assert (y1bar[53] == 20.0)
    assert (y1bar[54] == 19.0)
    assert (y1bar[55] == 18.0)
    assert (y1bar[56] == 17.0)
    assert (y1bar[57] == 16.0)
    assert (y1bar[58] == 15.0)
    assert (y1bar[59] == 14.0)
    assert (y1bar[60] == 13.0)
    assert (y1bar[61] == 12.0)
    assert (y1bar[62] == 11.0)
    assert (y1bar[63] == 10.0)
    assert (y1bar[64] == 9.0)
    assert (y1bar[65] == 8.0)
    assert (y1bar[66] == 7.0)
    assert (y1bar[67] == 6.0)
    assert (y1bar[68] == 5.0)
    assert (y1bar[69] == 4.0)

def test_CompositeManifold_operator_bracket():
    from lielab.domain import CompositeManifold

    x1 = _make_cmanifold()
    
    x10 = x1[0]
    x11 = x1[1]
    x12 = x1[2]
    x13 = x1[3]
    x14 = x1[4]
    x15 = x1[5]
    x16 = x1[6]
    x17 = x1[7]
    x18 = x1[8]
    x19 = x1[9]
    x110 = x1[10]
    x111 = x1[11]
    x112 = x1[12]
    x113 = x1[13]
    x114 = x1[14]
    x115 = x1[15]
    x1m16 = x1[-16]
    x1m15 = x1[-15]
    x1m14 = x1[-14]
    x1m13 = x1[-13]
    x1m12 = x1[-12]
    x1m11 = x1[-11]
    x1m10 = x1[-10]
    x1m9 = x1[-9]
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
    assert (x18.to_string() == "c^2")
    assert (x19.to_string() == "gl(2, C)")
    assert (x110.to_string() == "gl(2, R)")
    assert (x111.to_string() == "r^2")
    assert (x112.to_string() == "se(2)")
    assert (x113.to_string() == "so(3)")
    assert (x114.to_string() == "sp(2, R)")
    assert (x115.to_string() == "su(2)")
    assert (x1m16.to_string() == "C^2")
    assert (x1m15.to_string() == "GL(2, C)")
    assert (x1m14.to_string() == "GL(2, R)")
    assert (x1m13.to_string() == "R^2")
    assert (x1m12.to_string() == "SE(2)")
    assert (x1m11.to_string() == "SO(2)")
    assert (x1m10.to_string() == "SP(2, R)")
    assert (x1m9.to_string() == "SU(2)")
    assert (x1m8.to_string() == "c^2")
    assert (x1m7.to_string() == "gl(2, C)")
    assert (x1m6.to_string() == "gl(2, R)")
    assert (x1m5.to_string() == "r^2")
    assert (x1m4.to_string() == "se(2)")
    assert (x1m3.to_string() == "so(3)")
    assert (x1m2.to_string() == "sp(2, R)")
    assert (x1m1.to_string() == "su(2)")

    # Out of bounds
    with pytest.raises(RuntimeError):
        x116 = x1[16]
    
    with pytest.raises(RuntimeError):
        x1m17 = x1[-17]
