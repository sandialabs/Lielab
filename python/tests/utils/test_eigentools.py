
def test_concatenate():
    from lielab.utils import concatenate
    import numpy as np

    li = []
    ld = []

    assert (concatenate(li).size == 0)
    assert (concatenate(ld).size == 0)

    i1 = np.array([1, 2, 3])

    d1 = np.array([1.5, 2.5, 3.5])

    t1 = concatenate([i1])
    assert (t1.size == 3)
    assert (t1[0] == 1)
    assert (t1[1] == 2)
    assert (t1[2] == 3)

    t2 = concatenate([d1])
    assert (t2.size == 3)
    assert (t2[0] == 1.5)
    assert (t2[1] == 2.5)
    assert (t2[2] == 3.5)

    i2 = np.array([4, 5, 6, 7])

    d2 = np.array([4.5, 5.5, 6.5, 7.5])

    t5 = concatenate([i1, i2])
    assert (t5.size == 7)
    assert (t5[0] == 1)
    assert (t5[1] == 2)
    assert (t5[2] == 3)
    assert (t5[3] == 4)
    assert (t5[4] == 5)
    assert (t5[5] == 6)
    assert (t5[6] == 7)

    t6 = concatenate([d1, d2])
    assert (t6.size == 7)
    assert (t6[0] == 1.5)
    assert (t6[1] == 2.5)
    assert (t6[2] == 3.5)
    assert (t6[3] == 4.5)
    assert (t6[4] == 5.5)
    assert (t6[5] == 6.5)
    assert (t6[6] == 7.5)

def test_arange():
    from lielab.utils import arange

    assert (arange(0, 0).size == 0)
    assert (arange(0.0, 0.0).size == 0)
    assert (arange(0).size == 0)
    assert (arange(0.0).size == 0)

    assert (arange(1, 0).size == 0)
    assert (arange(1.0, 0.0).size == 0)

    t1 = arange(0, 1)
    assert (t1.size == 1)
    assert (t1[0] == 0)
    t2 = arange(0.5, 1.5)
    assert (t2.size == 1)
    assert (t2[0] == 0.5)
    t3 = arange(1)
    assert (t3.size == 1)
    assert (t3[0] == 0)
    t4 = arange(0.5)
    assert (t4.size == 1)
    assert (t4[0] == 0.0)

    t5 = arange(2, 6)
    assert (t5.size == 4)
    assert (t5[0] == 2)
    assert (t5[1] == 3)
    assert (t5[2] == 4)
    assert (t5[3] == 5)
    t6 = arange(2.5, 6.5)
    assert (t6.size == 4)
    assert (t6[0] == 2.5)
    assert (t6[1] == 3.5)
    assert (t6[2] == 4.5)
    assert (t6[3] == 5.5)
    t7 = arange(4)
    assert (t7.size == 4)
    assert (t7[0] == 0)
    assert (t7[1] == 1)
    assert (t7[2] == 2)
    assert (t7[3] == 3)
    t8 = arange(3.5)
    assert (t8.size == 4)
    assert (t8[0] == 0.0)
    assert (t8[1] == 1.0)
    assert (t8[2] == 2.0)
    assert (t8[3] == 3.0)

    t9 = arange(-6, -2)
    assert (t9.size == 4)
    assert (t9[0] == -6)
    assert (t9[1] == -5)
    assert (t9[2] == -4)
    assert (t9[3] == -3)
    t10 = arange(-6.5, -2.5)
    assert (t10.size == 4)
    assert (t10[0] == -6.5)
    assert (t10[1] == -5.5)
    assert (t10[2] == -4.5)
    assert (t10[3] == -3.5)
    t11 = arange(-4)
    assert (t11.size == 0)
    t12 = arange(-3.5)
    assert (t12.size == 0)

    t13 = arange(2.5, 6.5 + 1e-8)
    assert (t13.size == 5)
    assert (t13[0] == 2.5)
    assert (t13[1] == 3.5)
    assert (t13[2] == 4.5)
    assert (t13[3] == 5.5)
    assert (t13[4] == 6.5)
    t14 = arange(4.0 + 1e-8)
    assert (t14.size == 5)
    assert (t14[0] == 0.0)
    assert (t14[1] == 1.0)
    assert (t14[2] == 2.0)
    assert (t14[3] == 3.0)
    assert (t14[4] == 4.0)

def test_repeat():
    from lielab.utils import repeat
    import numpy as np

    assert (repeat([], -1).size == 0)
    assert (repeat([], -1).size == 0)

    assert (repeat([], 0).size == 0)
    assert (repeat([], 0).size == 0)

    assert (repeat([], 1).size == 0)
    assert (repeat([], 1).size == 0)

    assert (repeat([], 2).size == 0)
    assert (repeat([], 2).size == 0)

    t1 = repeat([1, 2, 3], -1)
    assert (t1.size == 0)
    t2 = repeat([1, 2, 3], 0)
    assert (t2.size == 0)
    t3 = repeat([1, 2, 3], 1)
    assert (t3.size == 3)
    assert (t3[0] == 1)
    assert (t3[1] == 2)
    assert (t3[2] == 3)
    t4 = repeat([1, 2, 3], 2)
    assert (t4.size == 6)
    assert (t4[0] == 1)
    assert (t4[1] == 1)
    assert (t4[2] == 2)
    assert (t4[3] == 2)
    assert (t4[4] == 3)
    assert (t4[5] == 3)

    v = np.array([1.5, 2.5, 3.5])

    t9 = repeat(v, -1)
    assert (t9.size == 0)
    t10 = repeat(v, 0)
    assert (t10.size == 0)
    t11 = repeat(v, 1)
    assert (t11.size == 3)
    assert (t11[0] == 1.5)
    assert (t11[1] == 2.5)
    assert (t11[2] == 3.5)
    t12 = repeat(v, 2)
    assert (t12.size == 6)
    assert (t12[0] == 1.5)
    assert (t12[1] == 1.5)
    assert (t12[2] == 2.5)
    assert (t12[3] == 2.5)
    assert (t12[4] == 3.5)
    assert (t12[5] == 3.5)


def test_tile():
    from lielab.utils import tile
    import numpy as np

    assert (tile([], -1).size == 0)
    assert (tile([], -1).size == 0)

    assert (tile([], 0).size == 0)
    assert (tile([], 0).size == 0)

    assert (tile([], 1).size == 0)
    assert (tile([], 1).size == 0)

    assert (tile([], 2).size == 0)
    assert (tile([], 2).size == 0)

    t1 = tile([1, 2, 3], -1)
    assert (t1.size == 0)
    t2 = tile([1, 2, 3], 0)
    assert (t2.size == 0)
    t3 = tile([1, 2, 3], 1)
    assert (t3.size == 3)
    assert (t3[0] == 1)
    assert (t3[1] == 2)
    assert (t3[2] == 3)
    t4 = tile([1, 2, 3], 2)
    assert (t4.size == 6)
    assert (t4[0] == 1)
    assert (t4[1] == 2)
    assert (t4[2] == 3)
    assert (t4[3] == 1)
    assert (t4[4] == 2)
    assert (t4[5] == 3)

    v = np.array([1.5, 2.5, 3.5])

    t9 = tile(v, -1)
    assert (t9.size == 0)
    t10 = tile(v, 0)
    assert (t10.size == 0)
    t11 = tile(v, 1)
    assert (t11.size == 3)
    assert (t11[0] == 1.5)
    assert (t11[1] == 2.5)
    assert (t11[2] == 3.5)
    t12 = tile(v, 2)
    assert (t12.size == 6)
    assert (t12[0] == 1.5)
    assert (t12[1] == 2.5)
    assert (t12[2] == 3.5)
    assert (t12[3] == 1.5)
    assert (t12[4] == 2.5)
    assert (t12[5] == 3.5)

def test_linspace():
    from lielab.utils import linspace

    t1 = linspace(1.0, 5.0, -1)
    assert (t1.size == 0)

    t2 = linspace(1.0, 5.0, 0)
    assert (t2.size == 0)

    t3 = linspace(1.0, 5.0, 1)
    assert (t3.size == 1)
    assert (t3[0] == 1.0)
    
    t4 = linspace(1.0, 5.0, 2)
    assert (t4.size == 2)
    assert (t4[0] == 1.0)
    assert (t4[1] == 5.0)

    t5 = linspace(1.0, 5.0, 3)
    assert (t5.size == 3)
    assert (t5[0] == 1.0)
    assert (t5[1] == 3.0)
    assert (t5[2] == 5.0)

    t6 = linspace(1.0, 5.0, 5)
    assert (t6.size == 5)
    assert (t6[0] == 1.0)
    assert (t6[1] == 2.0)
    assert (t6[2] == 3.0)
    assert (t6[3] == 4.0)
    assert (t6[4] == 5.0)
