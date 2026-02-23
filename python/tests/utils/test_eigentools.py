
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

# TODO: Test logspace

def test_column_stack():
    from lielab.utils import column_stack
    import numpy as np

    li = []
    ld = []

    assert (column_stack(li).size == 0)
    assert (column_stack(ld).size == 0)

    i1 = np.zeros(3, dtype=np.int32)
    i1[0] = 1
    i1[1] = 2
    i1[2] = 3

    i2 = np.zeros(2, dtype=np.int32)
    i2[0] = 4
    i2[1] = 5

    t1 = column_stack([i1])
    assert (t1.shape[0] == 3)
    assert (t1.shape[1] == 1)
    assert (t1[0, 0] == 1)
    assert (t1[1, 0] == 2)
    assert (t1[2, 0] == 3)

    t2 = column_stack([i1, i2])
    assert (t2.shape[0] == 2)
    assert (t2.shape[1] == 2)
    assert (t2[0, 0] == 1)
    assert (t2[1, 0] == 2)
    assert (t2[0, 1] == 4)
    assert (t2[1, 1] == 5)

    t3 = column_stack([i1, i2, i1])
    assert (t3.shape[0] == 2)
    assert (t3.shape[1] == 3)
    assert (t3[0, 0] == 1)
    assert (t3[1, 0] == 2)
    assert (t3[0, 1] == 4)
    assert (t3[1, 1] == 5)
    assert (t3[0, 2] == 1)
    assert (t3[1, 2] == 2)

    li.append(i1)
    t4 = column_stack(li)
    assert (t4.shape[0] == 3)
    assert (t4.shape[1] == 1)
    assert (t4[0, 0] == 1)
    assert (t4[1, 0] == 2)
    assert (t4[2, 0] == 3)

    li.append(i2)
    t5 = column_stack(li)
    assert (t5.shape[0] == 2)
    assert (t5.shape[1] == 2)
    assert (t5[0, 0] == 1)
    assert (t5[1, 0] == 2)
    assert (t5[0, 1] == 4)
    assert (t5[1, 1] == 5)

    li.append(i1)
    t6 = column_stack(li)
    assert (t6.shape[0] == 2)
    assert (t6.shape[1] == 3)
    assert (t6[0, 0] == 1)
    assert (t6[1, 0] == 2)
    assert (t6[0, 1] == 4)
    assert (t6[1, 1] == 5)
    assert (t6[0, 2] == 1)
    assert (t6[1, 2] == 2)

    d1 = np.zeros(3, dtype=np.float64)
    d1[0] = 1.5
    d1[1] = 2.5
    d1[2] = 3.5

    d2 = np.zeros(2, dtype=np.float64)
    d2[0] = 4.5
    d2[1] = 5.5

    t7 = column_stack([d1])
    assert (t7.shape[0] == 3)
    assert (t7.shape[1] == 1)
    assert (t7[0, 0] == 1.5)
    assert (t7[1, 0] == 2.5)
    assert (t7[2, 0] == 3.5)

    t8 = column_stack([d1, d2])
    assert (t8.shape[0] == 2)
    assert (t8.shape[1] == 2)
    assert (t8[0, 0] == 1.5)
    assert (t8[1, 0] == 2.5)
    assert (t8[0, 1] == 4.5)
    assert (t8[1, 1] == 5.5)

    t9 = column_stack([d1, d2, d1])
    assert (t9.shape[0] == 2)
    assert (t9.shape[1] == 3)
    assert (t9[0, 0] == 1.5)
    assert (t9[1, 0] == 2.5)
    assert (t9[0, 1] == 4.5)
    assert (t9[1, 1] == 5.5)
    assert (t9[0, 2] == 1.5)
    assert (t9[1, 2] == 2.5)

    ld.append(d1)
    t10 = column_stack(ld)
    assert (t10.shape[0] == 3)
    assert (t10.shape[1] == 1)
    assert (t10[0, 0] == 1.5)
    assert (t10[1, 0] == 2.5)
    assert (t10[2, 0] == 3.5)

    ld.append(d2)
    t11 = column_stack(ld)
    assert (t11.shape[0] == 2)
    assert (t11.shape[1] == 2)
    assert (t11[0, 0] == 1.5)
    assert (t11[1, 0] == 2.5)
    assert (t11[0, 1] == 4.5)
    assert (t11[1, 1] == 5.5)

    ld.append(d1)
    t12 = column_stack(ld)
    assert (t12.shape[0] == 2)
    assert (t12.shape[1] == 3)
    assert (t12[0, 0] == 1.5)
    assert (t12[1, 0] == 2.5)
    assert (t12[0, 1] == 4.5)
    assert (t12[1, 1] == 5.5)
    assert (t12[0, 2] == 1.5)
    assert (t12[1, 2] == 2.5)

def test_horizontal_stack():
    from lielab.utils import horizontal_stack
    import numpy as np

    d22 = np.array([[1.0, 2.0], [3.0, 4.0]])
    d32 = np.array([[1.1, 2.1], [3.1, 4.1], [5.1, 6.1]])
    d23 = np.array([[1.2, 2.2, 3.2], [4.2, 5.2, 6.2]])

    vlist = []

    t0 = horizontal_stack(vlist)
    assert (t0.shape[0] == 0)
    assert (t0.shape[1] == 0)

    vlist.append(d32)

    t1 = horizontal_stack(vlist)
    assert (t1.shape[0] == 3)
    assert (t1.shape[1] == 2)
    assert (t1[0,0] == 1.1)
    assert (t1[0,1] == 2.1)
    assert (t1[1,0] == 3.1)
    assert (t1[1,1] == 4.1)
    assert (t1[2,0] == 5.1)
    assert (t1[2,1] == 6.1)

    vlist.append(d23)

    t2 = horizontal_stack(vlist)
    assert (t2.shape[0] == 2)
    assert (t2.shape[1] == 5)
    assert (t2[0,0] == 1.1)
    assert (t2[0,1] == 2.1)
    assert (t2[0,2] == 1.2)
    assert (t2[0,3] == 2.2)
    assert (t2[0,4] == 3.2)
    assert (t2[1,0] == 3.1)
    assert (t2[1,1] == 4.1)
    assert (t2[1,2] == 4.2)
    assert (t2[1,3] == 5.2)
    assert (t2[1,4] == 6.2)

    vlist.append(d22)

    t3 = horizontal_stack(vlist)
    assert (t3.shape[0] == 2)
    assert (t3.shape[1] == 7)
    assert (t3[0,0] == 1.1)
    assert (t3[0,1] == 2.1)
    assert (t3[0,2] == 1.2)
    assert (t3[0,3] == 2.2)
    assert (t3[0,4] == 3.2)
    assert (t3[0,5] == 1.0)
    assert (t3[0,6] == 2.0)
    assert (t3[1,0] == 3.1)
    assert (t3[1,1] == 4.1)
    assert (t3[1,2] == 4.2)
    assert (t3[1,3] == 5.2)
    assert (t3[1,4] == 6.2)
    assert (t3[1,5] == 3.0)
    assert (t3[1,6] == 4.0)

def test_vertical_stack():
    from lielab.utils import vertical_stack
    import numpy as np

    d22 = np.array([[1.0, 2.0], [3.0, 4.0]])
    d23 = np.array([[1.1, 2.1, 3.1], [4.1, 5.1, 6.1]])
    d32 = np.array([[1.2, 2.2], [3.2, 4.2], [5.2, 6.2]])

    vlist = []

    t0 = vertical_stack(vlist)
    assert (t0.shape[0] == 0)
    assert (t0.shape[1] == 0)

    vlist.append(d23)

    t1 = vertical_stack(vlist)
    assert (t1.shape[0] == 2)
    assert (t1.shape[1] == 3)
    assert (t1[0,0] == 1.1)
    assert (t1[0,1] == 2.1)
    assert (t1[0,2] == 3.1)
    assert (t1[1,0] == 4.1)
    assert (t1[1,1] == 5.1)
    assert (t1[1,2] == 6.1)

    vlist.append(d32)

    t2 = vertical_stack(vlist)
    assert (t2.shape[0] == 5)
    assert (t2.shape[1] == 2)
    assert (t2[0,0] == 1.1)
    assert (t2[0,1] == 2.1)
    assert (t2[1,0] == 4.1)
    assert (t2[1,1] == 5.1)
    assert (t2[2,0] == 1.2)
    assert (t2[2,1] == 2.2)
    assert (t2[3,0] == 3.2)
    assert (t2[3,1] == 4.2)
    assert (t2[4,0] == 5.2)
    assert (t2[4,1] == 6.2)

    vlist.append(d22)

    t3 = vertical_stack(vlist)
    assert (t3.shape[0] == 7)
    assert (t3.shape[1] == 2)
    assert (t3[0,0] == 1.1)
    assert (t3[0,1] == 2.1)
    assert (t3[1,0] == 4.1)
    assert (t3[1,1] == 5.1)
    assert (t3[2,0] == 1.2)
    assert (t3[2,1] == 2.2)
    assert (t3[3,0] == 3.2)
    assert (t3[3,1] == 4.2)
    assert (t3[4,0] == 5.2)
    assert (t3[4,1] == 6.2)
    assert (t3[5,0] == 1.0)
    assert (t3[5,1] == 2.0)
    assert (t3[6,0] == 3.0)
    assert (t3[6,1] == 4.0)
