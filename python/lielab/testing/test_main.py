import numpy as np

TOL_FINE = 1e-6
TOL_COARSE = 1e-3

an_int = int(2)
a_double = float(2)
an_imag_int = 2j # TODO: I think this is actually just a double
an_imag_double = 2.0j

def assert_matrix(mat1: np.ndarray, mat2: np.ndarray) -> None:
    if isinstance(mat1, np.ndarray):
        pass
    else:
        mat1 = mat1.get_matrix()

    r1 = mat1.shape[0]
    if len(mat1.shape) == 1:
        c1 = 0
    else:
        c1 = mat1.shape[1]
    
    if isinstance(mat2, np.ndarray):
        pass
    else:
        mat2 = mat2.get_matrix()
    
    r2 = mat2.shape[0]
    if len(mat2.shape) == 1:
        c2 = 0
    else:
        c2 = mat2.shape[1]
    
    assert r1 == r2
    assert c1 == c2

    if (c1 != 0 and c2 != 0):
        for ii in range(r1):
            for jj in range(c1):
                assert abs(mat1[ii,jj] - mat2[ii,jj]) < TOL_FINE
    else:
        for ii in range(r1):
            assert abs(mat1[ii] - mat2[ii]) < TOL_FINE

def assert_domain(d1, d2):
    assert_matrix(d1.get_matrix(), d2.get_matrix())
