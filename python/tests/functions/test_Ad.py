def test_Ad():
    """
    Tests the Ad function.
    """

    from lielab.domain import so
    from lielab.functions import exp, Ad
    from lielab.testing import check_almost_equal_tol
    import numpy as np

    u = so.from_vector([1.0, 0.0, 0.0])
    v = so.from_vector([0.0, 1.0, 0.0])
    w = so.from_vector([0.0, 0.0, 1.0])

    Gso = exp(v)

    # GuG^-1
    ansso = Ad(Gso, u)
    truthso = np.array([[0, 0.841470984807896, 0],
                        [-0.841470984807897, 0, -0.540302305868140],
                        [0, 0.540302305868140, 0]])
    
    assert check_almost_equal_tol(ansso.get_matrix(), truthso)

    # GvG^-1 = v when G = exp(v)
    ansso = Ad(Gso, v)
    
    assert check_almost_equal_tol(ansso.get_matrix(), v.get_matrix())

    # GwG^-1
    ansso = Ad(Gso, w)
    truthso = np.array([[0, -0.540302305868140, 0],
                        [0.540302305868140, 0, -0.841470984807897],
                        [0, 0.841470984807897, 0]])
    
    assert check_almost_equal_tol(ansso.get_matrix(), truthso)
