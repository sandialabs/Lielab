from lielab.testing import *

def test_Ad():
    """
    Tests the Ad function.
    """

    from lielab.domain import so, SO
    from lielab.functions import exp, Ad

    u = so(3)
    v = so(3)
    w = so(3)
    ansso = so(3)
    Gso = SO(3)

    u.set_vector([1,0,0])
    v.set_vector([0,1,0])
    w.set_vector([0,0,1])
    Gso = exp(v)

    # GuG^-1
    ansso = Ad(Gso, u)
    truthso = np.array([[0, 0.841470984807896, 0],
                        [-0.841470984807897, 0, -0.540302305868140],
                        [0, 0.540302305868140, 0]])
    
    assert_matrix(ansso.get_matrix(), truthso)

    # GvG^-1 = v when G = exp(v)
    ansso = Ad(Gso, v)
    
    assert_matrix(ansso.get_matrix(), v.get_matrix())

    # GwG^-1
    ansso = Ad(Gso, w)
    truthso = np.array([[0, -0.540302305868140, 0],
                        [0.540302305868140, 0, -0.841470984807897],
                        [0, 0.841470984807897, 0]])
    
    assert_matrix(ansso.get_matrix(), truthso)
