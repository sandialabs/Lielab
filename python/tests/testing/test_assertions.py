
def test_check_topology():
    from lielab.domain import CompositeManifold, CompositeGroup, CompositeAlgebra, RN, SO, SE, SU, rn, so, se, su
    from lielab.testing import check_topology
    
    # CompositeManifold
    assert check_topology(CompositeManifold([]), CompositeManifold([]))
    assert check_topology(CompositeManifold([RN(), SO(), SE(), su()]), CompositeManifold([RN(), SO(), SE(), su()]))
    assert check_topology(CompositeManifold([SO(5), SE(3), su(4)]), CompositeManifold([SO(5), SE(3), su(4)]))

    assert not check_topology(CompositeManifold([SO(3), SE(3), su(4)]), CompositeManifold([SO(5), SE(3), su(4)]))
    assert not check_topology(CompositeManifold([SO(5), SE(3), su(4)]), CompositeManifold([SO(5), SE(3), su(10)]))
    assert not check_topology(CompositeManifold([SE(3), SO(5), su(4)]), CompositeManifold([SO(5), SE(3), su(4)]))
    assert not check_topology(CompositeManifold([SO(5), SE(3), su(4)]), CompositeManifold([SO(5), SE(3), SU(4)]))
    assert not check_topology(CompositeManifold([SO(5), su(4)]), CompositeManifold([SO(5), SE(3), SU(4)]))
    assert not check_topology(CompositeManifold([SO(5), SE(3), su(4)]), CompositeManifold([SU(4)]))
    assert not check_topology(CompositeManifold([]), CompositeManifold([SO(5), SE(3), SU(4)]))
    

    # CompositeGroup
    assert check_topology(CompositeGroup([]), CompositeGroup([]))
    assert check_topology(CompositeGroup([RN(), SO(), SE()]), CompositeGroup([RN(), SO(), SE()]))
    assert check_topology(CompositeGroup([SO(5), SE(3)]), CompositeGroup([SO(5), SE(3)]))

    assert not check_topology(CompositeGroup([SO(3), SE(3)]), CompositeGroup([SO(5), SE(3)]))
    assert not check_topology(CompositeGroup([SO(5), SE(3)]), CompositeGroup([SO(5), SE(10)]))
    assert not check_topology(CompositeGroup([SE(3), SO(5)]), CompositeGroup([SO(5), SE(3)]))
    assert not check_topology(CompositeGroup([SO(5), SE(3)]), CompositeGroup([SO(5), SE(3), SU(4)]))
    assert not check_topology(CompositeGroup([SO(5)]), CompositeGroup([SO(5), SE(3), SU(4)]))
    assert not check_topology(CompositeGroup([SO(5), SE(3)]), CompositeGroup([SU(4)]))
    assert not check_topology(CompositeGroup([]), CompositeGroup([SO(5), SE(3), SU(4)]))

    # CompositeAlgebra
    assert check_topology(CompositeAlgebra([]), CompositeAlgebra([]))
    assert check_topology(CompositeAlgebra([rn(), so(), se()]), CompositeAlgebra([rn(), so(), se()]))
    assert check_topology(CompositeAlgebra([so(5), se(3)]), CompositeAlgebra([so(5), se(3)]))

    assert not check_topology(CompositeAlgebra([so(3), se(3)]), CompositeAlgebra([so(5), se(3)]))
    assert not check_topology(CompositeAlgebra([so(5), se(3)]), CompositeAlgebra([so(5), se(10)]))
    assert not check_topology(CompositeAlgebra([se(3), so(5)]), CompositeAlgebra([so(5), se(3)]))
    assert not check_topology(CompositeAlgebra([so(5), se(3)]), CompositeAlgebra([so(5), se(3), su(4)]))
    assert not check_topology(CompositeAlgebra([so(5)]), CompositeAlgebra([so(5), se(3), su(4)]))
    assert not check_topology(CompositeAlgebra([so(5), se(3)]), CompositeAlgebra([su(4)]))
    assert not check_topology(CompositeAlgebra([]), CompositeAlgebra([so(5), se(3), su(4)]))


def test_check_almost_equal_tol():
    from lielab.testing import check_almost_equal_tol
    import numpy as np

    nan = np.nan

    assert check_almost_equal_tol(1.0e10, 1.00001e10, 1.0e-8, 1.0e-5)
    assert not check_almost_equal_tol(1.0e-7, 1.0e-8, 1.0e-8, 1.0e-5)
    assert not check_almost_equal_tol(1.0e10, 1.0001e10, 1.0e-8, 1.0e-5)
    assert check_almost_equal_tol(1.0e-8, 1.0e-9, 1.0e-8, 1.0e-5)
    assert check_almost_equal_tol(1.0, 1.0, 1.0e-8, 1.0e-5)
    assert not check_almost_equal_tol(nan, nan, 1.0e-8, 1.0e-5)
    assert check_almost_equal_tol(1.0e-8, 0.0, 1.0e-8, 1.0e-5)
    assert not check_almost_equal_tol(1.0e-7, 0.0, 1.0e-8, 1.0e-5)
    assert not check_almost_equal_tol(1.0e-100, 0.0, 0.0, 1.0e-5)
    assert not check_almost_equal_tol(1.0e-7, 0.0, 0.0, 1.0e-5)
    assert check_almost_equal_tol(1.0e-10, 1.0e-20, 1.0e-8, 1.0e-5)
    assert check_almost_equal_tol(1.0e-10, 0.0, 1.0e-8, 1.0e-5)
    assert not check_almost_equal_tol(1.0e-10, 1.0e-20, 0.0, 1.0e-5)
    assert check_almost_equal_tol(1.0e-10, 0.999999e-10, 0.0, 1.0e-5)


def test_check_almost_equal_nulp():
    from lielab.testing import check_almost_equal_nulp
    import numpy as np

    inf = np.inf
    nan = np.nan

    # Check NaNs always fail
    assert not check_almost_equal_nulp(1.0, nan)
    assert not check_almost_equal_nulp(nan, 1.0)
    assert not check_almost_equal_nulp(nan, nan)

    # Check infs succeed when their direction is the same
    assert check_almost_equal_nulp(inf, inf)
    assert check_almost_equal_nulp(-inf, -inf)
    assert not check_almost_equal_nulp(inf, -inf)
    assert not check_almost_equal_nulp(-inf, inf)
    assert not check_almost_equal_nulp(1.0, inf)
    assert not check_almost_equal_nulp(inf, 1.0)
    
    # Check around 1.0
    assert check_almost_equal_nulp(1.0, 1.0)
    assert check_almost_equal_nulp(1.0, 1.0, 0)

    assert not check_almost_equal_nulp(1.0, np.nextafter(1.0, inf), 0)
    assert not check_almost_equal_nulp(1.0, np.nextafter(1.0, -inf), 0)

    assert check_almost_equal_nulp(1.0, np.nextafter(1.0, inf), 1)
    assert check_almost_equal_nulp(1.0, np.nextafter(1.0, -inf), 1)

    assert not check_almost_equal_nulp(1.0, np.nextafter(np.nextafter(1.0, inf), inf), 1)
    assert not check_almost_equal_nulp(1.0, np.nextafter(np.nextafter(1.0, -inf), -inf), 1)
    assert check_almost_equal_nulp(1.0, np.nextafter(np.nextafter(1.0, inf), -inf), 1)
    assert check_almost_equal_nulp(1.0, np.nextafter(np.nextafter(1.0, -inf), inf), 1)

    assert check_almost_equal_nulp(1.0, np.nextafter(np.nextafter(1.0, inf), inf), 2)
    assert check_almost_equal_nulp(1.0, np.nextafter(np.nextafter(1.0, -inf), -inf), 2)
    assert not check_almost_equal_nulp(1.0, np.nextafter(np.nextafter(np.nextafter(1.0, inf), inf), inf), 2)
    assert not check_almost_equal_nulp(1.0, np.nextafter(np.nextafter(np.nextafter(1.0, -inf), -inf), -inf), 2)

    # Check around 0
    assert check_almost_equal_nulp(0.0, 0.0)
    assert check_almost_equal_nulp(0.0, 0.0, 0)
    assert check_almost_equal_nulp(0.0, -0.0, 0)

    assert not check_almost_equal_nulp(0.0, np.nextafter(0.0, inf), 0)
    assert not check_almost_equal_nulp(0.0, np.nextafter(0.0, -inf), 0)

    assert check_almost_equal_nulp(0.0, np.nextafter(0.0, inf), 1)
    assert check_almost_equal_nulp(0.0, np.nextafter(0.0, -inf), 1)

    assert not check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, inf), inf), 1)
    assert not check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, -inf), -inf), 1)
    assert check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, inf), -inf), 1)
    assert check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, -inf), inf), 1)

    assert check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, inf), inf), 2)
    assert check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, -inf), -inf), 2)
    assert not check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(np.nextafter(0.0, inf), inf), inf), 2)
    assert not check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(np.nextafter(0.0, -inf), -inf), -inf), 2)

    # Check around a big number
    bignum = 1.6e5
    assert check_almost_equal_nulp(bignum, bignum)
    assert check_almost_equal_nulp(bignum, bignum, 0)

    assert not check_almost_equal_nulp(bignum, np.nextafter(bignum, inf), 0)
    assert not check_almost_equal_nulp(bignum, np.nextafter(bignum, -inf), 0)

    assert check_almost_equal_nulp(bignum, np.nextafter(bignum, inf), 1)
    assert check_almost_equal_nulp(bignum, np.nextafter(bignum, -inf), 1)

    assert not check_almost_equal_nulp(bignum, np.nextafter(np.nextafter(bignum, inf), inf), 1)
    assert not check_almost_equal_nulp(bignum, np.nextafter(np.nextafter(bignum, -inf), -inf), 1)
    assert check_almost_equal_nulp(bignum, np.nextafter(np.nextafter(bignum, inf), -inf), 1)
    assert check_almost_equal_nulp(bignum, np.nextafter(np.nextafter(bignum, -inf), inf), 1)

    assert check_almost_equal_nulp(bignum, np.nextafter(np.nextafter(bignum, inf), inf), 2)
    assert check_almost_equal_nulp(bignum, np.nextafter(np.nextafter(bignum, -inf), -inf), 2)
    assert not check_almost_equal_nulp(bignum, np.nextafter(np.nextafter(np.nextafter(bignum, inf), inf), inf), 2)
    assert not check_almost_equal_nulp(bignum, np.nextafter(np.nextafter(np.nextafter(bignum, -inf), -inf), -inf), 2)

    # Check gating near 0
    assert check_almost_equal_nulp(0.0, 0.0, 1, True)
    assert check_almost_equal_nulp(0.0, 0.0, 0, True)

    assert not check_almost_equal_nulp(0.0, np.nextafter(1.0, inf) - 1.0, 0, True)
    assert not check_almost_equal_nulp(0.0, np.nextafter(1.0, -inf) - 1.0, 0, True)

    ulp_smaller = np.abs(np.nextafter(1.0, -inf) - 1.0)
    assert check_almost_equal_nulp(0.0, (1.0 + ulp_smaller) - 1.0, 1, True)
    assert not check_almost_equal_nulp(0.0, np.nextafter(1.0, inf) - 1.0, 1, True)
    assert check_almost_equal_nulp(0.0, np.nextafter(1.0, inf) - 1.0, 2, True)
    assert check_almost_equal_nulp(0.0, np.nextafter(1.0, -inf) - 1.0, 1, True)

    assert not check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(1.0, inf), inf) - 1.0, 1, True)
    assert not check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(1.0, -inf), -inf) - 1.0, 1, True)
    assert check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, inf), -inf), 1, True)
    assert check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, -inf), inf), 1, True)

    assert check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, inf), inf), 2, True)
    assert check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(0.0, -inf), -inf), 2, True)
    assert not check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(np.nextafter(1.0, inf), inf), inf) - 1.0, 2, True)
    assert not check_almost_equal_nulp(0.0, np.nextafter(np.nextafter(np.nextafter(1.0, -inf), -inf), -inf) - 1.0, 2, True)
