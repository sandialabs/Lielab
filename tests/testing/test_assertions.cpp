#include <catch2/catch_all.hpp>
#include <Eigen/Core>

#include <Lielab/domain.hpp>
#include <Lielab/testing.hpp>

#include <cmath>
#include <limits>

TEST_CASE("check_topology", "[testing]")
{
    using namespace Lielab::domain;
    using Lielab::testing::check_topology;
    
    // CompositeManifold
    CHECK(check_topology(CompositeManifold{}, CompositeManifold{}));
    CHECK(check_topology(CompositeManifold{RN{}, SO{}, SE{}, su{}}, CompositeManifold{RN{}, SO{}, SE{}, su{}}));
    CHECK(check_topology(CompositeManifold{SO(5), SE(3), su(4)}, CompositeManifold{SO(5), SE(3), su(4)}));

    CHECK(!check_topology(CompositeManifold{SO(3), SE(3), su(4)}, CompositeManifold{SO(5), SE(3), su(4)}));
    CHECK(!check_topology(CompositeManifold{SO(5), SE(3), su(4)}, CompositeManifold{SO(5), SE(3), su(10)}));
    CHECK(!check_topology(CompositeManifold{SE(3), SO(5), su(4)}, CompositeManifold{SO(5), SE(3), su(4)}));
    CHECK(!check_topology(CompositeManifold{SO(5), SE(3), su(4)}, CompositeManifold{SO(5), SE(3), SU(4)}));
    CHECK(!check_topology(CompositeManifold{SO(5), su(4)}, CompositeManifold{SO(5), SE(3), SU(4)}));
    CHECK(!check_topology(CompositeManifold{SO(5), SE(3), su(4)}, CompositeManifold{SU(4)}));
    CHECK(!check_topology(CompositeManifold{}, CompositeManifold{SO(5), SE(3), SU(4)}));
    

    // CompositeGroup
    CHECK(check_topology(CompositeGroup{}, CompositeGroup{}));
    CHECK(check_topology(CompositeGroup{RN{}, SO{}, SE{}}, CompositeGroup{RN{}, SO{}, SE{}}));
    CHECK(check_topology(CompositeGroup{SO(5), SE(3)}, CompositeGroup{SO(5), SE(3)}));

    CHECK(!check_topology(CompositeGroup{SO(3), SE(3)}, CompositeGroup{SO(5), SE(3)}));
    CHECK(!check_topology(CompositeGroup{SO(5), SE(3)}, CompositeGroup{SO(5), SE(10)}));
    CHECK(!check_topology(CompositeGroup{SE(3), SO(5)}, CompositeGroup{SO(5), SE(3)}));
    CHECK(!check_topology(CompositeGroup{SO(5), SE(3)}, CompositeGroup{SO(5), SE(3), SU(4)}));
    CHECK(!check_topology(CompositeGroup{SO(5)}, CompositeGroup{SO(5), SE(3), SU(4)}));
    CHECK(!check_topology(CompositeGroup{SO(5), SE(3)}, CompositeGroup{SU(4)}));
    CHECK(!check_topology(CompositeGroup{}, CompositeGroup{SO(5), SE(3), SU(4)}));

    // CompositeAlgebra
    CHECK(check_topology(CompositeAlgebra{}, CompositeAlgebra{}));
    CHECK(check_topology(CompositeAlgebra{rn{}, so{}, se{}}, CompositeAlgebra{rn{}, so{}, se{}}));
    CHECK(check_topology(CompositeAlgebra{so(5), se(3)}, CompositeAlgebra{so(5), se(3)}));

    CHECK(!check_topology(CompositeAlgebra{so(3), se(3)}, CompositeAlgebra{so(5), se(3)}));
    CHECK(!check_topology(CompositeAlgebra{so(5), se(3)}, CompositeAlgebra{so(5), se(10)}));
    CHECK(!check_topology(CompositeAlgebra{se(3), so(5)}, CompositeAlgebra{so(5), se(3)}));
    CHECK(!check_topology(CompositeAlgebra{so(5), se(3)}, CompositeAlgebra{so(5), se(3), su(4)}));
    CHECK(!check_topology(CompositeAlgebra{so(5)}, CompositeAlgebra{so(5), se(3), su(4)}));
    CHECK(!check_topology(CompositeAlgebra{so(5), se(3)}, CompositeAlgebra{su(4)}));
    CHECK(!check_topology(CompositeAlgebra{}, CompositeAlgebra{so(5), se(3), su(4)}));
}

TEST_CASE("check_almost_equal_tol", "[testing]")
{
    using Lielab::testing::check_almost_equal_tol;

    const double nan = std::numeric_limits<double>::quiet_NaN();

    // Tests from numpy.isclose
    CHECK(check_almost_equal_tol(1.0e10, 1.00001e10, 1.0e-8, 1.0e-5));
    CHECK(!check_almost_equal_tol(1.0e-7, 1.0e-8, 1.0e-8, 1.0e-5));
    CHECK(!check_almost_equal_tol(1.0e10, 1.0001e10, 1.0e-8, 1.0e-5));
    CHECK(check_almost_equal_tol(1.0e-8, 1.0e-9, 1.0e-8, 1.0e-5));
    CHECK(check_almost_equal_tol(1.0, 1.0, 1.0e-8, 1.0e-5));
    CHECK(!check_almost_equal_tol(nan, nan, 1.0e-8, 1.0e-5));
    CHECK(check_almost_equal_tol(1.0e-8, 0.0, 1.0e-8, 1.0e-5));
    CHECK(!check_almost_equal_tol(1.0e-7, 0.0, 1.0e-8, 1.0e-5));
    CHECK(!check_almost_equal_tol(1.0e-100, 0.0, 0.0, 1.0e-5));
    CHECK(!check_almost_equal_tol(1.0e-7, 0.0, 0.0, 1.0e-5));
    CHECK(check_almost_equal_tol(1.0e-10, 1.0e-20, 1.0e-8, 1.0e-5));
    CHECK(check_almost_equal_tol(1.0e-10, 0.0, 1.0e-8, 1.0e-5));
    CHECK(!check_almost_equal_tol(1.0e-10, 1.0e-20, 0.0, 1.0e-5));
    CHECK(check_almost_equal_tol(1.0e-10, 0.999999e-10, 0.0, 1.0e-5));
}

TEST_CASE("check_almost_equal_nulp", "[testing]")
{
    using Lielab::testing::check_almost_equal_nulp;
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // Check NaNs always fail
    CHECK(!check_almost_equal_nulp(1.0, nan));
    CHECK(!check_almost_equal_nulp(nan, 1.0));
    CHECK(!check_almost_equal_nulp(nan, nan));

    // Check infs succeed when their direction is the same
    CHECK(check_almost_equal_nulp(inf, inf));
    CHECK(check_almost_equal_nulp(-inf, -inf));
    CHECK(!check_almost_equal_nulp(inf, -inf));
    CHECK(!check_almost_equal_nulp(-inf, inf));
    CHECK(!check_almost_equal_nulp(1.0, inf));
    CHECK(!check_almost_equal_nulp(inf, 1.0));
    
    // Check around 1.0
    CHECK(check_almost_equal_nulp(1.0, 1.0));
    CHECK(check_almost_equal_nulp(1.0, 1.0, 0));

    CHECK(!check_almost_equal_nulp(1.0, std::nextafter(1.0, inf), 0));
    CHECK(!check_almost_equal_nulp(1.0, std::nextafter(1.0, -inf), 0));

    CHECK(check_almost_equal_nulp(1.0, std::nextafter(1.0, inf), 1));
    CHECK(check_almost_equal_nulp(1.0, std::nextafter(1.0, -inf), 1));

    CHECK(!check_almost_equal_nulp(1.0, std::nextafter(std::nextafter(1.0, inf), inf), 1));
    CHECK(!check_almost_equal_nulp(1.0, std::nextafter(std::nextafter(1.0, -inf), -inf), 1));
    CHECK(check_almost_equal_nulp(1.0, std::nextafter(std::nextafter(1.0, inf), -inf), 1));
    CHECK(check_almost_equal_nulp(1.0, std::nextafter(std::nextafter(1.0, -inf), inf), 1));

    CHECK(check_almost_equal_nulp(1.0, std::nextafter(std::nextafter(1.0, inf), inf), 2));
    CHECK(check_almost_equal_nulp(1.0, std::nextafter(std::nextafter(1.0, -inf), -inf), 2));
    CHECK(!check_almost_equal_nulp(1.0, std::nextafter(std::nextafter(std::nextafter(1.0, inf), inf), inf), 2));
    CHECK(!check_almost_equal_nulp(1.0, std::nextafter(std::nextafter(std::nextafter(1.0, -inf), -inf), -inf), 2));

    // Check around 0
    CHECK(check_almost_equal_nulp(0.0, 0.0));
    CHECK(check_almost_equal_nulp(0.0, 0.0, 0));
    CHECK(check_almost_equal_nulp(0.0, -0.0, 0));

    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(0.0, inf), 0));
    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(0.0, -inf), 0));

    CHECK(check_almost_equal_nulp(0.0, std::nextafter(0.0, inf), 1));
    CHECK(check_almost_equal_nulp(0.0, std::nextafter(0.0, -inf), 1));

    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, inf), inf), 1));
    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, -inf), -inf), 1));
    CHECK(check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, inf), -inf), 1));
    CHECK(check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, -inf), inf), 1));

    CHECK(check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, inf), inf), 2));
    CHECK(check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, -inf), -inf), 2));
    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(std::nextafter(0.0, inf), inf), inf), 2));
    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(std::nextafter(0.0, -inf), -inf), -inf), 2));

    // Check around a big number
    const double bignum = 1.6e5;
    CHECK(check_almost_equal_nulp(bignum, bignum));
    CHECK(check_almost_equal_nulp(bignum, bignum, 0));

    CHECK(!check_almost_equal_nulp(bignum, std::nextafter(bignum, inf), 0));
    CHECK(!check_almost_equal_nulp(bignum, std::nextafter(bignum, -inf), 0));

    CHECK(check_almost_equal_nulp(bignum, std::nextafter(bignum, inf), 1));
    CHECK(check_almost_equal_nulp(bignum, std::nextafter(bignum, -inf), 1));

    CHECK(!check_almost_equal_nulp(bignum, std::nextafter(std::nextafter(bignum, inf), inf), 1));
    CHECK(!check_almost_equal_nulp(bignum, std::nextafter(std::nextafter(bignum, -inf), -inf), 1));
    CHECK(check_almost_equal_nulp(bignum, std::nextafter(std::nextafter(bignum, inf), -inf), 1));
    CHECK(check_almost_equal_nulp(bignum, std::nextafter(std::nextafter(bignum, -inf), inf), 1));

    CHECK(check_almost_equal_nulp(bignum, std::nextafter(std::nextafter(bignum, inf), inf), 2));
    CHECK(check_almost_equal_nulp(bignum, std::nextafter(std::nextafter(bignum, -inf), -inf), 2));
    CHECK(!check_almost_equal_nulp(bignum, std::nextafter(std::nextafter(std::nextafter(bignum, inf), inf), inf), 2));
    CHECK(!check_almost_equal_nulp(bignum, std::nextafter(std::nextafter(std::nextafter(bignum, -inf), -inf), -inf), 2));

    // Check gating near 0
    CHECK(check_almost_equal_nulp(0.0, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(0.0, 0.0, 0, true));

    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(1.0, inf) - 1.0, 0, true));
    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(1.0, -inf) - 1.0, 0, true));

    const double ulp_smaller = std::abs(std::nextafter(1.0, -inf) - 1.0);
    CHECK(check_almost_equal_nulp(0.0, (1.0 + ulp_smaller) - 1.0, 1, true));
    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(1.0, inf) - 1.0, 1, true));
    CHECK(check_almost_equal_nulp(0.0, std::nextafter(1.0, inf) - 1.0, 2, true));
    CHECK(check_almost_equal_nulp(0.0, std::nextafter(1.0, -inf) - 1.0, 1, true));

    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(1.0, inf), inf) - 1.0, 1, true));
    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(1.0, -inf), -inf) - 1.0, 1, true));
    CHECK(check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, inf), -inf), 1, true));
    CHECK(check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, -inf), inf), 1, true));

    CHECK(check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, inf), inf), 2, true));
    CHECK(check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(0.0, -inf), -inf), 2, true));
    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(std::nextafter(1.0, inf), inf), inf) - 1.0, 2, true));
    CHECK(!check_almost_equal_nulp(0.0, std::nextafter(std::nextafter(std::nextafter(1.0, -inf), -inf), -inf) - 1.0, 2, true));
}
