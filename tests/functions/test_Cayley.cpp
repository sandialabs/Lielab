#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

TEST_CASE("cay", "[functions]")
{
    /*!
    * Tests the cay function
    */

    using Lielab::testing::check_almost_equal_tol;

    const Lielab::domain::so rx = Lielab::domain::so::from_vector({1.0, 0.0, 0.0});
    const Lielab::domain::so ry = Lielab::domain::so::from_vector({0.0, 1.0, 0.0});
    Eigen::MatrixXd ans(3, 3);

    // Values calculated by hand
    const Lielab::domain::SO ex1 = Lielab::functions::cay(rx);
    ans << 1.0, 0.0, 0.0,
           0.0, 0.6,-0.8,
           0.0, 0.8, 0.6;

    CHECK(check_almost_equal_tol(ex1.get_matrix(), ans));

    const Lielab::domain::SO ex2 = Lielab::functions::cay(ry);
    ans << 0.6, 0.0, 0.8,
           0.0, 1.0, 0.0,
           -0.8, 0.0, 0.6;

    CHECK(check_almost_equal_tol(ex2.get_matrix(), ans));
}

TEST_CASE("cay2", "[functions]")
{
    /*!
    * Tests the cay2 function.
    */

    using Lielab::testing::check_almost_equal_tol;

    const Lielab::domain::so rx = Lielab::domain::so::from_vector({1.0, 0.0, 0.0});
    const Lielab::domain::so ry = Lielab::domain::so::from_vector({0.0, 1.0, 0.0});
    Eigen::MatrixXd ans(3, 3);

    // Values calculated by hand
    const Lielab::domain::SO ex1 = Lielab::functions::cay2(rx);
    ans << 1.0, 0.0, 0.0,
           0.0, 0.6, -0.8,
           0.0, 0.8, 0.6;
    
    CHECK(check_almost_equal_tol(ex1.get_matrix(), ans));

    const Lielab::domain::SO ex2 = Lielab::functions::cay2(rx + 2*ry);
    ans << 0.0, 0.0, 1.0,
           0.8, 0.6, 0.0,
           -0.6, 0.8, 0.0;
    
    CHECK(check_almost_equal_tol(ex2.get_matrix(), ans));
}

TEST_CASE("cay 1 and 2", "[functions]")
{
    /*!
    * Tests cay and cay2 together with known identities.
    */

    using namespace Lielab::domain;
    using Lielab::functions::cay;
    using Lielab::functions::cay2;
    using Lielab::testing::check_almost_equal_nulp;

    // Identity Cayley = Cayley2 for all basis elements
    const int dim = Lielab::domain::so::basis(0,10).get_dimension();
    for (int ii = 0; ii < dim; ii++)
    {
        const Lielab::domain::so g = Lielab::domain::so::basis(ii, 10);
        CHECK(check_almost_equal_nulp(CompositeGroup{cay(g)}, CompositeGroup{cay2(g)}, 1, true));
    }
}

TEST_CASE("dcayinv", "[functions]")
{
    /*!
    * Tests the inverse of the dcay function
    */

    using Lielab::testing::check_almost_equal_tol;

    Lielab::domain::so u(3);
    Lielab::domain::so v(3);
    Lielab::domain::so ansso(3);

    Eigen::VectorXd xx(3);
    xx << 1.0, 0.0, 0.0;
    Eigen::VectorXd yy(3);
    yy << 0.0, 1.0, 0.0;
    Eigen::MatrixXd truthso(3,3);

    u.set_vector(xx);
    v.set_vector(yy);

    ansso = Lielab::functions::dcayinv(u, v);
    truthso << 0.0, 0.5, 1.0,
              -0.5, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    CHECK(check_almost_equal_tol(ansso.get_matrix(), truthso));

    ansso = Lielab::functions::dcayinv(v, u);
    truthso << 0.0,-0.5, 0.0,
               0.5, 0.0,-1.0,
               0.0, 1.0, 0.0;

    CHECK(check_almost_equal_tol(ansso.get_matrix(), truthso));
}
