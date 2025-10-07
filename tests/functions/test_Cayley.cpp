#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

TEST_CASE("Cayley", "[functions]")
{
    /*!
    * Tests the Cayley function
    */

    const Lielab::domain::so rx = Lielab::domain::so::from_vector({1.0, 0.0, 0.0});
    const Lielab::domain::so ry = Lielab::domain::so::from_vector({0.0, 1.0, 0.0});
    Eigen::MatrixXd ans(3, 3);

    // Values calculated by hand
    const Lielab::domain::SO ex1 = Lielab::functions::Cayley(rx);
    ans << 1.0, 0.0, 0.0,
           0.0, 0.6,-0.8,
           0.0, 0.8, 0.6;

    assert_matrix(ex1.get_matrix(), ans);

    const Lielab::domain::SO ex2 = Lielab::functions::Cayley(ry);
    ans << 0.6, 0.0, 0.8,
           0.0, 1.0, 0.0,
           -0.8, 0.0, 0.6;

    assert_matrix(ex2.get_matrix(), ans);
}

TEST_CASE("Cayley2", "[functions]")
{
    /*!
    * Tests the Cayley2 function.
    */

    const Lielab::domain::so rx = Lielab::domain::so::from_vector({1.0, 0.0, 0.0});
    const Lielab::domain::so ry = Lielab::domain::so::from_vector({0.0, 1.0, 0.0});
    Eigen::MatrixXd ans(3, 3);

    // Values calculated by hand
    const Lielab::domain::SO ex1 = Lielab::functions::Cayley2(rx);
    ans << 1.0, 0.0, 0.0,
           0.0, 0.6, -0.8,
           0.0, 0.8, 0.6;
    
    assert_matrix(ex1.get_matrix(), ans);

    const Lielab::domain::SO ex2 = Lielab::functions::Cayley2(rx + 2*ry);
    ans << 0.0, 0.0, 1.0,
           0.8, 0.6, 0.0,
           -0.6, 0.8, 0.0;
    
    assert_matrix(ex2.get_matrix(), ans);
}

TEST_CASE("cayley 1 and 2", "[functions]")
{
    /*!
    * Tests Cayley and Cayley2 together with known identities.
    */

    // Identity Cayley = Cayley2 for all basis elements
    const size_t dim = Lielab::domain::so::basis(0,10).get_dimension();
    for (size_t ii = 0; ii < dim; ii++)
    {
        const Lielab::domain::so g = Lielab::domain::so::basis(ii, 10);
        assert_domain(Lielab::functions::Cayley(g), Lielab::functions::Cayley2(g));
    }
}

TEST_CASE("dCayleyinv", "[functions]")
{
    /*!
    * Tests the inverse of the dCayley function
    */
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

    ansso = Lielab::functions::dCayleyinv(u, v);
    truthso << 0.0, 0.5, 1.0,
              -0.5, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    ansso = Lielab::functions::dCayleyinv(v, u);
    truthso << 0.0,-0.5, 0.0,
               0.5, 0.0,-1.0,
               0.0, 1.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);
}
