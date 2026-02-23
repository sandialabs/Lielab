#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

TEST_CASE("Ad", "[functions]")
{
    /*!
    * Tests the Ad function.
    */

    using Lielab::domain::so;
    using Lielab::domain::SO;
    using Lielab::functions::Ad;
    using Lielab::functions::exp;
    using Lielab::testing::check_almost_equal_tol;

    const so u = so::from_vector({1.0, 0.0, 0.0});
    const so v = so::from_vector({0.0, 1.0, 0.0});
    const so w = so::from_vector({0.0, 0.0, 1.0});

    const SO Gso = exp(v);

    // GuG^-1
    so ansso = Ad(Gso, u);
    Eigen::MatrixXd truthso(3,3);
    truthso << 0, 0.841470984807896, 0,
              -0.841470984807897, 0, -0.540302305868140,
               0, 0.540302305868140, 0;
    
    CHECK(check_almost_equal_tol(ansso.get_matrix(), truthso));

    // GvG^-1 = v when G = exp(v)
    ansso = Ad(Gso, v);
    
    CHECK(check_almost_equal_tol(ansso.get_matrix(), v.get_matrix()));

    // GwG^-1
    ansso = Ad(Gso, w);
    truthso << 0, -0.540302305868140, 0,
               0.540302305868140, 0, -0.841470984807897,
               0, 0.841470984807897, 0;
    
    CHECK(check_almost_equal_tol(ansso.get_matrix(), truthso));
}
