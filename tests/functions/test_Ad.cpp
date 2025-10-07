#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

TEST_CASE("Ad", "[functions]")
{
    /*!
    * Tests the Ad function.
    */
    Lielab::domain::so u(3);
    Lielab::domain::so v(3);
    Lielab::domain::so w(3);
    Lielab::domain::so ansso(3);
    Lielab::domain::SO Gso(3);

    Eigen::VectorXd xx(3);
    xx << 1.0, 0.0, 0.0;
    Eigen::VectorXd yy(3);
    yy << 0.0, 1.0, 0.0;
    Eigen::VectorXd zz(3);
    zz << 0.0, 0.0, 1.0;
    Eigen::MatrixXd truthso(3,3);

    u.set_vector(xx);
    v.set_vector(yy);
    w.set_vector(zz);
    Gso = Lielab::functions::exp(v);

    // GuG^-1
    ansso = Lielab::functions::Ad(Gso, u);
    truthso << 0, 0.841470984807896, 0,
              -0.841470984807897, 0, -0.540302305868140,
               0, 0.540302305868140, 0;
    
    assert_matrix(ansso.get_matrix(), truthso);

    // GvG^-1 = v when G = exp(v)
    ansso = Lielab::functions::Ad(Gso, v);
    
    assert_matrix(ansso.get_matrix(), v.get_matrix());

    // GwG^-1
    ansso = Lielab::functions::Ad(Gso, w);
    truthso << 0, -0.540302305868140, 0,
               0.540302305868140, 0, -0.841470984807897,
               0, 0.841470984807897, 0;
    
    assert_matrix(ansso.get_matrix(), truthso);
}
