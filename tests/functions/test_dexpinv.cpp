#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

TEST_CASE("dexpinv_numerical", "[functions]")
{
    /*!
    * Tests the dexpinv_numerical function.
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


    // order = 0
    ansso = Lielab::functions::dexpinv_numerical(u, v, 0);
    truthso << 0.0, 0.0, 1.0,
               0.0, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 1
    ansso = Lielab::functions::dexpinv_numerical(u, v, 1);
    truthso << 0.0, 0.5, 1.0,
              -0.5, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 3
    ansso = Lielab::functions::dexpinv_numerical(u, v, 3);
    truthso << 0.0, 0.5, 0.916666666666667,
              -0.5, 0.0, 0.0,
              -0.916666666666667, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 11
    ansso = Lielab::functions::dexpinv_numerical(u, v, 11);
    truthso << 0.0, 0.500000000000000, 0.915243861398375,
              -0.500000000000000, 0.0, 0.0,
              -0.915243861398375, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    Lielab::domain::rn x(3);
    Lielab::domain::rn y(3);
    Lielab::domain::rn ansrn(3);
    Eigen::MatrixXd truthrn(4, 4);
    x.set_vector(xx);
    y.set_vector(yy);

    // default order
    ansrn = Lielab::functions::dexpinv_numerical(x, y);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_matrix(), truthrn);

    // ridiculous order (checks abelian speedhack)
    ansrn = Lielab::functions::dexpinv_numerical(x, y, 999999999);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_matrix(), truthrn);
}
