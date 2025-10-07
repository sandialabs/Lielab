#include <catch2/catch_all.hpp>

#include <Lielab.hpp>
#include "../test_utils.hpp"

TEST_CASE("Killing", "[functions]")
{
    /*!
    * Tests the Killing function.
    */

    const Lielab::domain::so rx = Lielab::domain::so::from_vector({1.0, 0.0, 0.0});
    const Lielab::domain::so ry = Lielab::domain::so::from_vector({0.0, 1.0, 0.0});
    const Lielab::domain::so rz = Lielab::domain::so::from_vector({0.0, 0.0, 1.0});

    CHECK(std::abs(Lielab::functions::Killing(rx,rx) + 2.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rx,ry) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rx,rz) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(ry,rx) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(ry,ry) + 2.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(ry,rz) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rz,rx) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rz,ry) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rz,rz) + 2.0) <= TOL_FINE);
}

TEST_CASE("Killingform", "[functions]")
{
    /*!
    * Tests the killing form function.
    */

    const Lielab::domain::so rx = Lielab::domain::so::from_vector({1.0, 0.0, 0.0});

    Eigen::MatrixXd K = Lielab::functions::Killingform(rx);
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(rx.get_dimension(), rx.get_dimension());

    CHECK(std::abs(K.trace() + 6) <= TOL_FINE);
    assert_matrix(K*K.inverse(), Id);

    Lielab::domain::so so6 = Lielab::domain::so::basis(0,6);

    K = Lielab::functions::Killingform(so6);
    Id = Eigen::MatrixXd::Identity(so6.get_dimension(), so6.get_dimension());

    CHECK(std::abs(K.trace() + 120) <= TOL_FINE);
    assert_matrix(K*K.inverse(), Id);
}
