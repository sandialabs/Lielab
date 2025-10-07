#include "../test_utils.hpp"

#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <iostream>
#include <numbers>

TEST_CASE("log SO2", "[functions]")
{
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    // Identity and small valued tests.
    const double eps = 1.0e-15;

    const SO SO2Id = SO(2);
    const so zero = log(SO2Id);
    REQUIRE(zero.get_shape() == 2);
    const Eigen::VectorXd zerobar = zero.get_vector();
    REQUIRE(zerobar.size() == 1);
    CHECK(zerobar(0) == 0.0);

    const SO SO2_epsx = exp(eps*so::basis(0, 2));
    const so epsx = log(SO2_epsx);
    REQUIRE(epsx.get_shape() == 2);
    const Eigen::VectorXd epsxbar = epsx.get_vector();
    REQUIRE(epsxbar.size() == 1);
    CHECK(epsxbar(0) == eps);

    // Large-ish valued tests.
    const double m = 2.0;

    const SO SO2_mx = exp(m*so::basis(0, 2));
    const so mx = log(SO2_mx);
    REQUIRE(mx.get_shape() == 2);
    const Eigen::VectorXd mxbar = mx.get_vector();
    REQUIRE(mxbar.size() == 1);
    CHECK(mxbar(0) == m);
}

TEST_CASE("log SO3", "[functions]")
{
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    // Identity and small valued tests. Tests the divide by zero is correctly handled.
    const double eps = 1.0e-15;

    const SO SO3Id = SO(3);
    const so zero = log(SO3Id);
    REQUIRE(zero.get_shape() == 3);
    const Eigen::VectorXd zerobar = zero.get_vector();
    REQUIRE(zerobar.size() == 3);
    CHECK(zerobar(0) == 0.0);
    CHECK(zerobar(1) == 0.0);
    CHECK(zerobar(2) == 0.0);

    const SO SO3_epsx = exp(eps*so::basis(0, 3));
    const so epsx = log(SO3_epsx);
    REQUIRE(epsx.get_shape() == 3);
    const Eigen::VectorXd epsxbar = epsx.get_vector();
    REQUIRE(epsxbar.size() == 3);
    CHECK(epsxbar(0) == eps);
    CHECK(epsxbar(1) == 0.0);
    CHECK(epsxbar(2) == 0.0);

    const SO SO3_epsy = exp(eps*so::basis(1, 3));
    const so epsy = log(SO3_epsy);
    REQUIRE(epsy.get_shape() == 3);
    const Eigen::VectorXd epsybar = epsy.get_vector();
    REQUIRE(epsybar.size() == 3);
    CHECK(epsybar(0) == 0.0);
    CHECK(epsybar(1) == eps);
    CHECK(epsybar(2) == 0.0);

    const SO SO3_epsz = exp(eps*so::basis(2, 3));
    const so epsz = log(SO3_epsz);
    REQUIRE(epsz.get_shape() == 3);
    const Eigen::VectorXd epszbar = epsz.get_vector();
    REQUIRE(epszbar.size() == 3);
    CHECK(epszbar(0) == 0.0);
    CHECK(epszbar(1) == 0.0);
    CHECK(epszbar(2) == eps);

    const SO SO3_epsv = exp(eps*so::basis(0, 3) + 2*eps*so::basis(1, 3) - 1*eps*so::basis(2, 3));
    const so epsv = log(SO3_epsv);
    REQUIRE(epsv.get_shape() == 3);
    const Eigen::VectorXd epsvbar = epsv.get_vector();
    REQUIRE(epsvbar.size() == 3);
    CHECK(epsvbar(0) == eps);
    CHECK(epsvbar(1) == 2*eps);
    CHECK(epsvbar(2) == -eps);

    // Large-ish valued tests. Tests angle > pi/2 is correctly handled.
    const double m = 2.0;

    const SO SO3_mx = exp(m*so::basis(0, 3));
    const so mx = log(SO3_mx);
    REQUIRE(mx.get_shape() == 3);
    const Eigen::VectorXd mxbar = mx.get_vector();
    REQUIRE(mxbar.size() == 3);
    CHECK(mxbar(0) == m);
    CHECK(mxbar(1) == 0.0);
    CHECK(mxbar(2) == 0.0);

    const SO SO3_my = exp(m*so::basis(1, 3));
    const so my = log(SO3_my);
    REQUIRE(my.get_shape() == 3);
    const Eigen::VectorXd mybar = my.get_vector();
    REQUIRE(mybar.size() == 3);
    CHECK(mybar(0) == 0.0);
    CHECK(mybar(1) == m);
    CHECK(mybar(2) == 0.0);

    const SO SO3_mz = exp(m*so::basis(2, 3));
    const so mz = log(SO3_mz);
    REQUIRE(mz.get_shape() == 3);
    const Eigen::VectorXd mzbar = mz.get_vector();
    REQUIRE(mzbar.size() == 3);
    CHECK(mzbar(0) == 0.0);
    CHECK(mzbar(1) == 0.0);
    CHECK(mzbar(2) == m);

    const SO SO3_mv = exp(m/3.0*so::basis(0, 3) + 2*m/3.0*so::basis(1, 3) - 1*m/3.0*so::basis(2, 3));
    const so mv = log(SO3_mv);
    REQUIRE(mv.get_shape() == 3);
    const Eigen::VectorXd mvbar = mv.get_vector();
    REQUIRE(mvbar.size() == 3);
    CHECK_THAT(mvbar(0), Catch::Matchers::WithinAbs(m/3.0, 1e-14));
    CHECK_THAT(mvbar(1), Catch::Matchers::WithinAbs(2*m/3.0, 1e-14));
    CHECK_THAT(mvbar(2), Catch::Matchers::WithinAbs(-m/3.0, 1e-14));

    // Gimbal lock points
    const double pi = std::numbers::pi_v<double>;

    const SO SO3_pix = exp(pi*so::basis(0, 3));
    const so pix = log(SO3_pix);
    REQUIRE(pix.get_shape() == 3);
    const Eigen::VectorXd pixbar = pix.get_vector();
    REQUIRE(pixbar.size() == 3);
    CHECK_THAT(pixbar(0), Catch::Matchers::WithinAbs(pi, 1e-14));
    CHECK(pixbar(1) == 0.0);
    CHECK(pixbar(2) == 0.0);

    const SO SO3_piy = exp(pi*so::basis(1, 3));
    const so piy = log(SO3_piy);
    REQUIRE(piy.get_shape() == 3);
    const Eigen::VectorXd piybar = piy.get_vector();
    REQUIRE(piybar.size() == 3);
    CHECK(piybar(0) == 0.0);
    CHECK_THAT(piybar(1), Catch::Matchers::WithinAbs(pi, 1e-14));
    CHECK(piybar(2) == 0.0);

    const SO SO3_piz = exp(pi*so::basis(2, 3));
    const so piz = log(SO3_piz);
    REQUIRE(piz.get_shape() == 3);
    const Eigen::VectorXd pizbar = piz.get_vector();
    REQUIRE(pizbar.size() == 3);
    CHECK(pizbar(0) == 0.0);
    CHECK(pizbar(1) == 0.0);
    CHECK_THAT(pizbar(2), Catch::Matchers::WithinAbs(pi, 1e-14));
}
