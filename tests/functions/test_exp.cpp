#include "../test_utils.hpp"

#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <cmath>
#include <iostream>
#include <numbers>

TEST_CASE("exp so2", "[functions]")
{
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    // Identity
    const so so2zero = so(2);
    const SO SO2Id = exp(so2zero);
    REQUIRE(SO2Id.get_shape() == 2);
    CHECK(SO2Id(0, 0) == 1.0);
    CHECK(SO2Id(0, 1) == 0.0);
    CHECK(SO2Id(1, 0) == 0.0);
    CHECK(SO2Id(1, 1) == 1.0);

    // Other
    const double pi = std::numbers::pi_v<double>;
    const SO SO2_mx = exp(pi/4.0*so::basis(0, 2));
    REQUIRE(SO2_mx.get_shape() == 2);
    CHECK_THAT(SO2_mx(0, 0), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-14));
    CHECK_THAT(SO2_mx(0, 1), Catch::Matchers::WithinAbs(-std::sqrt(2.0)/2.0, 1e-14));
    CHECK_THAT(SO2_mx(1, 0), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-14));
    CHECK_THAT(SO2_mx(1, 1), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-14));
}

TEST_CASE("exp so3", "[functions]")
{
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    // Identity and small valued tests. Tests the divide by zero is correctly handled.
    const double eps = 1.0e-15;

    const so so3_zero = so(3);
    const SO SO3_Id = exp(so3_zero);
    REQUIRE(SO3_Id.get_shape() == 3);
    CHECK(SO3_Id(0, 0) == 1.0);
    CHECK(SO3_Id(0, 1) == 0.0);
    CHECK(SO3_Id(0, 2) == 0.0);
    CHECK(SO3_Id(1, 0) == 0.0);
    CHECK(SO3_Id(1, 1) == 1.0);
    CHECK(SO3_Id(1, 2) == 0.0);
    CHECK(SO3_Id(2, 0) == 0.0);
    CHECK(SO3_Id(2, 1) == 0.0);
    CHECK(SO3_Id(2, 2) == 1.0);

    const so so3_epsx = eps*so::basis(0, 3);
    const SO SO3_epsx = exp(so3_epsx);
    REQUIRE(SO3_epsx.get_shape() == 3);
    CHECK(SO3_epsx(0, 0) == 1.0);
    CHECK(SO3_epsx(0, 1) == 0.0);
    CHECK(SO3_epsx(0, 2) == 0.0);
    CHECK(SO3_epsx(1, 0) == 0.0);
    CHECK(SO3_epsx(1, 1) == 1.0);
    CHECK(SO3_epsx(1, 2) == -eps);
    CHECK(SO3_epsx(2, 0) == 0.0);
    CHECK(SO3_epsx(2, 1) == eps);
    CHECK(SO3_epsx(2, 2) == 1.0);

    const so so3_epsy = eps*so::basis(1, 3);
    const SO SO3_epsy = exp(so3_epsy);
    REQUIRE(SO3_epsy.get_shape() == 3);
    CHECK(SO3_epsy(0, 0) == 1.0);
    CHECK(SO3_epsy(0, 1) == 0.0);
    CHECK(SO3_epsy(0, 2) == eps);
    CHECK(SO3_epsy(1, 0) == 0.0);
    CHECK(SO3_epsy(1, 1) == 1.0);
    CHECK(SO3_epsy(1, 2) == 0.0);
    CHECK(SO3_epsy(2, 0) == -eps);
    CHECK(SO3_epsy(2, 1) == 0.0);
    CHECK(SO3_epsy(2, 2) == 1.0);

    const so so3_epsz = eps*so::basis(2, 3);
    const SO SO3_epsz = exp(so3_epsz);
    REQUIRE(SO3_epsz.get_shape() == 3);
    CHECK(SO3_epsz(0, 0) == 1.0);
    CHECK(SO3_epsz(0, 1) == -eps);
    CHECK(SO3_epsz(0, 2) == 0.0);
    CHECK(SO3_epsz(1, 0) == eps);
    CHECK(SO3_epsz(1, 1) == 1.0);
    CHECK(SO3_epsz(1, 2) == 0.0);
    CHECK(SO3_epsz(2, 0) == 0.0);
    CHECK(SO3_epsz(2, 1) == 0.0);
    CHECK(SO3_epsz(2, 2) == 1.0);

    // Other
    const double pi = std::numbers::pi_v<double>;
    const double m = pi/4.0;

    const so so3_mx = m*so::basis(0, 3);
    const SO SO3_mx = exp(so3_mx);
    REQUIRE(SO3_mx.get_shape() == 3);
    CHECK(SO3_mx(0, 0) == 1.0);
    CHECK(SO3_mx(0, 1) == 0.0);
    CHECK(SO3_mx(0, 2) == 0.0);
    CHECK(SO3_mx(1, 0) == 0.0);
    CHECK_THAT(SO3_mx(1, 1), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-14));
    CHECK_THAT(SO3_mx(1, 2), Catch::Matchers::WithinAbs(-std::sqrt(2.0)/2.0, 1e-14));
    CHECK(SO3_mx(2, 0) == 0.0);
    CHECK_THAT(SO3_mx(2, 1), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-14));
    CHECK_THAT(SO3_mx(2, 2), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-14));
}
