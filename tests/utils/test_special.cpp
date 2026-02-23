#include <catch2/catch_all.hpp>

#include <Lielab/utils/special.hpp>

#include <limits>

TEST_CASE("sign", "[utils]")
{
    using Lielab::utils::sign;

    const double inf = std::numeric_limits<double>::infinity();

    CHECK(sign(0.25) == sign(0.5));
    CHECK(sign(-0.25) == sign(-0.5));
    CHECK(sign(0.25) != sign(-0.5));
    CHECK(sign(-0.25) != sign(0.5));
    
    CHECK(sign(0.5) != sign(0.0));
    CHECK(sign(-0.5) != sign(0.0));
    CHECK(sign(0.0) == sign(0.0));

    CHECK(sign(inf) == sign(0.5));
    CHECK(sign(-inf) == sign(-0.5));
    CHECK(sign(inf) != sign(-0.5));
    CHECK(sign(-inf) != sign(0.5));
    CHECK(sign(inf) == sign(inf));
    CHECK(sign(-inf) == sign(-inf));
}
