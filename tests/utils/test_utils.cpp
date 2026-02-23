#include <catch2/catch_all.hpp>
#include <Eigen/Core>

#include <Lielab.hpp>
#include "../test_utils.hpp"

double myfun(double x)
{
    return x - 4.5;
}

// TEST_CASE("search_linearx", "[optim]")
// {
//     Lielab::optim::search_linearx search;
//     double x = 5.0;

//     search.lower = 4.0;
//     search.upper = 6.0;

//     search.init(x);

//     CHECK(std::abs(search.lower - 4.0) < TOL_FINE);
//     CHECK(std::abs(search.upper - 6.0) < TOL_FINE);

//     x = search.step(x, myfun(x));

//     CHECK(std::abs(x - 5.000006) < TOL_FINE);
//     CHECK(search.k == 3);

//     x = search.step(x, myfun(x));

//     CHECK(std::abs(x - 4.5) < TOL_FINE);
//     CHECK(search.k == 4);

//     x = search.step(x, myfun(x));

//     CHECK(std::abs(x - 4.5) < TOL_FINE);
//     CHECK(search.k == 8);

//     x = 5.0;

//     x = search(*myfun, x);

//     CHECK(std::abs(x - 4.5) < TOL_FINE);
//     CHECK(search.k == 8);

// }
